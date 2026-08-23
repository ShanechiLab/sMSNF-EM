function mdlS = fitSwitchMultiscaleEM(trainSet,paramSet,nIter,pathModel)
% FITSWITCHMULTISCALEEM Fit system parameters using switch multiscale EM
%  Learn parameters for a switching dynamical system with multiscale
%  observations using a switch EM framework.
% 
%  INPUTS:
%  trainSet: struct containing training data and pertinent info
%     FIELDS (required):
%     fs:       sampling frequency (Hz), double
%     obsNt:    Poisson observations, formatted as single matrix
%               N(# of channels) x T(# of samples) or cell array of
%               separate trials with consistent N but allows for varying T.
%     obsYt:    Gaussian observations, formatted the same as obsNt.
%               Unobserved time points should be specified as column of
%               NaN.
%     FIELDS (optional):
%     stateXt:  behavior observations, formatted similar to obsNt as 
%               d(# of dimensions) x either T or T+1. T+1 for when x0
%               provided. Does not impact learned parameters. Only used to
%               track decoding performance of learned parameters through a
%               similarity transform from decoded xt to stateXt.
%     stateSt:  regime observations formatted the same as obsNt with
%               columns representing either 1D true regime labels within
%               [1,M] or equivalent MD one hot represenations. Does not
%               impact learned parameters. Only used to track decoding
%               performance of learned parameters.
%  paramSet: struct containing learn settings
%     FIELDS (main):
%     dimXtEst: d, dimension of latent brain state of model
%     dimStEst: M, # of regimes of model
%     obsPred:  bool for whether to track observation prediction metrics
%     statObs:  bool for fixed observation parameters across regimes
%     statDyn:  bool for fixed dynamic parameters across regimes
%     verbose:  bool for printing tracked metrics
%     filtSelect: str for choosing stationary filter to embed in switch
%               filter, can choose 'MSF','MSNF' see EM helpers/FuncOrg.m'
%
%  nIter: # of EM iterations to run
%  pathModel: path to where learned model is saved. If provided, is also 
%             where progress is saved to stop/resume learning. Increasing
%             the # of iterations can also pick up from last saved
%             iteration
% 
%  OUTPUT:
%  mdlS: struct containing learned model parameters and any tracked metrics
%     FIELDS (guaranteed):
%     thetaCell:cell array 1 x nIter+1 of learned parameters from each EM
%               iteration with thetaCell{1} being the initial parameters
%               and thetaCell{end} being the final learned parameters
%     static:   struct containing learning parameters and misc info if EM
%               iterations unexpectedly fail
%     FIELDS (metrics):
%     ntPrdMet: PP = 2*AUC-1, One step ahead prediction of spikes nt using
%               P(nt|nt-1)
%     xtCC:     Correlation coefficient between provided stateXt and
%               decoded xt post similarity transform
%     stAcc:    Accuracy of decoded st against provided stateSt
%     ytPrdMet: NRMSE of N step ahead prediction of yt where N is the
%               period between subsequent yt
% 
%  [1] Kim, Unsupervised learning of multiscale switching dynamical system
%  models from multimodal neural data
% 
%  Author: DongKyu Kim, Apr 2026, dongkyuk(at)usc(dot)edu

    %% Setting setup
    fs = trainSet.fs;
    xtGiven = isfield(trainSet,'stateXt');
    stGiven = isfield(trainSet,'stateSt');
    if ~exist('paramSet','var') || isempty(paramSet)
        paramSet = struct;
    end
    if ~isfield(paramSet,'verbose')
        verbose = false;
        curFile = '';
    else
        verbose = paramSet.verbose;
        curFile = 'fitSwitchMultiscaleEM';
    end
    if isfield(paramSet,'fsYt') && isempty(paramSet.fsYt)
        paramSet.fsYt = fs;
    end
    varList = {'obsPred','fsYt','verbose','saveFreq'};
    defaultVals = {true,fs,false,1};
    paramSet = insertStructDefaults(paramSet,varList,defaultVals,curFile);

    if isfield(paramSet,'printName')
        msg = sprintf(', %s',paramSet.printName);
    else
        msg = '';
    end
    if verbose
        fprintf('%s%s',curFile,msg);
    end
    
    %% Metrics to track

    genFields = {'runTime'};
    genFSizes = 3;
    obsFields = {'ntPrdMet','ytPrdMet'};
    obsFSizes = [1,1];
    xFields = {'xtCC'};
    xFSizes = 1;
    sFields = {'stAcc'};
    sFSizes = 1;
    
    mdlFields = genFields;
    mdlFSizes = genFSizes;
    
    if paramSet.obsPred
        mdlFields = [mdlFields, obsFields];
        mdlFSizes = [mdlFSizes, obsFSizes];
    end
    if xtGiven
        mdlFields = [mdlFields, xFields];
        mdlFSizes = [mdlFSizes, xFSizes];
    end
    if stGiven
        mdlFields = [mdlFields, sFields];
        mdlFSizes = [mdlFSizes, sFSizes];
    end

    %% Checking for and Loading Existing Model
    if exist('pathModel','var') && ~isempty(pathModel)
        savePerIter = true;
        dirLearn = fileparts(pathModel);
        if ~exist(dirLearn,'dir')
            mkdir(dirLearn);
        end
    else
        savePerIter = false;
    end
    if savePerIter && exist(pathModel,'file')
        [mdlS,newMdl,toRun,startIter,newIter] = loadModel(pathModel,...
                                                  paramSet,nIter);
        if newIter ~= nIter
            nIter = newIter;
        end
    else
        newMdl = true;
        toRun = true;
        startIter = 1;
    end

    if toRun == false
        if verbose
            fprintf(' Existing completed detected,');
        end
    else
        %% Start prep for EM             
        if verbose
            fprintf('Start: %s',datetime);
        end
        
        if newMdl == true
            %% Initialization
            thetaInit = getInitGeneral(trainSet,paramSet);

            mdlS = struct;
            mdlS.thetaCell = cell(1,nIter+1);
            mdlS.static = struct;
            for fInd = 1:length(mdlFields)
                field = mdlFields{fInd};
                mdlS.(field) = zeros(mdlFSizes(fInd),nIter);
            end
            mdlS.thetaCell{1} = thetaInit;
            mdlS.static.paramSet = paramSet;
        end
        %% Prepping observations and objects
        [obsYtAll,obsNtAll,~,~,~] = prepInputObs(trainSet,paramSet);
        nSegs = length(obsYtAll);
        if paramSet.obsPred
            ytTot = catCell(obsYtAll);
            ntTot = catCell(obsNtAll);
        end
        if xtGiven == true
            xtWhole = prepXtCat(trainSet);
        else
            xtWhole = [];
        end
        if stGiven == true
            stWhole = catCell(trainSet.stateSt);
        end
        [sumClct,memMngr,funMngr] = initObjects(obsYtAll,...
                                                mdlS.thetaCell{startIter},...
                                                paramSet.obsPred,...
                                                xtGiven,stGiven,paramSet);

        %% EM
        bailedOut = false;
        iter = startIter;
        if verbose
            fprintf(', iter:');
        end
        while iter <= nIter && ~bailedOut
            if verbose
                fprintf('%d,',iter);
            end
            tic;
            thetaEst = mdlS.thetaCell{iter};
            
            funMngr.resetTheta(thetaEst);
            memMngr.updateTheta(thetaEst);
            sumClct.reset();
            %% E-Step
            for s = 1:nSegs
                obsYtSeg = obsYtAll{s}; obsNtSeg = obsNtAll{s};
                segBail = estepSEM(thetaEst,obsYtSeg,obsNtSeg,...
                                      memMngr,sumClct,funMngr);
                bailedOut = bailedOut || segBail;
            end
            mdlS.runTime(1,iter) = toc; tic;

            if bailedOut
                mdlS.static.bailedOut = true;
                mdlS.static.bailedOutInd = iter;
                mdlS.static.bailedOutTheta = thetaEst;
                iter = iter - 1;
            else
                %% M-Step
                thetaNew = mstepEM(sumClct,paramSet,thetaEst);
                mdlS.runTime(2,iter) = toc; tic;
                
                %% metrics
                if verbose
                    fprintf('(');
                end
                if paramSet.obsPred
                    inds = isfinite(sumClct.yPrdTot(1,:));                        
                    yNRMSE = averageCorrCoef(ytTot(:,inds),...
                            sumClct.yPrdTot(:,inds));
                    mdlS.ytPrdMet(:,iter) = yNRMSE;
                    if verbose
                        fprintf('y%0.3g',yNRMSE);
                    end
                    spPP = getPP(sumClct.PNtPrdTot,ntTot);
                    mdlS.ntPrdMet(:,iter) = spPP;
                    if verbose
                        fprintf('n%0.3g',spPP);
                    end
                end
                if xtGiven == true
                    indsTemp = false(size(sumClct.xDecTot(1,:)));
                    indsTemp(1:paramSet.Y_masking_rate:end) = true;
                    xTApp = [ones(1, size(sumClct.xDecTot(:,indsTemp), 2)); sumClct.xDecTot(:,indsTemp)];
                    simTran = learnProjSignal(xtWhole(:,indsTemp), xTApp);
                    mdlS.thetaCell{iter}.simTran = simTran;
                    xCC = averageCorrCoef(xtWhole(:,indsTemp),simTran*xTApp);
                    mdlS.xtCC(iter) = xCC;
                    if verbose
                        fprintf('x%0.3g',xCC);
                    end
                end
                if stGiven == true
                    stAcc = max(calcStAccuracy(sumClct.PStDecTot,stWhole));
                    mdlS.stAcc(iter) = stAcc;
                    if verbose
                        fprintf('s%0.3g',stAcc);
                    end
                end
                if verbose
                    fprintf('pSty%0.3g', mean(diag(thetaEst.sTran)));
                end

                mdlS.runTime(3,iter) = toc; tic;
                mdlS.thetaCell{iter + 1} = thetaNew;
                if verbose
                    fprintf('),');
                end
            end
            iter = iter + 1;
            if savePerIter
                saveCond = paramSet.saveFreq <= 1 ...
                           || mod(iter-1,paramSet.saveFreq)==1 ...
                           || ~(iter <= nIter);
                if saveCond
                    save(pathModel,'-struct','mdlS','-v7.3');
                end
            end
        end
        
        if bailedOut
            mdlS = trimModel(mdlS,iter - 1);
            if savePerIter
                save(pathModel,'-struct','mdlS','-v7.3');
            end
            if verbose
                fprintf(', bailed out. ');
            end
        end
    end
    if verbose
        fprintf(' End: %s \n', datetime);
    end
end


function [sumO,memO,funO] = initObjects(obsCell,theta,...
                                        obsPred,xtGiven,stGiven,lrnS)
    dimXt = size(theta.Acell{1},1);
    dimSt = length(theta.sInit);
    dimYt = size(theta.Ccell{1},1);
    dimNt = length(theta.alphaCell{1});
    nSegs = length(obsCell);
    segLens = zeros(1,nSegs);
    for s = 1:nSegs
        segLens(s) = size(obsCell{s},2);
    end
    sviseFlag = false;
    funO = FuncOrg(theta,lrnS);
    sumO = SumSEM(dimXt,dimSt,dimYt,dimNt,...
                  segLens,...
                  obsPred,xtGiven,stGiven,sviseFlag);
    memO = MemSEM(theta,dimXt,dimSt,dimYt,dimNt,max(segLens),lrnS);
end
