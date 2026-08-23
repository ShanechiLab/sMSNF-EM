function thetaInit = genInitThetaData(paramSet,trainSet)
% GENINITTHETADATA Parameter initializaion for EM based on data
%  Randomly generates initial system parameters based on provided data and
%  initialization settings.
%  INPUTS:
%  paramSet: struct containing initialization / loading settings
%     FIELDS (required):
%     dimXtEst:     d, dimension of latent brain state of model
%     dimStEst:     M, # of regimes of model
%     pStay:        diagonal values of regime transition matrix
%     sTran:        prespecified regime transition matrix (alternative to
%                   using dimStEst and pStay)
%     sTranEst:     same as sTran
%     FIELDS (optional):
%     Qparams:      1x2 double of mean and variance of Gaussian that Q's
%                   eigenvalues are drawn from
%     initDiagQ:    bool for initializing Q as diagonal matrices
%     statObs:      bool for fixed observation parameters across regimes
%     AdiagVal:     double for what to initialize diagonal values of A as,
%                   defaulted to 1 (identity matrices)
%
%  Author: DongKyu Kim, Apr 2026, dongkyuk(at)usc(dot)edu

    if isfield(paramSet,'dimXtEst')
        dimXt = paramSet.dimXtEst;
    else
        error('Need dimXtEst field in paramSet');
    end
    fs = trainSet.fs;
    if isfield(paramSet,'sTran')
        dimSt = size(paramSet.sTran,1);
        sTran = paramSet.sTran;
    elseif isfield(paramSet,'dimStEst') && isfield(paramSet,'pStay')
        dimSt = paramSet.dimStEst;
        sTran = generateMarkovParams(paramSet.pStay,dimSt);
    elseif isfield(paramSet,'sTranEst')
        dimSt = size(paramSet.sTranEst,1);
        sTran = paramSet.sTranEst;
    else
        error('Need dimStEst info');
    end
    
    fields = {'Qparams','initDiagQ','initDiagR',...
			  'statObs','AdiagVal'};
	defaults = {[0.025,0.002],true,true,...
				false,1};
    paramSet = insertStructDefaults(paramSet,fields,defaults);
    
    %% Prepping data

    [obsYtAll,obsNtAll,~,~,~] = prepInputObs(trainSet,paramSet);
    ytTot = catCell(obsYtAll);
    ntTot = catCell(obsNtAll);
    tTot = size(ntTot,2);

    dimYt = size(obsYtAll{1},1);
    dimNt = size(obsNtAll{1},1);
    ntSum = sum(ntTot,2);
    
    thetaInit.Acell = cell(dimSt,1);
    thetaInit.Ccell = cell(dimSt,1);
    thetaInit.Qcell = cell(dimSt,1);
    thetaInit.Rcell = cell(dimSt,1);
    thetaInit.alphaCell = cell(dimSt,1);
    thetaInit.betaCell = cell(dimSt,1);
    thetaInit.x0Mean = zeros(dimXt,1);
    thetaInit.x0Cov = eye(dimXt);
    thetaInit.sInit = (1/dimSt).*ones(dimSt,1);
    thetaInit.sTran = (1./sum(sTran,1)).*sTran;

    %% Generating A Matrices
    thetaInit.Acell(:) = {paramSet.AdiagVal.*eye(dimXt)};

    
    %% Generating xt Noise Matrix
    Qprms = paramSet.Qparams;
    for m = 1:dimSt
        thetaInit.Qcell{m} = getRandCov(dimXt,Qprms(1),Qprms(2),...
                                        'norm',paramSet.initDiagQ);
    end
    QvarSum = zeros(dimXt,1);
    for m = 1:dimSt
        QvarSum = QvarSum + diag(thetaInit.Qcell{m});
    end
    Qvar = 1/dimSt .* QvarSum;

    %% Generating C and R Matrices
        ytFinInds = ~isnan(ytTot(1,:));
        ytFinTot = ytTot(:,ytFinInds);

        ytVar = var(ytFinTot,[],2);
        signalVar = (1/8).*ytVar;
        noiseVar = (7/8).*ytVar;
        
        Xvar = 100*Qvar;
        if paramSet.statObs == true
            Cbase = getRandCInit(dimYt, dimXt);
            Cscale = sqrt( signalVar ./ (Cbase.^2 * Xvar));
            thetaInit.Ccell(:) = {Cscale.*Cbase};
        else
            for m = 1:dimSt
                Cbase = getRandCInit(dimYt, dimXt);
                Cscale = sqrt( signalVar ./ (Cbase.^2 * Xvar));
                thetaInit.Ccell{m} = Cscale.*Cbase;
            end
        end
        thetaInit.Rcell(:) = {diag(noiseVar)};

    %% Generating Poisson Parameters
        baseFires = ntSum ./ (tTot/fs); %avg firerate
        
        baseL = 0.8.*(0.5.*baseFires);
        baseH = 1.2.*(0.5.*baseFires);
        thetaInit.alphaCell(:) = {zeros(1,dimNt)};
        for c = 1:dimNt
            bases = getUnifRand(baseL(c),baseH(c),1,dimSt);
            for m = 1:dimSt
                a = log(bases(m)) - log(fs);
                thetaInit.alphaCell{m}(c) = a;
            end
        end
        
        xMax = 4*sqrt(1000*max(Qvar)); % V[xt] ~ 1/(1-decay)^2 * V[wt]
        betaRngs = abs(log(60./min(baseFires,59)))./(dimXt*xMax);
        
        if paramSet.statObs
            alphaTemp = thetaInit.alphaCell{1};
            thetaInit.alphaCell(:) = {alphaTemp};
            
            betaTemp = zeros(dimXt,dimNt);
            for c = 1:dimNt
                bL = -1*betaRngs(c);
                bH = betaRngs(c);
                betaTemp(:,c) = getUnifRand(bL,bH,dimXt,1);
            end
            thetaInit.betaCell(:) = {betaTemp};
        else
            for m = 1:dimSt
                betaTemp = zeros(dimXt,dimNt);
                for c = 1:dimNt
                    bL = -1*betaRngs(c);
                    bH = betaRngs(c);
                    betaTemp(:,c) = getUnifRand(bL,bH,dimXt,1);
                end
                thetaInit.betaCell{m} = betaTemp;
            end
        end
end