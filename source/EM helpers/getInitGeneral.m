function thetaInit = getInitGeneral(trainSet,paramSet)
% GETINITGENERAL Parameter initializaion for EM
%  General function for getting an intial set of model parameters (theta) 
%  for a stationary/switching dynamical system. Will try to load a provided
%  theta, otherwise will randomly generate using genInitThetaData.
%  INPUTS:
%  trainSet: struct containing training data and pertinent info
%     FIELDS: see fitSwitchMultiscalEM
%  paramSet: struct containing initialization / loading settings
%     FIELDS (main):
%     initGiven:    str for how to try loading initial theta. Choose
%                   between [] (default) for no loading, 'paramSet' for
%                   checking in paramSet, or 'saved' for looking in
%                   paramSet.initLocation
%     initLocation: path for saved initial theta if initGiven set to
%                   'saved'
%  OUTPUT:
%  thetaNew: struct containing initial system parameters
%
%  Author: DongKyu Kim, Apr 2026, dongkyuk(at)usc(dot)edu

    fields = {'initGiven','initLocation'};
    default = {[],[]};
    paramSet = insertStructDefaults(paramSet,fields,default);
    
    needInit = true;
    %% Try to load in a valid theta
    if ~isempty(paramSet.initGiven)
        loadS = [];
        switch lower(paramSet.initGiven)
            case 'true'
                if isfield(trainSet,'thetaTrue')
                    loadS = trainSet.thetaTrue;
                elseif isfield(paramSet,'thetaTrue')
                    loadS = paramSet.thetaTrue;
                end
            case 'paramset'
                if isfield(paramSet,'thetaInit')
                    loadS = paramSet.thetaInit;
                elseif isfield(paramSet,'theta')
                    loadS = paramSet.theta;
                end
            case 'saved'
                pathInit = paramSet.initLocation;
                if isfile(pathInit)
                    loadS = load(pathInit);
                end
        end
        
        if iscell(loadS)
            loadS = loadS{1};
        end
        if isfield(loadS,'Acell')
            thetaInit = loadS;
            needInit = false;
        elseif isfield(loadS,'thetaInit')
            if isfield(loadS.thetaInit,'Acell')
                thetaInit = loadS.thetaInit;
                needInit = false;
            end
        elseif isfield(loadS,'theta')
            if isfield(loadS.theta,'Acell')
                thetaInit = loadS.theta;
                needInit = false;
            end
        elseif isfield(loadS,'thetaCell')
            if isfield(loadS.thetaCell{1},'Acell')
                thetaInit = loadS.thetaCell{1};
                needInit = false;
            end
        end
    end
    
    %% If failed to load or init randomly
    if needInit
        baseS = paramSet;
        baseS.initGiven = [];
        baseS.verbose = false;

        thetaInit = genInitThetaData(baseS,trainSet);
    end

end