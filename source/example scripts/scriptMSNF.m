%% Setting for learning
    paramSet1 = struct;
    paramSet1.dimXtEst = 10;
    paramSet1.dimStEst = 2;
    paramSet1.statObs = true;
    paramSet1.statDyn = false; 
    paramSet1.verbose = true;
    paramSet1.filtSelect = 'MSNF';
    paramSet1.fs = 100;
    paramSet1.Y_masking_rate = 5; %masking rate for LFP
    paramSet1.fsYt = paramSet1.fs/paramSet1.Y_masking_rate;
    paramSet1.ytStartInd = 1;
    paramSet1.mDiagQ = 'none';
    paramSet1.mDiagR = 'all';
    paramSet1.mWeight = 1;
    paramSet1.nIter = 100;
    paramSet1.initGiven = []; paramSet1.initLocation = '';
    paramSet1.pStay = 0.995; paramSet1.AdiagVal = 0.90;
    paramSet1.Qparams = [0.025, 0.002];
    paramSet1.initDiagQ = false; paramSet1.initDiagR = true;
    paramSet1.saveFreq = 5;
    paramSet1.nPts = 5;
    paramSet1.nSamps = 1;
    paramSet1.sf = 1; %scaling factor for LFP
%% MSNF-EM
trainSet = load("./experimental/train.mat"); 
init;
NUM_ITER = 300;
mdlEM = fitSwitchMultiscaleEM(trainSet, paramSet1, NUM_ITER);%sMSNF-EM
 
%% MSF-EM
paramSet2 = paramSet1;
paramSet2.filtSelect = 'MSF';
mdlEM2 = fitSwitchMultiscaleEM(trainSet_small, paramSet2, NUM_ITER);%sMSF-EM
