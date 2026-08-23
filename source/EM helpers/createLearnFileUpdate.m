function saveS = createLearnFileUpdate(pathLearn,loadS)
    saveS = struct;
    
    %% Named Parameters
    saveS.dimXtEst = 6; saveS.dimStEst = 1;
    saveS.fsYt = 50;
    saveS.ytStartInd = 1;
    saveS.statObs = false; saveS.statDyn = false;
    saveS.mDiagQ = 'none'; saveS.mDiagR = 'none';
    saveS.mWeight = 1;

    %% Non-named parameters
    saveS.nIter = 150;
    saveS.initGiven = []; %[],true,paramSet,saved
    saveS.initLocation = '';
    saveS.pStay = 0.95; saveS.AdiagVal = 0.95;
    saveS.Qparams = [0.025,0.002]; % normal~[mean,var]
    saveS.initDiagQ = false; saveS.initDiagR = true;
    saveS.obsPred = true;
    saveS.saveFreq = 1; saveS.verbose = true;
    saveS.nPts = 5; saveS.nSamps = 1;
    if exist('loadS','var') && isstruct(loadS)
        owFields = fieldnames(loadS);
        for fInd = 1:length(owFields)
            field = owFields{fInd};
            if isfield(saveS,field)
                saveS.(field) = loadS.(field);
            end
        end
    end
    if exist('pathLearn','var') && ~isempty(pathLearn)
        if ~exist(fileparts(pathLearn),'dir')
            mkdir(fileparts(pathLearn));
        end
        save(pathLearn,'-struct','saveS');
    end
end

