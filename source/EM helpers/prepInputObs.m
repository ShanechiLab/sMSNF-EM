function [obsYtAll,obsNtAll,ytFinInds,ytOut,ntOut,dimYt,dimNt] = prepInputObs(trainSet,paramSet,noCell)
    if ~exist('paramSet','var')
        paramSet = struct;
    end
    if ~exist('noCell','var')
        noCell = false;
    end
    
    varList = {'fsYt','ytStartInd'};
    defaultVals = {trainSet.fs,1};    

    paramSet = insertStructDefaults(paramSet,varList,defaultVals);
    
    if isempty(paramSet.fsYt)
        paramSet.fsYt = trainSet.fs;
    end
    
    obsHandle = [];
    ytIsCell = false;
    ntIsCell = false;
    if isfield(trainSet,'obsYt') && ~isempty(trainSet.obsYt)
        if iscell(trainSet.obsYt)
            finiteFound = false;
            ind = 1;
            while ind <= length(trainSet.obsYt) && ~finiteFound
                if any(isfinite(trainSet.obsYt{ind}(:)))
                    finiteFound = true;
                end
                ind = ind + 1;
            end
        else
            finiteFound = any(isfinite(trainSet.obsYt(:)));
        end
        if finiteFound == true
            ytFound = true;
            if iscell(trainSet.obsYt)
                ytIsCell = true;
            end
            obsHandle = trainSet.obsYt;
        else
            ytFound = false;
        end
    else
        ytFound = false;
    end
    if isfield(trainSet,'obsNt') && ~isempty(trainSet.obsNt)
        if iscell(trainSet.obsNt)
            finiteFound = false;
            ind = 1;
            while ind <= length(trainSet.obsNt) && ~finiteFound
                if any(isfinite(trainSet.obsNt{ind}(:)))
                    finiteFound = true;
                end
                ind = ind + 1;
            end
        else
            finiteFound = any(isfinite(trainSet.obsNt(:)));
        end
        if finiteFound == true
            ntFound = true;
            if iscell(trainSet.obsNt)
                ntIsCell = true;
            end
            if ytFound && (ytIsCell ~= ntIsCell)
                error(['Multiscale observation formatting '...
                       'must be consistent']);
            end
            if isempty(obsHandle)
                obsHandle = trainSet.obsNt;
            end
        else
            ntFound = false;
        end
    else
        ntFound = false;
    end
    if (~ytFound) && (~ntFound)
        error(['training struct must contain obsYt and/or '...
               'obsNt']);
    end
    if ~ytFound
        if ntIsCell
            nSegs = length(obsHandle);
            trainSet.obsYt = cell(1,nSegs);
            for k = 1:nSegs
                trainSet.obsYt{k} = nan(1,size(obsHandle{k},2));
            end
        else
            trainSet.obsYt = nan(1,size(obsHandle,2));
        end
    end
    if ~ntFound
        if ytIsCell
            nSegs = length(obsHandle);
            trainSet.obsNt = cell(1,nSegs);
            for k = 1:nSegs
                trainSet.obsNt{k} = nan(1,size(obsHandle{k},2));
            end
        else
            trainSet.obsNt = nan(1,size(obsHandle,2));
        end
    end
    obsCell = ytIsCell || ntIsCell;
    if obsCell == false
        ytOrig = {trainSet.obsYt};
        ntOrig = {trainSet.obsNt};

        dimYt = size(ytOrig{1},1);
        dimNt = size(ntOrig{1},1);
    else
        ytOrig = trainSet.obsYt;
        ntOrig = trainSet.obsNt;
        if length(ytOrig) ~= length(ntOrig)
            error(['Segmented observations must have'...
                   ' the same number of segments']);
        end
        dimYt = size(ytOrig{1},1);
        dimNt = size(ntOrig{1},1);
    end
    
    nSegs = length(ytOrig);
    obsYtAll = cell(1,nSegs);
    ytFinInds = cell(1,nSegs);
    obsNtAll = cell(1,nSegs);
    ytOut = ytFound;
    ntOut = ntFound;
    for sInd = 1:nSegs
        if ytOut
            obsYt = ytOrig{sInd};
            if any(isnan(obsYt(:)))
                obsYtAll{sInd} = obsYt;
                ytFinInds{sInd} = find(sum(isfinite(obsYt)));
            else
                ds = trainSet.fs/paramSet.fsYt;
                startInd = paramSet.ytStartInd;
                [ytDS,finInds] = downsampCont(obsYt,ds,startInd);
                obsYtAll{sInd} = ytDS;
                ytFinInds{sInd} = finInds;
            end
        else
            obsYtAll{sInd} = nan(size(ytOrig{sInd}));
        end

        if ntOut
            obsNtAll{sInd} = ntOrig{sInd};
        else
            obsNtAll{sInd} = nan(size(ntOrig{sInd}));
        end
    end
    if nSegs == 1 && noCell
        obsYtAll = obsYtAll{1};
        obsNtAll = obsNtAll{1};
        ytFinInds = ytFinInds{1};
    end
end


function [ytDS,ytFinInds] = downsampCont(ytOrig,ds,ytStartInd)
    if ds > 1 || ytStartInd ~= 1
        ytDS = nan(size(ytOrig));
        ytFinInds = round(ytStartInd:ds:length(ytOrig));
        ytDS(:,ytFinInds) = ytOrig(:,ytFinInds);
    else
        ytDS = ytOrig;
        ytFinInds = ytStartInd:length(ytOrig);
    end
end