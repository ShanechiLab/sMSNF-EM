function [mdlS,needNewMdl,needToRun,startIter,newIter] = loadModel(pathModel,paramSet,nIter)
    
    if ~exist('paramSet','var')
        paramSet = struct;
    end
    if isfield(paramSet,'verbose')
        verbose = paramSet.verbose;
    else
        verbose = false;
    end

    mdlS = [];
    needNewMdl = true;
    needToRun = true;
    startIter = 1;
    newIter = nIter;
    if isfile(pathModel)
        loaded = false;
        try
            allMdl = load(pathModel);
            loadMdl = allMdl.mdlEM;
            if isstruct(loadMdl) && isfield(loadMdl,'thetaCell') ...
                                 && isfield(loadMdl,'static')
                loaded = true;
            end
        catch
            if verbose
                fprintf(' failed to load ');
            end
        end
        if loaded
            [done,finIt,newIt] = checkLoadProgress(loadMdl,nIter);
            needNewMdl = false;
            if done
                needToRun = false;
                ldIter = length(loadMdl.thetaCell)-1;
                if finIt < ldIter
                    mdlS = trimModel(loadMdl,finIt);
                    save(pathModel,'-struct','mdlS','-v7.3');
                else
                    mdlS = loadMdl;
                end
            else
                startIter = finIt + 1;
                if ~isempty(newIt)
                    nIter = newIt;
                    mdlS = expandModel(loadMdl,nIter);
                else
                    mdlS = loadMdl;
                end
            end
        end
        if needNewMdl
            [dirPath,filename,ext] = fileparts(pathModel);
            newName = [filename, '_old'];
            backupPath = fullfile(dirPath,[newName, ext]);
            try
                movefile(pathModel,backupPath);
            catch
                if verbose
                    fprintf(' failed to backup ');
                end
            end
        end
    else
        warning('model file doesnt exist');
    end
end

