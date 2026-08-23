function C = getRandCInit(dimYt,dimXt,paramSet)
    if ~exist('paramSet','var')
        paramSet = struct;
    end
    variableList = {'outMagMin','noNorm'};
    defaultVals = {1,false};
    paramSet = insertStructDefaults(paramSet,variableList,defaultVals);
    
    if ~isfield(paramSet,'outMagMax')
        paramSet.outMagMax = paramSet.outMagMin + 1;
    end
    Cdir = mvnrnd(zeros(dimYt,dimXt),eye(dimXt));
    if paramSet.noNorm == true
        C = Cdir;
    else
        low = paramSet.outMagMin;
        high = paramSet.outMagMax;
        Cmags = getUnifRand(low,high,dimYt,1);
        C = (Cmags./myVecnorm(Cdir,2)).*Cdir;
    end
end