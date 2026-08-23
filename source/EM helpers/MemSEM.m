classdef MemSEM < handle
    
    properties
        theta
        Ricell
        logDetRcell
        CRiCcell
        CRicell
        
        xMats
        KMats
        PStDec
        xTmpBch
        KTmpBch
        xTmpBch2
        KTmpBch2
        xPrdMixY
        xPrdY
        PStPrdY
        yTmpPrd
        PNtCondSt
        xSmtBch
        KSmtBch
    end
    
    methods
        function obj = MemSEM(theta,dimXt,dimSt,dimYt,dimNt,...
                              maxLen,paramSet)
            tlen = maxLen;
            obj.xMats = zeros(dimXt,tlen,dimSt,2);
            obj.KMats = zeros(dimXt,dimXt,tlen,dimSt,3);
            obj.PStDec = zeros(dimSt,tlen);
            obj.xTmpBch = zeros(dimXt,dimSt);
            obj.KTmpBch = zeros(dimXt,dimXt,dimSt);
            obj.xTmpBch2 = zeros(dimXt,dimSt^2);
            obj.KTmpBch2 = zeros(dimXt,dimXt,dimSt^2);
            obj.xPrdMixY = zeros(dimXt,dimSt);
            obj.xPrdY = zeros(dimXt,dimSt);
            obj.PStPrdY = zeros(dimSt,dimSt);
            obj.yTmpPrd = zeros(dimYt,dimSt);
            obj.PNtCondSt = zeros(dimNt,dimSt);
            obj.xSmtBch = zeros(dimXt,dimSt);
            obj.KSmtBch = zeros(dimXt,dimXt,dimSt);
            
            if ~exist('paramSet','var') || isempty(paramSet)
                paramSet = struct;
            end
            
            fields = {};
            defaults = {};
            paramSet = insertStructDefaults(paramSet,fields,defaults);
            
            obj.updateTheta(theta);
        end
        
        function obj = updateTheta(obj,theta)
            obj.theta = theta;
            obj.Ricell = cell(size(theta.Rcell));
            for k = 1:length(obj.Ricell)
                R = theta.Rcell{k};
                C = theta.Ccell{k};
                obj.Ricell{k} = theta.Rcell{k}^-1;
                try
                    obj.logDetRcell{k} = logdet(R,'chol');
                catch
                    obj.logDetRcell{k} = logdet(R);
                end
                CRi = C'*obj.Ricell{k};
                obj.CRicell{k} = CRi;
                CRiC = CRi*C;
                obj.CRiCcell{k} = 0.5*(CRiC + CRiC');
            end
        end
    end 
end