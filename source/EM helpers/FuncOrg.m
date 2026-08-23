classdef FuncOrg < handle
    properties
        filtSelect
        theta
        Ricell
        logDetRcell
        CRicell
        CRiCcell
        dimSt
        wts = []
        cubeDsig = []
        wtsPP = []
        ptsPP = []
        sf
    end
    
    methods
        function obj = FuncOrg(theta,paramSet)
            if ~exist('paramSet','var') || isempty(paramSet)
                paramSet = struct;
            end
            fields = {'filtSelect', 'nPts', 'sf'};
            defaults = {'MSNF',5, 1};
            paramSet = insertStructDefaults(paramSet,fields,defaults);
            dimSt = length(theta.sInit);
            obj.dimSt = dimSt;
            obj.resetTheta(theta);
            dimXt = length(theta.Acell{1});
            obj.sf = paramSet.sf;
            if contains(paramSet.filtSelect,'MSNF')
                [wts,dsig] = getCubePts(dimXt,paramSet.nPts);
                obj.wts = wts;
                obj.cubeDsig = dsig;
            end
            obj.filtSelect = paramSet.filtSelect;
            [obj.wtsPP,obj.ptsPP] = getCubePts(dimXt,5);
        end
                
        function obj = resetTheta(obj,theta)
            obj.theta = theta;
            obj.Ricell = cell(size(theta.Rcell));
            for k = 1:obj.dimSt
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
        
        function cif = getCif(obj,s,x)
            cif = exp(obj.theta.alphaCell{s}' ...
                      + obj.theta.betaCell{s}'*x);
        end

        function over = checkOver(obj,s,xSmt,KSmt)
            over = false;
            axb = obj.theta.alphaCell{s}' + ...
                  obj.theta.betaCell{s}'*xSmt;
              
            bxb = 0.5* sum(obj.theta.betaCell{s}'.* ...
                           obj.theta.betaCell{s}'*KSmt,2);
            if any(axb + bxb > 100)
                over = true;
            end
        end

        function [xUpd,KUpd,xPrd,KPrd,KPrdi] = applyFilt(obj,s,x,K,y,n)
            A = obj.theta.Acell{s}; C = obj.theta.Ccell{s};
            Q = obj.theta.Qcell{s}; CRiC = obj.CRiCcell{s};
            CRi = obj.CRicell{s}; a = obj.theta.alphaCell{s};
            b = obj.theta.betaCell{s};
            switch obj.filtSelect
                case 'MSF'
                    [xUpd,KUpd,xPrd,KPrd,KPrdi] = filtMSF(...
                        A,C,Q,CRiC,CRi,a,b,x,K,y,n,obj.sf);
                case 'MSNF'
                    [xUpd,KUpd,xPrd,KPrd,KPrdi] = filtMSNF(...
                        x,K,y,n,obj.wts,obj.cubeDsig,...
                        A,Q,a,b,C,CRiC,CRi,obj.sf);
            end
        end

        function logPStUpd = calcLogPStUpd(obj,PStPrd,s,xP,xU,KPi,KU,y,N)
            if PStPrd ~= 0
                logDetPart = 0.5*(logdet(KU) + logdet(KPi));
                logYtPart = 0;
                logNtPart = 0;
                if all(isfinite(N))
                    cifDel = obj.getCif(s,xU); 
                    logTemp = N'*log(cifDel) - sum(cifDel);
                    if isfinitereal(logTemp)
                        logNtPart = logTemp;
                    end
                end
                if all(isfinite(y))
                    yxD = y - obj.theta.Ccell{s}*xU;
                    logDetPart = logDetPart ...
                                 - 0.5*obj.logDetRcell{s} * obj.sf;
                    logYtPart = -0.5 * yxD'*obj.Ricell{s}*yxD * obj.sf;
                end
                logPrePart = -0.5*(xU-xP)'*KPi*(xU-xP);
                logPStUpd = logDetPart + logYtPart ...
                            + logNtPart + logPrePart ...
                            + log(PStPrd);
            else
                logPStUpd = log(PStPrd);
            end
        end
        function PnoSpk = noSpikeProb(obj,xPrd,KPrd,j)
            cSpk = length(obj.theta.alphaCell{1});
            PnoSpk = zeros(cSpk,1);
            as = obj.theta.alphaCell{j};
            bs = obj.theta.betaCell{j};
            switch obj.filtSelect
	            case 'MSF'
	            otherwise
		            L = chol(KPrd,'lower');
                sigPts = L*obj.ptsPP + xPrd;
                expabx = exp(as' + bs'*sigPts); %cSpk x nPts
                expexpabx = exp(-expabx); %cSpk x nPts
                PnoSpk(:) = nansum(expexpabx.*obj.wtsPP,2);
            end
        end
    end
    
end

