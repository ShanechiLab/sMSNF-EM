function [xUpd,KUpd,xPrd,KPrd,KPrdi] = filtMSNF(x,K,y,n,wts,dsig,...
                                                   A,Q,as,bs,C,CRiC,CRi,sf)
% FILTMSNF Decode latent states using multiscale numerical integration
% filter (MSNF)
%  One step of the MSNF [1].
%  INPUTS:
%  x: previous updated latent state
%  K: covariance of previous estimate
%  y: current Gaussian observation
%  n: current Poisson observation
%  wts: cubature weights (generated in getCubePts.m)
%  dsig: cubatutre points (generated in getCubePts.m)
%  A: dynamics matrix
%  Q: latent noise covariance
%  as: baseline firing alphas
%  bs: modulation depth betas
%  sf: scaling factor for multiscale
%
%  [1] Kim et al, "Unsupervised learning of multiscale switching dynamical
%      system models from multimodal neural data",  Journal of Neural
%      Engineering, Apr. 2026
%
%  Author: DongKyu Kim, Apr 2026, dongkyuk(at)usc(dot)edu
%
    xPrd = A*x;
    Ktmp = A*K*A' + Q;
    KPrd = 0.5*(Ktmp + Ktmp');
    KPrdi = KPrd^-1;
    
    if ~any(isfinite([y;n]))
        xUpd = xPrd;
        KUpd = KPrd;
    else       
        Kesum = zeros(size(Q));
        xsum = zeros(size(x));
        
        if all(isfinite(n))
            L = chol(KPrd,'lower');
            sigPts = L*dsig + xPrd;
            cifdels = exp(as' + bs'*sigPts);  
            muNcX = cifdels;
            covNcX = cifdels; 
            muNcN = muNcX*wts';
            KXN = (wts.*(sigPts-xPrd)) * (muNcX - muNcN)';
            Ct = KPrdi*KXN;
            nRi = diag(1./(covNcX*wts'));
            G = Ct*nRi;
            Ktmp = G*Ct';
            Kesum = Kesum + 0.5*(Ktmp+Ktmp');
            xsum = G*(n-muNcN);
        end
        if all(isfinite(y))
            Kesum = Kesum + sf * CRiC;
            ytUpd = CRi*(y-C*xPrd);
        else
            ytUpd = zeros(size(x));
        end
        
        KUpdi = KPrdi + Kesum;
        Ktmp = KUpdi^-1;
        KUpd = 0.5*(Ktmp+Ktmp'); 
        xUpd = xPrd + KUpd*(sf* ytUpd + xsum);
        
        try
            chol(KUpd,'lower');
        catch
            [xUpd,KUpd,xPrd,KPrd,KPrdi] = filtLocalLap_mix(...
                        A,C,Q,CRiC,CRi,as,bs,x,K,y,n,sf);
        end
    end
end