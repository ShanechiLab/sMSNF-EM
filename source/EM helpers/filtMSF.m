function [xUpd,KUpd,xPred,KPred,KPredi] = filtMSF(A,C,Q,CRiC,CRi,alpha,beta,x,K,y,N,sf)
% FILTMSF Decode latent states using multiscale filter (MSF)
%  One step of the MSF [1].
%  INPUTS:
%  x: previous updated latent state
%  K: covariance of previous estimate
%  y: current Gaussian observation
%  N: current Poisson observation
%  A: dynamics matrix
%  Q: latent noise covariance
%  alpha: baseline firing alphas
%  beta: modulation depth betas
%  sf: scaling factor for multiscale
%
%  [1] Hsieh et al, "Multiscale modeling and decoding algorithms for
%      spike-field activity", Journal of Neural Engineering, 2019
%
%  Author: DongKyu Kim, Apr 2026, dongkyuk(at)usc(dot)edu
%
    xPred = A*x;
    Ktmp = A*K*A' + Q;
    KPred = 0.5*(Ktmp + Ktmp');
    KPredi = KPred^-1;
    
    Kesum = zeros(size(Q));
    xsum = zeros(size(x));
    if all(isfinite(N))
        cifdel = exp(alpha' + beta'*xPred);
        xtemp =  beta * (N - cifdel);
        Ktemp = (cifdel' .* beta) * beta';
        if isfinitereal(xtemp) && isfinitereal(Ktemp)
            xsum = xtemp;
            Kesum = Ktemp;
        end
    end
    if all(isfinite(y))
        Kesum = Kesum + sf * CRiC;
        ytUpd = CRi*(y-C*xPred);
    else
        ytUpd = zeros(size(x));
    end
    
    KupdTemp = (KPredi + Kesum)^-1;     
    
    if any(isnan(KupdTemp(:)))
        KUpd = KPred;
    else
        KUpd = 0.5*(KupdTemp + KupdTemp');
    end
    
    xUpd = xPred + KUpd*(sf * ytUpd + xsum);

end

