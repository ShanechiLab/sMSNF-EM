function [nrmse,rmse,den] = getNRMSE(x,xhat,mode)
% getNRMSE: get NRMSE values = norm(est-true)/norm(mean-true)
% mode: vector for RMSE across vector, no mode for RMSE for each row
    dimXt = size(x,1);
    error = x - xhat;
    xbar = mean(x,2);
    if ~exist('mode','var')
        mode = '';
    end
    
    if any(strcmpi(mode,{'vector','vec','v'}))
        squaredNormError = sum(abs(error).^2,1);
        rmse = sqrt(mean(squaredNormError));
        den = sqrt(mean(sum(abs(x-xbar).^2,1)));
        nrmse = rmse/den;
    else
        nrmse = zeros(dimXt,1);
        rmse = zeros(dimXt,1);
        den = zeros(dimXt,1);
        for k = 1:dimXt
            rmse(k) = sqrt(mean(error(k,:).^2));
            den(k) = sqrt(mean((x(k,:)-xbar(k,:)).^2));
            nrmse(k) = rmse(k)/den(k);
        end
    end
end

