function normx = myVecnorm(x,dim)
    if ~exist('dim','var')
        dim = 1;
    end
    if dim == 2
        x = x';
    end
    
    x2 = x.^2;
    normx = sqrt(sum(x2,1));
    
    if dim == 2
        normx = normx';
    end
end

