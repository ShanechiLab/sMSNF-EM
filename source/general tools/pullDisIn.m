function stNew = pullDisIn(st,pullFactor)
    dimSt = size(st,1);
    stNew = (1-pullFactor).*st + pullFactor.*((1/dimSt).*ones(size(st)));
end

