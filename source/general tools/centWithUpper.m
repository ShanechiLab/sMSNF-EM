function out = centWithUpper(in,upper)
    maxIn = max(in);
    meanIn = mean(in);
    
    if (0 - meanIn) > (upper - maxIn)
        out = (in - maxIn) + upper;
    else
        out = in - meanIn;
    end
end