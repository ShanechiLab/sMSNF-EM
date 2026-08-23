function bKb = elemVecMult(b,Kt)

    [dimXt,nTimes] = size(b);
    [ord1,ord2,tlen] = size(Kt);
    
    if ord1 ~= dimXt || ord2 ~= dimXt
        error('b and Kt must have correct dim');
    end

    if nTimes > 1
        bKb = zeros(nTimes,tlen);
        for k = 1:nTimes
            bk = b(:,k);
            bKb(k,:) = elemVecMult(bk,Kt);
        end
    else
        Kb = reshape(b' * reshape(Kt,dimXt,[]),dimXt,[]);
        bKb = b'*Kb;
    end
    
end