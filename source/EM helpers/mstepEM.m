function thetaNew = mstepEM(sumO,paramSet,thetaOld)
% MSTEPSEM M-step of switch multiscale EM, smooth brain and regime states
%  Estimate brain and regime states for a swiching dynamical system for
%  switch EM parameterized by theta from multiscale observations. Requires
%  3 objects initialized in fitSwitchMultiscaleEM.
%
%  INPUTS:
%  sumO: SumSEM object with collected expected value terms from E-step.
%  prmS: struct of estimation settings
%     FIELDS (main):
%     statObs:  bool for fixed observation parameters across regimes
%     statDyn:  bool for fixed dynamic parameters across regimes
%     mWeight:  double for how much to mix in previous theta into new
%               theta. 1(default) for no mixing, 0 for no new theta
%     mDiagQ:   str for whether to contrain Q matrices to diagonals. Choose
%               between 'none'(default) or 'all'.
%  thetaOld: struct containing system parameters from previous iteration
%  OUTPUT:
%  thetaNew: struct containing new system parameters
% 
%
%  Author: DongKyu Kim, Apr 2026, dongkyuk(at)usc(dot)edu
    if ~exist('paramSet','var')
        paramSet = struct;
    end
    fields = {'mDiagQ','mDiagR','statObs',...
              'mWeight','statDyn','rescaleInit'};
    defaults = {'none','none',false,...
                1,false,1};
    paramSet = insertStructDefaults(paramSet,fields,defaults);
    
    dimXt = size(thetaOld.Acell{1},1);
    dimSt = length(thetaOld.sInit);
    dimYt = size(thetaOld.Ccell{1},1);
    dimNt = size(thetaOld.betaCell{1},2);
    
    xDiffSum = sumO.xDiffSum;
    xAuM1Sum = sumO.xAuM1Sum;
    xAutoSum = sumO.xAutoSum;
    PStSmtSum = sumO.PStSmtSum;
    x0SmtSum = sumO.x0SmtSum;
    K0SmtSum = sumO.K0SmtSum;
    sTranNumSum = sumO.sTranNumSum;
    sTranDenSum = sumO.sTranDenSum;
    sInitSum = sumO.sInitSum;
    nSegs = sumO.nSegs;
    
    if ~isempty(sumO.yxSum)
        yxSum = sumO.yxSum;
        xAutoYSum = sumO.xAutoYSum;
        yProdSum = sumO.yProdSum;
        PStSmtSumY = sumO.PStSmtSumY;
    end
    
    if ~isempty(sumO.xTot)
        xTot = sumO.xTot;
        KTot = sumO.KTot;
        PTot = sumO.PTot;
        ntTot = sumO.ntTot;

        alphaCur = thetaOld.alphaCell;
        betaCur = thetaOld.betaCell;

        maximizeFunc = @maximizePoiLikeli;
        optimOpt = optimoptions('fminunc',...
                                'Algorithm','trust-region',... 
                                'SpecifyObjectiveGradient',true,...
                                'HessianFcn','objective',... 
                                'Display','off',...
                                'OptimalityTolerance',1e-5,...
                                'StepTolerance',1e-5,...
                                'MaxIterations',200);
    end

    Acell = cell(dimSt,1); Qcell = cell(dimSt,1); 
    Ccell = cell(dimSt,1); Rcell = cell(dimSt,1);
    alphaCell = cell(dimSt,1); betaCell = cell(dimSt,1);

    if paramSet.statDyn
        A = sum(xDiffSum,3)*(sum(xAuM1Sum,3)^-1);
        Acell(:) = {A};
        
        Qnum = sum(xAutoSum,3) - A*sum(xDiffSum,3)' ...
               - sum(xDiffSum,3)*A' + A*sum(xAuM1Sum,3)*A';
        Qden = sum(PStSmtSum);
        Qtemp = Qnum * Qden^-1;
        switch paramSet.mDiagQ
            case 'none'
                Q = 0.5*(Qtemp + Qtemp');
            otherwise
                Q = diag(diag(Qtemp));
        end
        Qcell(:) = {Q};
    else
        for k = 1:dimSt
            A = xDiffSum(:,:,k) * xAuM1Sum(:,:,k)^-1;
            Acell{k} = A;

            Qnum = xAutoSum(:,:,k) - A*xDiffSum(:,:,k)' ...
                   - xDiffSum(:,:,k)*A' + A*xAuM1Sum(:,:,k)*A';
            Qden = PStSmtSum(k);
            Qtemp = Qnum * Qden^-1;
            switch paramSet.mDiagQ
                case 'all'
                    Qcell{k} = diag(diag(Qtemp));
                otherwise
                    Qcell{k} = 0.5*(Qtemp + Qtemp');
            end
        end
    end
    if paramSet.statObs == true
        if all(xAutoYSum(:)==0)
            Ccell(:) = {sum(yxSum,3)};
        else
            Ccell(:) = {sum(yxSum,3)*sum(xAutoYSum,3)^-1};
        end
        C = Ccell{1};
        Rnum = sum(yProdSum,3) - C*sum(yxSum,3)' ...
               - sum(yxSum,3)*C' + C*sum(xAutoYSum,3)*C';
        Rden = sum(PStSmtSumY);
        if all(Rnum(:)==0)
            Rtemp = {eye(dimYt)};
        else
            Rtemp = Rnum*Rden^-1;
        end
        if ~strcmp(paramSet.mDiagR,'none')
            Rcell(:) = {diag(diag(Rtemp))};
        else
            Rcell(:) = {0.5*(Rtemp+Rtemp')};
        end
    else
        for k = 1:dimSt
            if all(all(xAutoYSum(:,:,k) == 0))
                Ccell{k} = yxSum(:,:,k);
            else
                Ccell{k} = yxSum(:,:,k) * ...
                           xAutoYSum(:,:,k)^-1;
            end
            C = Ccell{k};
            Rnum = yProdSum(:,:,k) - C*yxSum(:,:,k)' ...
                 - yxSum(:,:,k)*C' + C*xAutoYSum(:,:,k)*C';
            Rden = PStSmtSumY(k);
            if all(Rnum(:) == 0)
                Rtemp = eye(dimYt);
            else
                Rtemp = Rnum * Rden^-1;
            end
            switch paramSet.mDiagR
                case 'all'
                    Rcell{k} = diag(diag(Rtemp));
                otherwise
                    Rcell{k} = 0.5*(Rtemp + Rtemp');
            end
        end
    end
    if paramSet.statObs == false
        for k = 1:dimSt
            inds = sumO.incInds & (sumO.regInds==k);
            [aNew,bNew] = maximizeFunc(ntTot(:,inds),...
                                xTot(:,inds),PTot(:,inds),...
                                KTot(:,:,inds),...
                                alphaCur{k},betaCur{k},...
                                optimOpt);
            alphaCell{k} = aNew;
            betaCell{k} = bNew;
        end
    else
        inds = sumO.incInds;
        [aNew,bNew] = maximizeFunc(ntTot(:,inds),...
                                xTot(:,inds),PTot(:,inds),...
                                KTot(:,:,inds),...
                                alphaCur{1},betaCur{1},...
                                optimOpt);
        alphaCell(:) = {aNew};
        betaCell(:) = {bNew};
    end
    
    sTranPre = sTranNumSum .*(1./sTranDenSum');
    sInitPre = (1/nSegs).*sInitSum;
    x0CovPre = (1/nSegs).*K0SmtSum;

    thetaNew = struct;
    thetaNew.Acell = Acell;
    thetaNew.Ccell = Ccell;
    thetaNew.Qcell = Qcell;
    thetaNew.Rcell = Rcell;
    thetaNew.alphaCell = alphaCell;
    thetaNew.betaCell = betaCell;
    thetaNew.sTran = (1./sum(sTranPre,1)).* sTranPre;
    
    if paramSet.rescaleInit ~= 1
        tmp = (1./sum(sInitPre)) .* sInitPre;
        factor = paramSet.rescaleInit;
        thetaNew.sInit = rescaleDist(tmp,factor);
    else
        thetaNew.sInit = (1./sum(sInitPre)) .* sInitPre;
    end
    
    thetaNew.x0Mean = (1/nSegs).*x0SmtSum;
    thetaNew.x0Cov = 0.5*(x0CovPre + x0CovPre');
    

    if isfield(paramSet,'mWeight') && exist('thetaOld','var')
        thetaNew = mixTheta(thetaOld,thetaNew,paramSet.mWeight);
    end
end