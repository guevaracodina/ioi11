function F = ioi_nlinfit(x,d,color,c1,include_flow)
warning('off')
%nonlinear fit - two steps
%fit the hemodynamic increase (decrease for HbR)
if strcmp(color.eng(c1),color.HbR)
    d = -d;
    p1 = [5.8 0.3];
    p2 = [10 0.5 1/8];
else
    p1 = [4.3 0.5];
    p2 = [9 0.5 1/8];
end
if ~strcmp(color.eng(c1),color.flow) || include_flow
    try
        %options = statset('Robust','on','DerivStep',eps^(1/3),'TolX',1e-8,'Display','Iter','TolFun',1e-8);
        [beta1,R1,J1,COVB1,mse1] = nlinfit(x,d/sum(d),@ioi_hrf1,p1); %,options);
        F.R1 = R1;
        F.J1 = J1;
        F.COVB1 = COVB1;
        F.mse1 = mse1;
        F.FitIncreaseOK = 1;
    catch
        %could not fit
        beta1 = p1;
        F.FitIncreaseOK = 0;
    end
    %fit the undershoot
    %p   = [beta(1) 10 beta(2) 1 6];
    d2 = d/sum(d) - ioi_hrf1(beta1,x);
    
    try
        %options = statset('Robust','on','DerivStep',eps^(1/3),'TolX',1e-8,'Display','Iter','TolFun',1e-8);
        [beta2,R2,J2,COVB2,mse2] = nlinfit(x,d2,@ioi_hrf2,p2); %,options);
        F.R2 = R2;
        F.J2 = J2;
        F.COVB2 = COVB2;
        F.mse2 = mse2;
        F.FitDecreaseOK = 1;
    catch
        beta2 = p2;
        F.FitDecreaseOK = 0;
    end
    % figure; plot(d/sum(d),'k'); hold on; plot( ioi_hrf1(beta1,x),'r'); hold on; plot( ioi_hrf2(beta2,x),'g'); hold on; plot( ioi_hrf1(beta1,x)+ioi_hrf2(beta2,x),'b');
    
    p = [beta1(1) beta2(1) beta1(2) beta2(2) beta2(3)];
    %Then try to fit both
    %together
    try
        %options = statset('Robust','on','DerivStep',eps^(1/3),'TolX',1e-8,'Display','Iter','TolFun',1e-8);
        [beta,R,J,COVB,mse] = nlinfit(x,d/sum(d),@ioi_hrf,p); %,options);
        %confidence intervals
        [ypred dlt] = nlpredci(@ioi_hrf,max(x),beta,R,'Covar',COVB);
        F.R = R;
        F.J = J;
        F.COVB = COVB;
        F.mse = mse;
        F.FitBothOK = 1;
    catch
        beta = p;
        F.FitBothOK = 0;
    end
    %Fill F
    F.beta = beta;
    F.beta1 = beta1;
    F.beta2 = beta2;
    %Calculate the fitting curve:
    if strcmp(color.eng(c1),color.HbR)
        F.yp = -sum(d)*(ioi_hrf(beta,x));
    else
        F.yp = sum(d)*(ioi_hrf(beta,x));
    end
else
    F = [];
end
warning('on')
end