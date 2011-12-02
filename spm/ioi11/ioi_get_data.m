function M = ioi_get_data(ROI,M,r1,s1)
%extract data for up to 3 modalities for session s1 at region-of-interest r1
xY = M.PS.xY;
cHbO = xY.cHbO;
cHbR = xY.cHbR;
cFlow = xY.cFlow;
includeHbR = xY.includeHbR;
includeHbT = xY.includeHbT;
includeFlow = xY.includeFlow;
%number of modalities
l = includeHbR + includeHbT + includeFlow;
%data length
ns = length(ROI{r1}{s1,cHbO});
Y = zeros(ns,l);
cl = 0;
HPF = M.HPF;
LPF = M.LPF;
if includeHbR
    cl = cl+1;
    ty = ROI{r1}{s1,cHbR};
    if HPF.hpf_butter_On
        ty = ButterHPF(1/M.TR,HPF.hpf_butter_freq,HPF.hpf_butter_order,ty);
    end
    
    %try, ty = ty/abs(median(ty))-1; end
    Y(:,cl) = ty;
    %Y(:,cl) = ROI{r1}{s1,cHbR}/median(ROI{r1}{s1,cHbR});
end
if includeHbT
    cl = cl+1;
    ty1 = ROI{r1}{s1,cHbR};
    ty2 = ROI{r1}{s1,cHbO};
    if HPF.hpf_butter_On
        ty1 = ButterHPF(1/M.TR,HPF.hpf_butter_freq,HPF.hpf_butter_order,ty1);
        ty2 = ButterHPF(1/M.TR,HPF.hpf_butter_freq,HPF.hpf_butter_order,ty2);
    end
    ty = ty1+ty2;
    %ty = ty/mean(ty)-1;
    Y(:,cl) =ty;
    %Y(:,cl) = ROI{r1}{s1,cHbR}/median(ROI{r1}{s1,cHbR})+ROI{r1}{s1,cHbO}/median(ROI{r1}{s1,cHbO});
end
if includeFlow
    cl = cl+1;
    ty = ROI{r1}{s1,cFlow};
    if HPF.hpf_butter_On
        ty = ButterHPF(1/M.TR,HPF.hpf_butter_freq,HPF.hpf_butter_order,ty);
    end
    %ty = ty/mean(ty)-1;
    Y(:,cl) = ty;
end
% HPF = M.HPF;
% if HPF.hpf_butter_On
%     Y = ButterHPF(1/M.TR,HPF.hpf_butter_freq,HPF.hpf_butter_order,Y);
% end
if LPF.lpf_gauss_On
    K = get_K(1:size(Y,1),LPF.fwhm1,M.TR);
    for i=1:size(Y,2)
        y = Y(:,i)';
        %forward
        y = spm_filter_HPF_LPF_WMDL(K,y')';
        %backward
        y = y(end:-1:1);
        y = spm_filter_HPF_LPF_WMDL(K,y')';
        y = y(end:-1:1);
        Y(:,i) = y;
    end
end
M.Y = Y;
end