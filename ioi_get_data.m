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
if includeHbR
    cl = cl+1;
    Y(:,cl) = ROI{r1}{s1,cHbR};
end
if includeHbT
    cl = cl+1;
    Y(:,cl) = ROI{r1}{s1,cHbR}+ROI{r1}{s1,cHbO};
end
if includeFlow
    cl = cl+1;
    Y(:,cl) = ROI{r1}{s1,cFlow};
end
HPF = M.HPF;
if HPF.hpf_butter_On
    Y = ButterHPF(1/M.TR,HPF.hpf_butter_freq,HPF.hpf_butter_order,Y);
end
M.Y = Y;
end