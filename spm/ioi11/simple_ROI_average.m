yO = ROI{1}{5,5};
yR = ROI{1}{5,6};
U = IOI.Sess(5).U;
A = IOI.sess_res{5}.parameters{1}/10000;
wA = 0; %Boolean to weigh by LFP amplitudes
ons = U.ons(diff(U.ons)>10);
wb = 2;
wa = 8;
off = 0;
sf = 5;
aO = zeros(1,sf*wa);
bO = zeros(1,sf*wb);
ons = ons-off;
for i0=1:length(ons)
     
end

