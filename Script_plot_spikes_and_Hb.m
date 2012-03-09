figure; r=1; c=5; s=5; sf = 5; ht = 10; d = ROI{r}{s,5}+ROI{r}{s,6}; ROI{r}{s,c}; n= length(d); 
d = ButterLPF(5,0.25,3,d);
d = ButterHPF(5,0.01,3,d);
[pkhMax pkMax] = findpeaks(d,'minpeakdistance',5*sf);
[pkhMin pkMin] = findpeaks(-d,'minpeakdistance',5*sf);
if (r==1 && c == 6) || (r==1 && c==5), pmax = 5; else pmax = 6; end
pkhMax = pkhMax(pmax:end-1)/sf;
pkMax = pkMax(pmax:end-1)/sf;
if r==2 && c == 6, pmin = 5; else pmin = 6; end
pkhMin = pkhMin(pmin:end-1)/sf;
pkMin = pkMin(pmin:end-1)/sf;
plot(linspace(0,n/sf,n),d,'k'); hold on; 
pk = IOI.sess_res{s}.onsets{1}; 
pk = pk([1:16 18:33 35:39 41:51 53]);
stem(pk,ht* ones(1,length(pk)),'r'); 
stem(pkMax,ht*ones(1,length(pkMax)),'g'); 
stem(pkMin,ht*ones(1,length(pkMin)),'b'); 
dp1 = pkMax-pk; std(dp1)
dp2 = pk-pkMin; std(dp2)
dpk = diff(pk); std(dpk)
%
dpk = diff(pk);
dph = pkhMax(1:end-1)-pkhMin(2:end);
figure; scatter(dpk,dph,'k'); hold on; scatter(dpk(1:end-1),dph(2:end),'r'); scatter(dpk(2:end),dph(1:end-1),'b'); 