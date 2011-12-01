function hrf = ioi_hrf(p,u)
% Karl Friston
% $Id: spm_hrf.m 3716 2010-02-08 13:58:09Z karl $
dt = 1;%
hrf = spm_Gpdf(u,p(1)/p(3),dt/p(3)) - spm_Gpdf(u,p(2)/p(4),dt/p(4))/p(5);
hrf = hrf/sum(hrf);
end