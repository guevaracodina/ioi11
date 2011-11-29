function hrf = ioi_hrf1(p,u)
% Karl Friston
% $Id: spm_hrf.m 3716 2010-02-08 13:58:09Z karl $

% modelled hemodynamic response function - only the increase
%--------------------------------------------------------------------------
dt = 1;%
hrf = spm_Gpdf(u,p(1)/p(2),dt/p(2));
hrf = hrf/sum(hrf);
end