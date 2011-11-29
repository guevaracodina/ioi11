function hrf = ioi_hrf2(p,u)
% Karl Friston
% $Id: spm_hrf.m 3716 2010-02-08 13:58:09Z karl $
% modelled hemodynamic response function - only the undershoot
%--------------------------------------------------------------------------
%dt  = u(2)-u(1); % RT/fMRI_T;
%u   = [0:(p(7)/dt)] - p(6)/dt;
dt = 1;%
hrf = -spm_Gpdf(u,p(1)/p(2),dt/p(2))/p(3);
%hrf = hrf/sum(hrf);
end