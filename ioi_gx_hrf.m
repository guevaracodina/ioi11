function y = ioi_gx_hrf(x,u,p,M)
dt = 1; %M.Y.dt;
y = p(6)*(spm_Gpdf(x,p(1)/p(3),dt/p(3)) - spm_Gpdf(x,p(2)/p(4),dt/p(4))/p(5));
