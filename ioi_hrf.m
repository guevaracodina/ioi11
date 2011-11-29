function hrf = ioi_hrf(p,u)
% Karl Friston
% $Id: spm_hrf.m 3716 2010-02-08 13:58:09Z karl $
%persistent p0

%dp = p- [6 16 1 1 6]; %[4 9 3 3 6];
% if dp * dp' > 1
%     a = 1;
% end
% %Set hard boundaries
% if p(1) <= 0
%     p(1) = 1;
% end
% if p(2) <= 0
%     p(2) = 1;
% end
% if p(3) <= 0
%     p(3) = 1;
% end
% if p(4) < 0
%     p(4) = -p(4);
%     p(5) = p(5)/2;
% end
% if p(5) <= 0
%     p(5) = 1;
% end
% if p(1) >= 15
%     p(1) = 5;
% end
% if p(2) >= 30
%     p(2) = 10;
% end
% if p(3) >= 10
%     p(3) = 3;
% end
% if p(4) >= 10
%     p(4) = 3;
% end
% if p(5) >= 20
%     p(5) = 10;
% end
% modelled hemodynamic response function - {mixture of Gammas}
%--------------------------------------------------------------------------
%dt  = u(2)-u(1); % RT/fMRI_T;
%u   = [0:(p(7)/dt)] - p(6)/dt;
dt = 1;%
hrf = spm_Gpdf(u,p(1)/p(3),dt/p(3)) - spm_Gpdf(u,p(2)/p(4),dt/p(4))/p(5);
%hrf = spm_Gpdf(u,p(1)/p(2),dt/p(2));
%hrf = hrf([0:(p(7)/RT)]*fMRI_T + 1);
hrf = hrf/sum(hrf);
end