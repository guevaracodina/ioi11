function M = ioi_nlsi(M)
% nonlinear system identification of a MIMO system
% FORMAT [Ep,Cp,Eh,K0,K1,K2,M0,M1,L1,L2] = spm_nlsi(M,U,Y)
% FORMAT [K0,K1,K2,M0,M1,L1,L2]          = spm_nlsi(M)
%
% Model specification
%--------------------------------------------------------------------------
% M.f     - dx/dt = f(x,u,P,M)  {function string or m-file}
% M.g     - y     = g(x,u,P,M)  {function string or m-file}
%
% M.pE    - (p x 1)   Prior expectation of p model parameters
% M.pC    - (p x p)   Prior covariance for p model parameters
%
% M.x     - (n x 1)   intial state x(0)
% M.m     - m         number of inputs
% M.n     - n         number of states
% M.l     - l         number of outputs
% M.N     -           kernel depth
% M.dt    -           kernel resolution {secs}
%
% System inputs
%--------------------------------------------------------------------------
% U.u   - (v x m)   m inputs
% U.dt  -           sampling interval for inputs
%
% System outputs
%--------------------------------------------------------------------------
% Y.y   - (v x l)   l outputs
% Y.X0  - (v x c)   Confounds or null space
% Y.dt  -           sampling interval for outputs
% Y.Q   -           observation error precision components
%
% Model Parameter estimates - conditional moments
%--------------------------------------------------------------------------
% Ep    - (p x 1)           conditional expectation  E{P|y}
% Cp    - (p x p)           conditional covariance   Cov{P|y}
% Eh    - (v x v)           conditional log-precision
%
% System identification     - Volterra kernels
%--------------------------------------------------------------------------
% K0    - (l x 1)             = k0(t)         = y(t)
% K1    - (N x l x m)         = k1i(t,s1)     = dy(t)/dui(t - s1)
% K2    - (N x N x l x m x m) = k2ij(t,s1,s2) = d2y(t)/dui(t - s1)duj(t - s2)
%
% System identification     - Bilinear approximation
%--------------------------------------------------------------------------
% M0    - (n x n)           df/dq
% M1    - {m}(n x n)        d2f/dqdu
% L1    - (l x n)           dg/dq
% L2    - {l}(n x n)        d2g/dqdq
%
%__________________________________________________________________________
%
% Returns the moments of the posterior p.d.f. of the parameters of a 
% nonlinear MIMO model under Gaussian assumptions
%
%              dx/dt  = f(x,u,P)
%                y    = g(x,u,P) + e                               (1)
%
% evaluated at x(0) = x0, using a Bayesian estimation scheme with priors
% on the model parameters P, specified in terms of expectations and 
% covariance. The estimation uses a Gauss-Newton method with MAP point 
% estimators at each iteration.  Both Volterra kernel and state-space 
% representations of the Bilinear approximation are provided.
% The Bilinear approximation to (1), evaluated at x(0) = x and u = 0 is:
%
%       dq/dt = M0*q + u(1)*M1{1}*q + u(2)*M1{2}*q + ....
%        y(i) = L1(i,:)*q + q'*L2{i}*q;
%
% where the states are augmented with a constant
%
%        q(t) = [1; x(t) - x(0)]
%
% The associated kernels are derived using closed form expressions based
% on the bilinear approximation.
%
%--------------------------------------------------------------------------
% If the inputs U and outputs Y are not specified the model is simply
% characterised in terms of its Volterra kernels and Bilinear
% approximation expanding around M.pE
%
% see also
% spm_nlsi_GN:   Bayesian parameter estimation using an EM/Gauss-Newton method
% spm_bireduce:  Reduction of a fully nonlinear MIMO system to Bilinear form
% spm_kernels:   Returns global Volterra kernels for a MIMO Bilinear system
%
% SEE NOTES AT THE END OF THIS SCRIPT FOR EXAMPLES
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_nlsi.m 3764 2010-03-08 20:18:10Z guillaume $

% Gauss-Newton/Bayesian/EM estimation
%======================================================================
[Ep,Cp,Eh,F,k,MSE,MSE_HbR] = ioi_nlsi_GN(M,M.U,M.Y);

% Bilinear representation
%==========================================================================
[M0,M1,L1,L2] = spm_bireduce(M,Ep);

% Volterra kernels
%==========================================================================
% time bins (if not specified)
%--------------------------------------------------------------------------

dt   = M.dt;
N    = M.N;

% get kernels
%--------------------------------------------------------------------------
[K0,K1,K2]   = spm_kernels(M0,M1,L1,L2,N,dt);
[dummy,H1] = spm_kernels(M0,M1,M.N,M.dt);

M.Ep = Ep;
M.Cp = Cp;
%M.Ce = Ce;
M.Eh = Eh;
M.K0 = K0;
M.K1 = K1;
M.K2 = K2;
M.M0 = M0;
M.M1 = M1;
M.L1 = L1;
M.L2 = L2;
M.F = F;
M.H1 = H1;
M.k = k;
M.MSE = MSE;
M.MSE_HbR = MSE_HbR;
ioi_HDM_display(M);
