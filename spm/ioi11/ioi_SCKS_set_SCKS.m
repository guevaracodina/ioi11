function SCKS = ioi_SCKS_set_SCKS(SCKS,U)
%data and onsets
if ~isempty(U)
    SCKS.pU.v{1} = U.ons; %?
else
    SCKS.pU.v{1} = [];
end
SCKS.pU.r = SCKS.Y; %r;
SCKS.N=64; %default 64
SCKS.IS='spm_int';    
%Create structure M, required for call to spm_DEM_set, and for SCKS
M(1).PS = SCKS.PS;
M(1).n = SCKS.n;
M(1).l = size(SCKS.Y,2);
M(1).E.linear = 0;                          % linear model
M(1).E.s      = 1;                          % smoothness
M(1).E.dt     = SCKS.TR;
M(1).A = SCKSparams.State_annealing; %0.9995;
M(1).Ap = SCKSparams.Parameter_annealing;
% level 1
%------------------------------------------------------------------
% prior expectation
M(1).W  = exp(blkdiag(5,6,9,9)); %exp(12);        % error precision on states?
M(1).W  = exp(blkdiag(5,5,5,5)-4); %exp(12);        % error precision on states?

% level 2
%------------------------------------------------------------------
M(2).l  = exp(0);                                % inputs
M(2).V  = exp(0);                                % with shrinkage priors (on inputs (sV))
M(2).PS.pC = 1;
M(2).PS.pE = 0;
% free parameters
%--------------------------------------------------------------------------
P       = SCKS.PS.pE;                                % true parameters
ip      = [1:length(P)];                          % free parameters
ip      = [];                          % free parameters
if length(ip)==7
    ip(end-1)=[];% remove logsignal;
end
np      = length(P);

%max_vol = size(r,2); % 150;
%r = r(:,1:max_vol);
%v = min(0.3,v); %truncate for display
nD = 1;

M(1).E.dt = M(1).E.dt; %temporarily
M(1).E.dt     = M(1).E.dt/nD; %just for data generation

%SCKS.temps = temps';
M(1,1).E.nN = 30;
M(1,1).E.nD = nD;

cb(1:6,1)= .6; %low bound
cb(1:6,2)= 1.5; %high bound

M(1).ip = ip;  % indices of model parameters to be estimated
M(1).cb = cb;  % option to specify constrain on parameters values [min max]
M(2).v  = 0;   % input initial condition
M(2).V  = 50;   % input noise precison (fixed) %ini 20% plus gros plus petites barres d'erreur sur U
% log-évidence plus élevée. Fréquence de U plus lent. Moins de risque de
% singularité
M(1).V  = exp(3); %observation noise?
M(1).xP = blkdiag(1e-3^2,1e-2^2,1e-1^2,1e-2^2); %eye(4)*1e-3^2;   % state error covariance matrix
M(1).uP = eye(1)*1e-1^2;   % input error covariance matrix
M(1).wP = eye(np)*1e-4^2;  % parameter error covariance matrix % not used if Q=[];
% SCKS.M(1).pC = diag([1e-5; 1e-5; 1e-6; 1e-8; 1e-8; 1e-8; 1e-8]);  % covarinace matrix of paramters noise
M(1).f  = 'ioi_fx_SCKS';  % state equations rewriten for matrix operations
M(1).g  = 'ioi_gx_SCKS';  % observation equations rewriten for matrix operations
M(1).Q  = {speye(M(1).l,M(1).l)}; % if Q is specified then algorithm performs
% estimation of measurement noise covariance
%  SCKS.M(1).Q  = [];     % if presion on measurement noise is known then Q = [];
% disp('MOdif sdgfegwergwergwergrwergwerg qerge qergwergweg')
M(1).Qf      = 'all';  % form of estimation of measurement noise covariance
% (after online VB estimatin); options: [auto,all,min,mean]
% SCKS.M(1).E.nN    = 15;    % max number of iteration of SCKF-SCKS algorithm
M(1).E.nN    = 8;    % max number of iteration of SCKF-SCKS algorithm
% % if isfield(SCKS.M(1),'SCKSScript') && strcmpi(SCKS.M(1).SCKSScript,'SPM_SCKS_phil')
% %     
% %     SCKS.M(1).A=0.9995;
% %     SCKS.M(1).Ap=0.95;
% %     SCKS.M(1).E.nN    = 20;    % max number of iteration of SCKF-SCKS algorithm
% %     % %Add parameters
% %     ip      = [1:length(P)];
% %     SCKS.M(1).ip = ip;
% %     SCKS.M(1).pE = SCKS.M(1).pE(ip);
% %     SCKS.M(1).pC = SCKS.M(1).pC(ip);
% %     SCKS.M(1).wP = SCKS.M(1).wP(ip,ip);
% %     SCKS.M(1).cb = SCKS.M(1).cb(ip,:);
% % end

M(1).E.Itol  = 1e-8;  % convergence tolerance value for SCKF_SCKS algorithm
M(1).E.RM    = [1e2 1e6];  % scaling parmater for Robbins-Monro approximation of
% parameter noise covariance [scaling parameter, max-limit]
M(1).VB.N    = 10;      % max number of VB iteration during one SCKF-SCKS run
M(1).VB.Itol = 1e-6;    % convergence tolerance value for VB algorithm
M(1).VB.l    = .95;  
SCKS.M = M;