function M = ioi_HDM_priors_Hu(M,h)
% returns priors for a hemodynamic dynamic causal model
% FORMAT [pE,pC] = spm_hdm_priors(m,[h])
% m   - number of inputs
% h   - number of hemodynamic modes (default = 3)
%
% pE  - prior expectations
% pC  - prior covariances
%

%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_hdm_priors.m 2050 2008-09-05 19:15:50Z klaas $

% biophysical parameters with prior expectation
%---------------------------------------------------------------------------
if ~isfield(M,'name')
    M.name={'ks'    'kr'    'kx'    'km'    'Vw0'    'beta'    'Ra0'    'effCMRO'    'effFlow'};
end

M.nameIni={'ks'    'kr'    'kx'    'km'    'Vw0'    'beta'    'Ra0'    'effCMRO'    'effFlow'};

%%% fixed constant
a.a=23400;%fixe
a.b=150;%fixe
a.Hn=1.39;%fixe
a.rho=1.04;
a.omega=44.6429;

%ajustable variables
a.alphap=3.9e-5;
a.alphat=1.18e-4;
a.beta=[2.33  ];
a.delta=0.5;
a.Fin0=.0148; % normally dérived from modelderived normaly 
a.HGB=.16;%fixe peut changer entre les sujets

a.ks=0.69; % 0.69;
a.kr=0.27;% .17;
a.kx=0.96;%; %gros retour lent    kx et km gros(montée rapide)
a.km=0.27;%gros retour rapide
% a.kr1=0.9;
% a.km1=0.5;

a.M0=0.0387;
a.PaO2=75; % 75 ou 100 à voir
a.P0=47;
a.PtO2=16;
a.Ra0=4.23e3;
a.Rw0=2.53e3;
a.Vt=0.9483;
a.Vw0=0.0388 ;%Vvasc*75%
a.Va0=0.0388 /3; %Vvasc*25% 
a.zeta1=1;
a.d0=300;
a.d1=.5;

a.effFlow=[.4 0.2];
a.effCMRO=[.051 .006];

a.HbT0=0.107;
a.phiw0=0.75;
 a.scFlow=1 ; 
 a.scHbR=1; 
 
 if isfield(M,'scaleGx') && ~isnumeric(M.scaleGx)
    a.scFlow=M.scaleGx.Flow;
     a.scHbR=M.scaleGx.HbR;
 end
 
a.doKm1=0;
a.normalise=1;
a.noCompliance=0;
if a.noCompliance==1
    a.effCMRO=[0.00000 0];
    a.ks=069;
    a.kr=170;
    a.Vw0=a.Fin0;
end

% put value to 0 (see IOI_fx_hdm_Hu2.m)
if isfield(M,'kxkm') && M.kxkm==0
    a.kx=0; a.km=0;
end
if isfield(M,'krkm') && M.krkm==0
    a.kr=0; a.km=0;
end


%tests
a.alphat=a.alphap;%%%%555 pas physique  disp('remettre alpaht')
a.effCMRO=.10 ;

%form Variable in M
b=fieldnames(a);
for i1=1:length(b)
    M.var.(b{i1}) =a.(b{i1})(1);
end


%form covariance
if a.normalise==1  %normalisé
    for i1=1:length(M.name)
        var=a.( M.name{i1});
        pE(i1,1)=(1);
        if length(var)>=2
            pC(i1,i1)=(var(2)/var(1))^2;
        elseif isfield(M,'covarianceIni') && M.covarianceIni>0
               pC(i1,i1)= M.covarianceIni;
        else
            pC(i1,i1)=.2^2*1^2;
        end
    end
elseif  1 %variance donnée par lles valeurs par défault
    for i1=1:length(M.name)
        var=a.(M.name{i1});
        pE(i1,1)=var(1);
        if length(var)>=2
            pC(i1,i1)=var(2);
        else
            pC(i1,i1)=.2^2*var(1)^2;
        end
    end
elseif 0  %variance sommaire    
    pC   = .4^2*pE'*pE;
    pC=diag(pE.*pE);    
end


M.pE=pE;
M.pC=pC;



