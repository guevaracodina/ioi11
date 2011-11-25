function M=ioi_set_priors(M)
%load priors
switch M.PS.PhysioModel_Choice
    case 0 %Buxton-Friston
        [pE,pC] = ioi_SCKS_priors_base_case(M.m,5);
        pE(end)=0.025; % put a better approximation around effFlow to extract good estimate of covariance
        M.name={'Signal_decay', 'Feedback', 'Transit_time',...
            'Exponent', 'Extraction',  'log_signal_ratio' 'eff'};
        %M = ioi_SCKS_priors_Buxton(M,5);
    case 1 %Zheng-Mayhew
        [pE,pC] = ioi_priors_ZM(M.m);
        M.name={'Signal_decay', 'Feedback', 'Transit_time',...
            'Exponent', 'Extraction','Vascular_tone' 'Gain_parameter', 'log_signal_ratio' 'eff'};
        pE(pE==0)=min(pE(pE~=0));%obligatoire que cela vaille pas z�ro
    case 2 %Huppert1
        M.name={'ks'    'kr'    'kx'    'km'    'Vw0'    'beta'    'Ra0'    'effCMRO'    'effFlow'};
        %%% fixed physical constants
        BH.a=23400;%fixed - from dissociation curve of hemoglobin
        BH.b=150;%fixed - from dissociation curve
        BH.Hn=1.39;%fixed
        BH.rho=1.04; %brain density in g/mL
        BH.omega=44.6429;
        
        %ajustable variables
        BH.alphap=3.9e-5;
        BH.alphat=1.18e-4;
        BH.beta=[2.33  ];
        BH.delta=0.5;
        BH.Fin0=.0148; %Initial flow, normally derived from model
        BH.HGB=.16;%fixed, but depends on subject
        
        BH.ks=0.69; % 0.69;
        BH.kr=0.27;% .17;
        BH.kx=0.96;%; %gros retour lent    kx et km gros(mont�e rapide)
        BH.km=0.27;%gros retour rapide
        % a.kr1=0.9;
        % a.km1=0.5;
        
        BH.M0=0.0387;
        BH.PaO2=75; % 75 ou 100 � voir
        BH.P0=47;
        BH.PtO2=16;
        BH.Ra0=4.23e3;
        BH.Rw0=2.53e3;
        BH.Vt=0.9483;
        BH.Vw0=0.0388 ;%Vvasc*75%
        BH.Va0=0.0388 /3; %Vvasc*25%
        BH.zeta1=1;
        BH.d0=300;
        BH.d1=.5;
        
        BH.effFlow= 0.4; %[.4 0.2];
        BH.effCMRO=0.051; %[.051 .006];
        
        BH.HbT0=0.107;
        BH.phiw0=0.75;
        BH.scFlow=1 ;
        BH.scHbR=1;
        
        BH.noCompliance=0;
        if BH.noCompliance==1
            BH.effCMRO=[0.00000 0];
            BH.ks=069;
            BH.kr=170;
            BH.Vw0=BH.Fin0;
        end
        for i1=1:length(M.name)
            var=BH.(M.name{i1});
            pE(i1,1)=var(1);
            if length(var)>=2
                pC(i1,i1)=var(2);
            else
                pC(i1,i1)=.2^2*var(1)^2;
            end
        end
        
       
        M.PS.BH = BH;
    otherwise
end
M.PS.pE=pE;
M.PS.pC=pC;
M.pE = pE;
M.pC = pC;