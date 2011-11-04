function [y Y] = ioi_HDM_gx(x,u,P,M)
% Simulated BOLD response to input.
% FORMAT [y] = spm_gx_hdm(x,u,P,M)
% y    - IOI response response (%)
% x    - state vector     (see spm_fx_dcm)
% P    - Parameter vector (see spm_fx_dcm)
%__________________________________________________________________________
%
% This function implements the BOLD signal model described in:
%
% Stephan KE, Weiskopf N, Drysdale PM, Robinson PA, Friston KJ (2007)
% Comparing hemodynamic models with DCM. NeuroImage 38: 387-401.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston & Klaas Enno Stephan
% $Id: spm_gx_hdm.m 3812 2010-04-07 16:52:05Z karl $


% biophysical constants for 1.5 T:
%==========================================================================


% exponentiation of hemodynamic state variables
%--------------------------------------------------------------------------

x     = exp(x);

% BOLD signal
%--------------------------------------------------------------------------
if ~exist('M','var')
    M=loadM;
    M=M(1);
end


[a void2 P]=IOI_extraitVar(M,P,x,u);

switch M.Model_Choice
    case {0,4, 5 6 8} %Buxton, 
%         x(3)=(x(3)-1)*1.1+1; testt
%         x(4)=(x(4)-1)*1.1+1;test
        v     = x(3);        q     = x(4);
        fv       = (x(3).^(1/P(4)));%Flow Out
      
        out.HbT = (v-1).*1; %Assume v = HbT = HbO+ HbR
        out.HbR = (q-1).*1; %Assume q = HbR 
%           out.HbT = (v-1)*.9; %Assume v = HbT = HbO+ HbR %test
%          out.HbR = (q-1)*.9; %Assume q = HbR %test
%         Flow = cf*(x(2)-1); %*(1+ y(1)); %speckle=flow in
        out.Flow= ( (x(2)+fv)./2-1); % speckle=mean (flow in , flow out)
        %          Flow     = V0*(k1.*(1 - q) + k2.*(1 - (q./v)) + k3.*(1 - v));  %Bold
        out.HbO=(out.HbT*100e-6-out.HbR*40e-6)./60e-6;
        %         FMRI
        epsi=1; E0= P(5); TE = 0.04;  V0    = 100*0.08;  r0    = 100;  nu0   = 80.6;
        k1    = 4.3.*nu0.*E0.*TE;
        k2    = epsi.*r0.*E0.*TE;
        k3    = 1 - epsi;
        out.bold     = V0*(k1.*(1 - q) + k2.*(1 - (q./v)) + k3.*(1 - v));  
         if  any(M.Model_Choice==[4])
            out.decay=x(6);
         elseif  any(M.Model_Choice==[5])
                    out.decay=x(5);
         end
        
        
    case {1,3} %, Zheng Mayhew and Zheng Decay
        v     = x(3);        q     = x(4);
        fv       = (x(3)^(1/P(4)))/x(5);%Flow Out       
        out.HbT = (v-1); %Assume v = HbT = HbO+ HbR
        out.HbR = (q-1); %Assume q = HbR
%         Flow = (x(2)-1); %*(1+ y(1)); %flow in
        out.Flow= ( (x(2)+fv)/2-1); % speckle=mean (flow in , flow out)
        %             Bold     = V0*(k1.*(1 - q) + k2.*(1 - (q./v)) + k3.*(1 - v));  %Bold
        out.HbO=(out.HbT*100e-6-out.HbR*40e-6)./60e-6;
        out.vascTone=x(5)-1;
        if  M.Model_Choice==3
            out.decay=x(7);
        end
%       
    case {2,7} %Huppert1
        
        %from dissosciation curve
        out.SwO2=a.S ;
        out.SwO2Ratio=(a.S -a.SwO20)/a.SwO20;
        out.Flow =(((a.fin+a.fv)/2)-1); %Flow
        %           Flow = x(3)*((a.fin+a.fv)/2)-1; %Flow
         if isfield(M,'newMeasurement') && M.newMeasurement==1 % better formulation, but not used yet. But almost no diffrence
            out.HbT=(a.Vw0*(x(3)-1 )                           -.5*a.Va0             *(x(2)-1))/(a.Vw0+a.Va0); %hbTratio
            out.HbO=(a.Vw0*(     a.S*x(3) -     a.SwO20 *1  )  -.5*a.Va0 *   a.SaO2  *(x(2)-1))/(a.Vw0*   a.S + a.Va0*   a.SaO2 );
            out.HbR=(a.Vw0*(  (1-a.S)*x(3) - (1-a.SwO20)*1  )  -.5*a.Va0 *(1-a.SaO2) *(x(2)-1))/(a.Vw0*(1-a.S)+ a.Va0*(1-a.SaO2)) ;
         else %old model with fixe arterial volume
             out.HbT=a.Vw0/(a.Vw0+a.Va0)*(x(3)-1); %hbTratio
             out.HbO=a.Vw0/(a.Vw0+a.Va0)*(  a.S*x(3)-a.SwO20*1  )/a.SwO20;
             out.HbR=a.Vw0/(a.Vw0+a.Va0)*(  (1-a.S)*x(3) - (1-a.SwO20)*1  )/(1-a.SwO20);
             
         end
        
        out.HbTConc=a.HbT0*(1+a.phiw0*(x(3)-1));
        
        out.CMRO=(x(5)-1);  %m  CMRO2 signal
        out.PtO2Ratio=(x(6)-1);  %Pression O2 tissus
        out.PtO2=x(6)*a.PtO2;  %Pression O2 tissus
        out.PwO2Ratio=(x(7)-1);  %Pression O2 winkesel
        out.PwO2=x(7)*a.P0;  %Pression O2 winkesel
        out.diffuWT0=a.Psi*( a.Gamma0 -  1) ;  %diffusion of O2 from winkessel to tissus
        out.diffuWT= a.Psi*( a.Gamma0*x(7) -  x(6) );
        out.diffuWTRatio=out.diffuWT/out.diffuWT0-1;
        %autre de Eq44
        out.f1=(a.fin+a.fv-2)*a.delta/(1-a.delta);
        out.a1=out.f1*a.Chi;
        out.a2=-out.f1*a.S;
        out.a3=-out.f1*a.Upsilon*(x(7)-1);
        out.a4=-a.Psi*a.Gamma0*(x(7)-1);
        out.a5=a.Psi*(x(6)-1);
%         out.HbR=out.HbR*1.1; %test
%          out.HbT=out.HbT*1.1; %test
    case 6
        out.Flow=x(2)-1;
end
% scale output 
if isfield(M,'scaleGx') && isstruct(M.scaleGx)
    out.Flow=out.Flow*M.scaleGx.Flow;
    out.HbR=out.HbR*M.scaleGx.HbR;
elseif  isfield(M,'scaleGx') && M.scaleGx==1
     out.Flow=out.Flow*a.scFlow;
    out.HbR=out.HbR*a.scHbR;
else
    %nothing is done %défault and not fitting correctly
end
if M.courbeToFit==99 % reconstruct R V J
    persistent espilonLEff
    if isempty(espilonLEff)
        [void2 espilonLEff]=loadEpsilon( [1 2 3],handAff,0);
    end    
    recon=espilonLEff*[out.HbO*60e-6 out.HbR*40e-6]';
    recon=exp(-recon)-1;
    out.V=recon(1); out.J=recon(2); out.R=recon(3);
    out.Flow=out.Flow/10;
  
end

if 0 %test flow est moins important dans la minimisation
    disp('enlever modif(Flow)')
    Flow=Flow/2;
end
if isfield(M,'output')
    y=[];Y.name={};
    for i1=1:length(M.output)
        try %données HbO
            y(end+1,:)=out.(M.output{i1});
            Y.name{end+1}=M.output{i1};
        catch
            try %données x(1)
                 y(end+1,:)=eval(M.output{i1});
                 Y.name{end+1}=M.output{i1};
            end
        end
    end
    
elseif isfield(M,'YName')
    for i1=1:length(M.YName)
        y(i1)=out.(M.YName{i1});
         Y.name{i1}=M.YName{i1};
    end  
else % determine which variable is present
   y=[]; Y.name=fieldnames(out);
end
y = y(:);
