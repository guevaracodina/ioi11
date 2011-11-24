function y = ioi_gx_SCKS(x,u,PS)
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

% exponentiation of hemodynamic state variables
%--------------------------------------------------------------------------
% if isstruct(PS)
P = PS.pE(:,1); % .* PS.pEreal;
% else
%     P = PS;
%     clear PS;
%     PS.PhysioModel_Choice = 0;
%     PS.xY.includeHbR = 1;
%     PS.xY.includeHbT = 1;
%     PS.xY.includeFlow = 1;
% end
x     = exp(x);
switch PS.PhysioModel_Choice
    case 0 %Buxton, 
        v     = x(3,:);        
        q     = x(4,:);
        fv       = (x(3,:).^(1/P(4)));%Flow Out      
        out.HbT = (v-1).*1; %Assume v = HbT = HbO+ HbR
        out.HbR = (q-1).*1; %Assume q = HbR 
        out.Flow= ( (x(2,:)+fv)/2-1); % speckle=mean (flow in , flow out)
        %out.HbO=(out.HbT*100e-6-out.HbR*40e-6)./60e-6;            
    case {1,3} %, Zheng Mayhew and Zheng Decay
        v     = x(3);        q     = x(4);
        fv       = (x(3)^(1/P(4)))/x(5);%Flow Out       
        out.HbT = (v-1); %Assume v = HbT = HbO+ HbR
        out.HbR = (q-1); %Assume q = HbR
        out.Flow= ( (x(2)+fv)/2-1); % speckle=mean (flow in , flow out)
        out.HbO=(out.HbT*100e-6-out.HbR*40e-6)./60e-6;
        out.vascTone=x(5)-1;
        if  M.Model_Choice==3
            out.decay=x(7);
        end
    case {2,7} %Huppert1        
        %from dissociation curve
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
end
% % % scale output 
% % if isfield(M,'scaleGx') && isstruct(M.scaleGx)
% %     out.Flow=out.Flow*M.scaleGx.Flow;
% %     out.HbR=out.HbR*M.scaleGx.HbR;
% % elseif  isfield(M,'scaleGx') && M.scaleGx==1
% %      out.Flow=out.Flow*a.scFlow;
% %     out.HbR=out.HbR*a.scHbR;
% % else
% %     %nothing is done %défault and not fitting correctly
% % end
y = [];
if PS.xY.includeHbR, y = [y; out.HbR]; end
if PS.xY.includeHbT, y = [y; out.HbT]; end   
if PS.xY.includeFlow, y = [y; out.Flow]; end
%y = y(:);