function y = ioi_gx(x,u,P,M)
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
%P = PS.pE(:,1); % .* PS.pEreal;
% else
%     P = PS;
%     clear PS;
%     PS.PhysioModel_Choice = 0;
%     PS.xY.includeHbR = 1;
%     PS.xY.includeHbT = 1;
%     PS.xY.includeFlow = 1;
% end
x     = exp(x);
switch M.PS.PhysioModel_Choice
    case 0 %Buxton, 
        v     = x(3,:);        
        q     = x(4,:);
        fv       = (x(3,:).^(1/P(4)));%Flow Out      
        out.HbT = (v-1).*1; %Assume v = HbT = HbO+ HbR
        out.HbR = (q-1).*1; %Assume q = HbR 
        out.Flow= ( (x(2,:)+fv)/2-1); % speckle=mean (flow in , flow out)
        %out.HbO=(out.HbT*100e-6-out.HbR*40e-6)./60e-6;            
    case 1 %, Zheng Mayhew and Zheng Decay
        v     = x(3,:);        
        q     = x(4,:);
        fv       = (x(3,:)^(1/P(4)))/x(5,:);%Flow Out       
        out.HbT = (v-1); %Assume v = HbT = HbO+ HbR
        out.HbR = (q-1); %Assume q = HbR
        out.Flow= ( (x(2,:)+fv)/2-1); % speckle=mean (flow in , flow out)
        %out.HbO=(out.HbT*100e-6-out.HbR*40e-6)./60e-6;
        out.vascTone=x(5,:)-1;
    case 2 %Huppert1        
        BH = M.PS.BH;
        BH = ioi_fillBH(x,BH);
        %from dissociation curve
% % %         out.SwO2=BH.S ;
% % %         out.SwO2Ratio=(BH.S -BH.SwO20)/BH.SwO20;
        out.Flow =(((BH.fin+BH.fv)/2)-1); %Flow
        %           Flow = x(3)*((BH.fin+BH.fv)/2)-1; %Flow
%          if isfield(M,'newMeasurement') && M.newMeasurement==1 % better formulation, but not used yet. But almost no diffrence
            out.HbT=(BH.Vw0*(x(3,:)-1 )                           -.5*BH.Va0             *(x(2,:)-1))/(BH.Vw0+BH.Va0); %hbTratio
            %out.HbO=(BH.Vw0*(     BH.S*x(3) -     BH.SwO20 *1  )  -.5*BH.Va0 *   BH.SaO2  *(x(2)-1))/(BH.Vw0*   BH.S + BH.Va0*   BH.SaO2 );
            out.HbR=(BH.Vw0*(  (1-BH.S)*x(3,:) - (1-BH.SwO20)*1  )  -.5*BH.Va0 *(1-BH.SaO2) *(x(2,:)-1))/(BH.Vw0*(1-BH.S)+ BH.Va0*(1-BH.SaO2)) ;
%          else %old model with fixe arterial volume
%              out.HbT=BH.Vw0/(BH.Vw0+BH.Va0)*(x(3)-1); %hbTratio
%              %out.HbO=BH.Vw0/(BH.Vw0+BH.Va0)*(  BH.S*x(3)-BH.SwO20*1  )/BH.SwO20;
%              out.HbR=BH.Vw0/(BH.Vw0+BH.Va0)*(  (1-BH.S)*x(3) - (1-BH.SwO20)*1  )/(1-BH.SwO20);
%              
%          end
        
% % %         out.HbTConc=BH.HbT0*(1+BH.phiw0*(x(3,:)-1));
% % %         out.CMRO=(x(5,:)-1);  %m  CMRO2 signal
% % %         out.PtO2Ratio=(x(6,:)-1);  %Pression O2 tissus
% % %         out.PtO2=x(6,:)*BH.PtO2;  %Pression O2 tissus
% % %         out.PwO2Ratio=(x(7,:)-1);  %Pression O2 winkesel
% % %         out.PwO2=x(7,:)*BH.P0;  %Pression O2 winkesel
% % %         out.diffuWT0=BH.Psi*( BH.Gamma0 -  1) ;  %diffusion of O2 from winkessel to tissus
% % %         out.diffuWT= BH.Psi*( BH.Gamma0*x(7,:) -  x(6,:) );
% % %         out.diffuWTRatio=out.diffuWT/out.diffuWT0-1;
% % %         %autre de Eq44
% % %         out.f1=(BH.fin+BH.fv-2)*BH.delta/(1-BH.delta);
% % %         out.a1=out.f1*BH.Chi;
% % %         out.a2=-out.f1*BH.S;
% % %         out.a3=-out.f1*BH.Upsilon*(x(7,:)-1);
% % %         out.a4=-BH.Psi*BH.Gamma0*(x(7,:)-1);
% % %         out.a5=BH.Psi*(x(6,:)-1);
end
% % % scale output 
% % if isfield(M,'scaleGx') && isstruct(M.scaleGx)
% %     out.Flow=out.Flow*M.scaleGx.Flow;
% %     out.HbR=out.HbR*M.scaleGx.HbR;
% % elseif  isfield(M,'scaleGx') && M.scaleGx==1
% %      out.Flow=out.Flow*BH.scFlow;
% %     out.HbR=out.HbR*BH.scHbR;
% % else
% %     %nothing is done %défault and not fitting correctly
% % end
y = [];
if M.PS.xY.includeHbR, y = [y; out.HbR]; end
if M.PS.xY.includeHbT, y = [y; out.HbT]; end   
if M.PS.xY.includeFlow, y = [y; out.Flow]; end
%y = y(:);