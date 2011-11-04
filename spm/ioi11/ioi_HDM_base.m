function ioi_HDM_base
% load data
%[Y U M temps tempsStim void1 ]=loadIOI1(M);

% if length(Y.y)>4000
% elseif length(Y.y)>400 %undersampling des output faster processing
%     [Y,U,temps]=undersample(Y,U,temps,400);
% end
% if 0%M.indexJob==313;
%     [Y,U,temps]=undersample(Y,U,temps,100);
% end

% [Y,U,M]=setOptionalSetting(Y,U,M);
% 
% if false % Répéter les données (pour des tests seulement)
%     nn=1;
%     Y.y=repmat(Y.y,nn,1);temps=repmat(temps,nn,1);
% end

% if false % Mettre des confounds  % à retravailler. par défaut pas de confounds
%     Y=doconfound(Y,M.nX0);
% end

% if M.dispDo(1)==1   %apercu des données et des stims
%     figure()
%     U.y=U.u/max(max(U.u))*max(max(Y.y));
%     plot_IOI_SPM(M,temps,Y,temps,U,[],[],1);
%     %     printPDF(gcf,['A:\SPMResume\job' num2str(M.indexJob) '\' M.exp 'ini.png'],'png',[],[8 4])
%     %     return
% end


% if  M.dispDo(2)==1  % test pour voir differentes réponses
%     figure
%     MM=M; MM.output={'HbR','HbT','HbO','CMRO','PwO2Ratio','Flow','PtO2Ratio','diffuWTRatio'};
%     MM=M; MM.output={'HbR','HbT','HbO','CMRO','PwO2Ratio','Flow','PtO2Ratio','diffuWTRatio','x(2)-1','x(5)-1'};
%     MM=M; MM.output={'HbR','HbT','HbO','R','V','J'};
%     MM=M; MM.output={'x(1)-1','x(2)-1' 'x(5)-1','x(6)-1','x(7)-1'};
%     
%     Y33=display22(MM,U,Y,temps,MM.pE);
%     %     figure
%     %     MM=M; MM.output={'a1' 'a2' 'a3' 'a4' 'a5'};
%     %     Y33=display22(MM,U,Y,temps,MM.pE);
%     % return
%     
% end
% 
% 
% if 0 % try to change each prior by 25% and see effect
%     doDistrib(M,[],U,Y,temps) %%try with priors' value
%     return
% end

% if isempty(Y.y)
%     return
% end
% nonlinear system identification
%--------------------------------------------------------------------------
[Ep,Cp,Ce,K0,K1,K2,M0,M1] = spm_nlsi(M,U,Y);
HDM=spm_NLSI_GN_save;

% if 0 % plot graph for report
%     reportArticle
% end


% get values not normalised
[void pE1]=IOI_extraitVar(M,M.pE,ones(M.n,1),ones(M.m,1)); %prior ini
[M.varEst Ep1]=IOI_extraitVar(M,Ep,ones(M.n,1),ones(M.m,1));
Cp1= Cp.* (pE1*pE1');

% if 0 % try to change each prior by 25% and see effect
%     doDistrib(M,Ep,U,Y,temps,pE1) % try with estimated value
% end


if M.dispDo(3)==1  % to see different response on a figure
    figure()
    MM=M; MM.output={'HbR','HbT','HbO','CMRO','PwO2Ratio','Flow','PtO2Ratio','diffuWTRatio'};
    MM=M; MM.output={'x(1)-1','x(2)-1' 'x(5)-1','x(6)'};
    display22(MM,U,Y,temps,Ep)  ; title('from Ep')
end
% 
% if 0 %try different initial state to si if convergent to same resutl
%     testInitialState(M,U,Y,Ep)
% end


if M.dispDo(4)==1% try different prior to see if convergent to same result
    testdifferentpriors(M,U,Y,Ep,Cp,temps)
end

% 
% if 0 % save simulated data pour pouvoir les réutiliser
%     Ysim.y = feval(M.IS,Ep,M,U);
%     loadIOI1(M,0,Ysim,Y) %save simul data if necessary let here
% end

%%-display results
%==========================================================================
Fhdm    = spm_figure;
set(Fhdm ,'name','Hemodynamic Modeling')


%load second set of data if necessary
Y2=Y;

[Y2 U2 M2 temps2 tempsStim void1 ]=loadIOI1(M,-1) ;


if isfield(M,'simulVarName')  %%  Fit new set of data with different parameters
    MM=M;
    %extract parameters from first set of data
    [void1 void2 ind]=intersect(M.simulVarName,MM.name);
    MM.name=MM.name(ind);
    MM.pC= MM.pC(ind,ind)    ;
    MM.pE= ones(size(ind)) ;
    MM.var=MM.varEst;
    [void M.varVctEst]=IOI_extraitVar(M,Ep,ones(MM.n,1));
    [MM ]=t2323(MM,M);
    [Ep3,Cp3,Ce1,K0,K1,K2,M0,M1] = spm_nlsi(MM,U2,Y2);
    [MM.varEst MM.varVctEst]=IOI_extraitVar(MM,Ep3,ones(MM.n,1));
    Ep(ind)=Ep(ind).*MM.varVctEst./M.varVctEst(ind);
    Ep1(ind)=Ep1(ind).*MM.varVctEst./M.varVctEst(ind);
    figure(Fhdm)   ; %retourner bonne figure
end

% display result and simulate
subplot(3,2,4)
title('estimate response','FontSize',9)
Ysim=Y2;
M.d=2;
M=rmfield(M,'d');
U2=scaleStim(U2,Y2,M);

[Ysim.y] = feval(M.IS,Ep,M,U2);
plot_IOI_SPM(M,temps2,Y2,temps2,Ysim,temps2,U2);
set(gca,'xlim',[temps(1) temps(end)])
ylabel('Experimental and estimated '), axis tight


% display extra simulated data
subplot(3,2,6)
title('Simulated response','FontSize',9)
MM=M;
Ysim4=Y2;
MM.output={'vascTone','CMRO','PwO2Ratio','PtO2Ratio','diffuWTRatio','decay'};
MM.output={'x(1)-1' 'x(2)-1' ,'x(5)-1', 'x(6)-1', 'x(7)-1'   'x(8)-1'  'x(9)-1'};
[void,Y4]=IOI_gx_hdm(MM.x,U.u,MM.pE,MM);
MM.output=Y4.name; Ysim4.name=MM.output; MM.l=length(Y4.name);
if length(MM.output)~=0
    [Ysim4.y] = feval(MM.IS,Ep,MM,U2);
    plot_IOI_SPM(M,temps2,Ysim4,[],[],[],[],1);
    set(gca,'xlim',[temps(1) temps(end)])
    ylabel('Simulate '), axis tight
end



% display input parameters
%--------------------------------------------------------------------------
subplot(2,2,1)
if M.mAff>=1
    n1=length(M.pE)-M.mAff+1;
    
    P     = Ep1(n1:end)*1000;
    C     = diag(Cp1(n1:end,n1:end))*1000^2;
    PP     = pE1(n1:end)*1000;
    [~, j] = max(abs(P));
    spm_barh(P,C,PP)
    axis tight
    title({'stimulus efficacy'; 'with 90% confidence intervals'},'FontSize',10)
    set(gca,'Ytick',[1:M.mAff],'YTickLabel',U.name{1},'FontSize',8)
    str = {};
    for i = 1:(M.mAff)
        str{end + 1} = M.name{end-M.mAff+i};
        str{end + 1} = sprintf('mean = %0.2f',P(i));
        str{end + 1} = '';
    end
    set(gca,'Ytick',[1:M.mAff*3]/3 + 1/2,'YTickLabel',str)
    xlabel('relative efficacy per event/sec [*1000]')
end

% display hemodynamic parameters
%---------------------------------------------------------------------------
subplot(2,2,3)
n2=length(M.pE)-M.mAff;
title({ 'hemodynamic parameters'},'FontSize',10)
P     = Ep(1:n2);
pE    = M.pE(1:n2);
C     = diag(Cp(1:n2,1:n2));
spm_barh(P,C,pE)
nomErrorBar(M,P,Ep1)

subplot(3,2,2)
F2=sum(HDM.FProgress);
plot(5:length(F2),F2(5:end),length(F2),HDM.logEv,'o'),title ('Log évidence'),xlabel('Step')


% % display 1 order kernel
doSubplotIOI

% %error calculation
[Ysim.y] = feval(M.IS,Ep,M,U2);
error=mean((Ysim.y-Y2.y(:,1:size(Ysim.y,2))).^2);
% find derivative
if 0
    [dfdpVct f] = spm_diff(M.IS,Ep,M,U,1);
    if iscell(dfdpVct)
        for i1=1:length(dfdpVct);
            dfdp2(i1,:)=mean((dfdpVct{i1}));
            errormat(i1,:)= mean((Y.y-(f+dfdpVct{i1})).^2);
            %     figure(66),plot(dfdpVct{i1}), hold on
        end
    end
else
    dfdp2=0;errormat=0;
end
if 0 % second derivative
    [dfdppVct f] = spm_diff(M.IS,Ep,M,U,[1 1]);
    for i1=1:length(dfdppVct);
        for i2=1:length(dfdppVct{i1});
            dfdpp2(i1,i2)=mean(mean(dfdppVct{i1}{i2}));
        end
    end
else
    dfdpp2=0;
end



%Saving results
% HDM.xY = xY;

M.varVct=pE1;
%find stim amplitude
M=stimArea(M,Ep,U);
HDM.M = M;  HDM.U = U; HDM.Y = Y;
HDM.P = P; HDM.C = C; HDM.EpNorm = Ep;
HDM.Ep=Ep1; HDM.CpNorm = Cp;  HDM.Cp=Cp1;
HDM.Ce = Ce; HDM.K0 = K0; HDM.K1 = K1; HDM.K2 = K2;
HDM.H0 = H0; HDM.H1 = H1; HDM.M0 = M0; HDM.M1 = M1;
HDM=normEp(HDM,U);
HDM.dfdp2=dfdp2;
HDM.dfdpp2=dfdpp2;
HDM.error=mean(error);
HDM.error2=error;
HDM.errormat=errormat;
HDM=orderfields(HDM);
try HDM.ma = ma; catch; end
assignin('base','HDM',HDM)
save('HDM','HDM');
savehdm(HDM,M)

%-Reset title
%--------------------------------------------------------------------------
% spm('FigName',header);
spm('Pointer','Arrow')
% sauver M info
aa=findall(gcf,'label','&SPM Figure');
M.handles=[];
M.U=U;
M.U2=U2;
uimenu('parent',aa(1),'label','display Job','callback',@(hObject,ev)dispFigJob(hObject,ev),'userdata',M)
out = HDM;
out.finish=1;
return

function MenuInfo(hObject,eventdata)
dispFigJob(hObject,eventdata) %laisser pour vieilli figure

function gille2(hObject,eventdata)
dispFigJob(hObject,eventdata)%laisser pour vieilli figure

function [Y,U,temps]=undersample(Y,U,temps,N)
%undersampling des output faster processing
n1=1;
n2=length(Y.y);
dn=round(length(Y.y)/N);
Y.y=Y.y(n1:dn:n2,:);
U.u=U.u(n1:dn:n2,:);
U.dt=Y.dt;
temps=temps(n1:dn:n2);
Y.dt=Y.dt*dn;

function U2=scaleStim(U2,Y2,M)
% scale stim for display
if isfield(M,'source2')
    index= strmatch('U',M.source2);
     U2.y=U2.u;
     if ~isempty(index)       
        U2.y(:,index)=U2.y(:,index)/max(max(U2.y(:,index)))*max(max(Y2.y));   
    end
else
    U2.y=U2.u/max(max(U2.u))*max(max(Y2.y));
end