function ioi_plotSensi(infoGld,infoDtc,vctGld,vctDtc,stim,S,names)
% plot sensitivity graph
%speci (e.g. infoSB):
%vctDtc (e.g. vctB)
%stim
%S
%name (e.g. {'Stim','Bwd'}): names for plot

if isempty (infoDtc)
    gld = infoGld.vct.gld2;
    dtc = infoGld.vct.dtc2;
    bck =zeros(size(dtc)); %background
else
    n=max(length(infoDtc.vct.dtc2),length(infoGld.vct.dtc2));
    bck=zeros(1,n); gld=zeros(1,n);  dtc=zeros(1,n);
    %     bck (infoGld.vct.gld2==0)=1; %background
    bck(1:length(infoDtc.vct.mask)) =infoDtc.vct.mask; %background to remove ( 1=remove, 0 keep)
    
    gld(1:length(infoGld.vct.dtc2))= infoGld.vct.dtc2;
    dtc(1:length(infoDtc.vct.dtc2)) = infoDtc.vct.dtc2;
    enlarge=2;
end
%consider conjunction with gld==0: what is happening outside of
%stimulations
%Considering series A as the truth, to which B is compared
TP = ~bck & gld~=0 & dtc~=0;
TN = ~bck & gld==0 & dtc==0;
FP = ~bck & gld==0 & dtc~=0;
FN = ~bck & gld~=0 & dtc==0;
%output values:
TP1 = sum(TP); TN1 = sum(TN); FP1 = sum(FP); FN1 = sum(FN);
sensitivity=TP1/(TP1+FN1); specificity=1-FP1/(FP1+TN1);
posPred=TP1/(TP1+FP1); negPred=TN1/(TN1+FN1);
%disp(['TP: ' int2str(TP1) ' FP: ' int2str(FP1) ' FN: ' int2str(FN1) ' TN: ' int2str(TN1)]);
%disp(['Sensitivity: ' num2str(sensitivity,'%1.2f') ' Specificity: ' ...
%    num2str(specificity,'%1.2f') ' Positive Predictive value: ' num2str(posPred,'%1.2f') ...
%    ' Negative Predictive value: ' num2str(negPred,'%1.2f')])

%plot data series
offset=S.offset;
temps=(1:length(vctDtc))/S.freq+offset;
h0 = plot(temps,vctDtc,'displayName',names{2}); axis tight
if isempty (infoDtc)
    ylabel([names{2}]);
else
    ylabel([names{2} ' against ' names{1}])
    hold on
    temps3=(1:length(vctGld))/S.freq+offset;
    if enlarge~=0
        enlarge= [enlarge/(max(vctGld)-min(vctGld)) -1];
        yy=get(gca,'ylim');
        enlarge= [diff(yy)/(max(vctGld)-min(vctGld)) diff(yy) ];
        h00 = plot(temps3,vctGld*enlarge(1)-enlarge(2),'k.-', 'displayname',names{1} ); axis tight
    else
        h00 = plot(temps3,vctGld,'k.-', 'displayname',names{1} );
    end
    
end
if S.set_limits
    set(gca,'xlim',S.xlim)
end
xlabel('Time (s)')
hold on
%get vertical limits, to use to set TP, FP, FN, TN labels and stim bars
yy=get(gca,'ylim');
in_vpos = yy(1)-0.05*diff(yy);

%set 4s windows boxes
% temps2 = S.dtwindow*(1:length(gld))-S.dtwindow/2; %center of the window
temps2 = S.dtwindow*(1:length(gld))-S.dtwindow/2+S.offset; %center of the window
%plot(temps2(TP),in_vpos,'bo',temps2(FP),in_vpos,'rs',temps2(FN),in_vpos,'gd',temps2(TN),in_vpos,'kx')
try h1 = plot(temps2(TP),temps2(TP)*0+in_vpos','bo','displayName','TP'); catch, h1 = plot(temps2,temps2*0+in_vpos,'','displayName','TP'); end
try h2 = plot(temps2(FP),temps2(FP)*0+in_vpos,'rs','displayName','FP'); catch, h2 = plot(temps2,temps2*0+in_vpos,'','displayName','FP'); end
try h3 = plot(temps2(FN),temps2(FN)*0+in_vpos,'gd','displayName','FN'); catch, h3 = plot(temps2,temps2*0+in_vpos,'','displayName','FN'); end
try h4 = plot(temps2(TN),temps2(TN)*0+in_vpos,'kx','displayName','TN'); catch, h4 = plot(temps2,temps2*0+in_vpos,'','displayName','TN'); end
try
    in_vpos2 = yy(2)+0.05*diff(yy);
    t=offset+ stim.ind/S.freq;
    %t=offset + S.stimInd/S.freq; %center of the window - there is a problem in
    %the timing of these measures
    hs=plot(t,t*0+in_vpos2,'k.','markersize',10,'displayname','Stim'); %ones(size(t))*
end
set(gca,'ylim',[yy(1)-0.1*diff(yy) yy(2)+0.1*diff(yy)])
if  isempty (infoDtc)
    if ~isempty(infoGld.s_th_ampli)
        hTh = plot(temps2,infoGld.s_th_ampli*ones(size(temps2)),'k-','displayname','Thld'); % Threshold line
        title(['Detection for a threshold corresponding to a sensitivity of ' num2str(infoGld.s_th,'%1.2f') ])
    end
elseif nargin == 7
    
    if ~isempty(infoDtc.s_th_ampli)
        hTh = plot(temps2,infoDtc.s_th_ampli*ones(size(temps2)),'b-'); %
    end
    if ~isempty(infoGld.s_th_ampli)
        if length(enlarge)==2
            hTh2 = plot(temps2,(infoGld.s_th_ampli*ones(size(temps2))) *enlarge(1)-enlarge(2),'k.-','displayname',['Thld ' names{1} num2str(infoGld.s_th,'%1.2f')  ]); % Threshold line
        else
            hTh2 = plot(temps2,infoGld.s_th_ampli*ones(size(temps2)),'k.-','displayname',['Thld ' names{1} num2str(infoGld.s_th,'%1.2f')  ]); % Threshold line
        end
    end
    %     legend([h0 h00 h1(1) h2(1) h3(1) h4(1) hTh(1) hTh2(1) hs],[names{2} names{1} {'TP','FP','FN','TN',['Thld ' names{2} num2str(infoDtc.s_th,'%1.2f') ] ['Thld ' names{1}
    
    title(['Detection for a threshold corresponding to a given sensitivity  ' ])
    
end
legend('show')