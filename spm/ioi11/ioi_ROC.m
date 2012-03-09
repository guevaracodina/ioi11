function ROC=ioi_ROC(SCKS,Options)
% to execute : sensitivitySCKS(SCKS,physio,info1)
%S: ROC structureaf
%ROCmethod:
%2: Set threshold on amplitude
%3: Set threshold on area
S.ROC_method = 4; %this is set for the frequency of 0.2 Hz. So changing this mode has no effect
S.UseFullStims = 1; %Boolean: 0: true positives only for beginning of stimulations; 1: true positives for whole duration of stimulation
S.ShiftedPower = 0; %Boolean: 0: no change; 1: shifts tested series to positive values, then take square to get power
S.offset=0;
S.shift = 0; %0.33;
S.set_limits = 1; %Boolean: 1: use S.xlim for x-axis limits, 0: display whole series
S.xlim=[1 200]; %Display limits
S.xlim=[1 2000]; %Display limits
S.decim=3; %Low pass filter frequency, in Hz
S.freq=5; %Acquisition frequency in Hz
S.dtwindow=4; %detection window length in seconds  % need to be a multiple of the acquisition time 0.2s
S.s_th = 0.85; %0.95; %Sensitivity threshold
S.sk_npts2 = 1; %becomes an offset??? - point 1 is time zero
S.sk_npts = 10; %5; %number of points set to zero at beginning of backward pass SCKS deconvolved inputs
S.dtFreq=round(S.dtwindow*S.freq);
S.winpts = round(S.dtwindow*S.freq);

%FIND STIM index
stim0 = SCKS.pU.v{1};
time0 = (1:length(stim0))*SCKS.dt;
S.stimInd=find(stim0>0.05)'; 

namea={'SCKS.qU.v{1}''' 'SCKS.qU.v{2}'''  'SCKS.qU.C{1}''' };
nameb={'vctB' 'vctF' 'vctN'};
% form segment sv
if isfield(SCKS,'sv') %if SCKS done with a part of all signal
    for i1=1:length(namea)
        temp=eval(namea{i1});
        temp(1:S.sk_npts)=0;
        r.(nameb{i1})=temp;
    end
    vctB=SCKS.qU.v{1}'; vctB(1:S.sk_npts)=0; %SCKS (backward pass) deconvolved inputs; first 5 data points set to zero due to spurious effects
    vctF=SCKS.qU.v{2}'; vctF(1:S.sk_npts)=0; %SCKS (backward pass) deconvolved inputs; first 5 data points set to zero due to spurious effects
    vctN=SCKS.qU.C{1}'; vctN(1:S.sk_npts)=0; %SCKS (ba
    S.offset=(SCKS.sv(1)-1)/S.freq;
else %if SCKS done with all signal (sv is not present)
    se=min(handAff([],'display.tFin')*S.freq,200000);
    si = max(handAff([],'display.tIni')*S.freq,50);
    se = min(se,size(SCKS.pU.r,2));se = min(se,length(time0));
    sv = si:se; %start and end indices - ignore first few data points
    SCKS.pU.r = SCKS.pU.r(:,sv); %Modalities: HbR, HbT, Flow
    SCKS.sv = sv;
    for i1=1:length(namea)
        temp=eval(namea{i1});
        temp(1:S.sk_npts)=0;
        r.(nameb{i1})=temp(SCKS.sv);
    end
    S.offset=(si-1)/S.freq;
end



stim0=stim0(SCKS.sv); time0=time0(SCKS.sv);
S.stimInd=S.stimInd-si+1;
S.xlim=[min(time0) max(time0)];
S.n = size(SCKS.pU.r,2);
S.s_th = 0.5; %Sensitivity threshold for electro

%Select segment used


%load data
if isfield(physio,'fftHigh')
    %         1        2        3           4         5           6     7        8      9      10        11
    name={'temps' 'fftHigh' 'fftLow' 'fftLFPLow' 'fftLFPHigh' 'MUA4' 'MUA4b' 'MUA6' 'MUA4b' 'fftLow' 'fftLFPHighPw' 'fftLFPHigh'};
    name2={'SCKS' 'MUAini'     'LFP'     'LFPLow'     'LFPHigh' 'MUA4' 'MUA4b' 'MUA' 'MUA4b' 'LFPF' 'LFPPw ' 'pseudorandom'};
else
    name={'temps' };
    name2={'SCKS' };
end
% vctH =extPhysio(physio,'fftHigh'); %Multiunit action potential
for i1=1:length(name);
    q(i1).name=name{i1};
    q(i1).name2=name2{i1};
    temp=extPhysio(physio,name{i1}); %Local field potential
    %     temp=max(min(std(temp)*10,temp),-std(temp)*10);
    %     temp=max(min(std(temp)*10,temp),-std(temp)*10);
    q(i1).vct=temp(SCKS.sv);
end

% 5create pseudoramdom vec
for i1=1:size(info1.stim1,1)-1
    try
        ind=find((time0>(info1.stim1(i1,3) -20)& time0<(info1.stim1(i1,3) )));
        ind2=find((time0>(info1.stim1(i1+1,3) -20)& time0<(info1.stim1(i1+1,3) )));
        q(12).vct(ind(1:100))=q(5).vct(ind2(1:100));
    catch
        1;
    end
end

if Options.filter % low pass     % to have same treatment as filter done in SCKS_run
    for i1=2:length(q);
        q(i1).vct=filtreHighSCKS( q(i1).vct);
    end
end
if 1 % low pass
    for i1=1:length(q)
        q(i1).vct  =filtredecim(q(i1).vct,S.decim);
    end
    r.vctF  =filtredecim(r.vctF,S.decim);
    %     vctH=filtredecim(vctH,S.decim);
    %     vctL=filtredecim(vctL,S.decim);
end
if 1 % high pass - PP
    [b1,a1] = butter(2,0.005,'high');
    q(1).vct=filtfilt(b1,a1,q(1).vct);
    r.vctF=filtfilt(b1,a1,r.vctF);
    
end
q(1).vct=r.vctB;
if 1 %highpass 2
    [b2,a2] = butter(4,.2/(5/2),'high');
    q(end).vct=filtfilt(b2,a2,q(end).vct);
end

%downsample
for i1=1:length(q)
    [q(i1).amp  q(i1).ind q(i1).indMax] =formResample(S,q(i1).vct);
end


%Find stimulation times
if 0 %old
    stim=findStim(r.vctS,S);
else
    stim.ind= find(stim0>0.05)';
    stim.amp= ones(size(stim.ind));
    % stim stimulation time of stim 1 one second stimulation
    aa=round(info1.stim1(info1.stim1(:,7)==1,3)*5-50)';
    bb=-20:1:20;
    bb=(repmat(aa',1,size(bb,2)) +repmat(bb,size(aa,2),1) );
    stim1.ind=intersect(stim.ind,bb(1:end));
    stim1.amp= ones(size(stim1.ind));
    aa=round(info1.stim1(info1.stim1(:,7)==2,3)*5-50)';
    bb=-20:1:20;
    bb=(repmat(aa',1,size(bb,2)) +repmat(bb,size(aa,2),1) );
    stim2.ind=intersect(stim.ind,bb(1:end));
    stim2.amp= ones(size(stim2.ind));
    aa=round(info1.stim1(info1.stim1(:,7)==3,3)*5-50)';
    bb=-20:1:20;
    bb=(repmat(aa',1,size(bb,2)) +repmat(bb,size(aa,2),1) );
    stim3.ind=intersect(stim.ind,bb(1:end));
    stim3.amp= ones(size(stim3.ind));
    % stim stimulation time of stim 5 and 15  second stimulation
    aa=round(info1.stim1((info1.stim1(:,7)==2 |info1.stim1(:,7)==3)  ,3)*5-50)';
    bb=-10:1:85;
    bb=(repmat(aa',1,size(bb,2)) +repmat(bb,size(aa,2),1) );
    stim515s.ind=intersect(stim.ind,bb(1:end));
    stim515s.amp= ones(size(stim515s.ind));
    
end

[SCKSResult.EQM SCKSResult.EQM2  ] =StatsSCKS(SCKS,[]); %find quality of fit

SCKSResult.S=S;
SCKSResult=findEQM(q(1:end-1),q(end).vct,SCKSResult,time0);

if 0 %plot des fit
    U.u=q(1).vct;    U.dt=median(diff(time0));
    U2.u=filtfilt(b1,a1,q(4).vct);    U2.dt=median(diff(time0));
    [Ysim.y] = feval(SCKS.M(1).IS,ones(6,1),SCKS.M(1),U);
    [Ysim2.y] = feval(SCKS.M(1).IS,ones(6,1),SCKS.M(1),U2);
    subplot(2,1,1);plot(time0,SCKS.pU.r')
    subplot(2,1,2);plot(time0,Ysim.y)
    subplot(2,1,2);plot(time0,Ysim2.y)
end




q2=q;
offset1=3; % time to enlarge stim to create mask
offset1=round(offset1/median(diff(time0)));
offset1=-offset1:1:offset1;
stimAll.ind=unique(repmat(stim.ind',1,size(offset1,2)) +repmat(offset1,size(stim.ind,2),1) );
stimAll.amp= ones(size( stimAll.ind));
stim3515s.ind=unique(repmat(stim515s.ind,1,size(offset1,2)) +repmat(offset1,1,size(stim515s.ind,2)) );
stim3515s.amp= ones(size( stim3515s.ind));
[mask1 mask2 mask3 mask1NS mask2NS mask3NS]=maskSensi(time0,info1,q(1),stimAll);
%stimAll : remove stim with an offset +-3 s
%stim3515s : remove stim 3 5 with an offset +-3 s
%Mask1 : take just the protocole 1
%Mask1 : take just the protocole 1 without the stim
%Find ROC goldstandard = stim
for i1= 1:length(q)
    q(i1).Stim    =roc(stim,q(i1),[],S);
    q(i1).Stim1s    =roc(stim1,q(i1),stim3515s,S); % only stim of 1 sec %old
    q(i1).Stim1    =roc(stim1,q(i1),mask1,S); % only stim of 1
    q(i1).Stim2    =roc(stim2,q(i1),mask2,S); % only stim of 2
    q(i1).Stim3    =roc(stim3,q(i1),mask3,S); % only stim of 3
    q(i1).ROCNoStim =roc(stim,q(i1),stimAll,S); %
    q(i1).ROCNoStim1 =roc(stim,q(i1),mask1NS,S); %
    q(i1).ROCNoStim2 =roc(stim,q(i1),mask2NS,S); %
    q(i1).ROCNoStim3 =roc(stim,q(i1),mask3NS,S); %
    q(i1).gldStim  =extractDtc(stim,q(i1),[],S,q(i1).Stim.roc,'sensi',.85);
    q(i1).gldStimB  =extractDtc(stim,q(i1),[],S,q(i1).ROCNoStim.roc,'speci',.5); %
    q(i1).gldStimCNoStim  =extractDtc(stim,q(i1),stimAll,S,[],'std',1); %
    q(i1).gldStimCNoStim1  =extractDtc(stim1,q(i1),stim1,S,[],'std',1); %
    q(i1).gldStimCNoStim2  =extractDtc(stim2,q(i1),stim2,S,[],'std',1); %
    q(i1).gldStimCNoStim3  =extractDtc(stim3,q(i1),stim3,S,[],'std',1); %
%     q(i1).gldStimDNoStim  =extractDtc(stim,q(i1),stimAll,S,[],'std',2); %
end


% find threshold electro to have sufficient number of stim outside
% stimulation period


noMask=stim;noMask.amp=noMask.amp(1); noMask.ind=noMask.ind(1);
noMask=stim;noMask.amp=[]; noMask.ind=[];
for i1=1:length(q)
    q(i1).elecNoStim     =roc(q(i1),q(1),stim,S, q(i1).gldStimB.s_th_ampli);
    q(i1).gldElecSNoStim  =extractDtc(q(i1),q(1),stim,S,q(i1).elecNoStim.roc,'sensi',.85);
    q(i1).elecAll=roc(q(i1),q(1),noMask,S,q(i1).gldStimB.s_th_ampli);
    q(i1).gldElecSAll  =extractDtc(q(i1),q(1),noMask,S,q(i1).elecAll.roc,'sensi',.85);
    
    q(i1).elecAllC=roc(q(i1),q(1),noMask,S,q(i1).gldStimCNoStim.s_th_ampli);
    q(i1).gldElecSAllC  =extractDtc(q(i1),q(1),noMask,S,q(i1).elecAllC.roc,'sensi',.85);
    q(i1).elecAll1C=roc(q(i1),q(1),mask1,S,q(i1).gldStimCNoStim.s_th_ampli);
    q(i1).gldElecSAll1C  =extractDtc(q(i1),q(1),mask1,S,q(i1).elecAll1C.roc,'sensi',.85);
    q(i1).elecAll2C=roc(q(i1),q(1),mask2,S,q(i1).gldStimCNoStim.s_th_ampli);
    q(i1).gldElecSAll2C  =extractDtc(q(i1),q(1),mask2,S,q(i1).elecAll2C.roc,'sensi',.85);
    q(i1).elecAll3C=roc(q(i1),q(1),mask3,S,q(i1).gldStimCNoStim.s_th_ampli);
    q(i1).gldElecSAll3C  =extractDtc(q(i1),q(1),mask3,S,q(i1).elecAll3C.roc,'sensi',.85);
    
    q(i1).elecNoStimB     =roc(q(i1),q(1),stimAll,S, q(i1).gldStimB.s_th_ampli);
    q(i1).gldElecSNoStimB  =extractDtc(q(i1),q(1),stimAll.ind,S,q(i1).elecNoStimB.roc,'sensi',.85);
    q(i1).elecNoStimC     =roc(q(i1),q(1),stimAll,S,  q(i1).gldStimCNoStim.s_th_ampli);
    q(i1).gldElecSNoStimC  =extractDtc(q(i1),q(1),stimAll.ind,S,q(i1).elecNoStimC.roc,'sensi',.85);
    q(i1).elecNoStim1C     =roc(q(i1),q(1),mask1NS,S,  q(i1).gldStimCNoStim.s_th_ampli);
    q(i1).gldElecSNoStim1C  =extractDtc(q(i1),q(1),mask1NS.ind,S,q(i1).elecNoStimC.roc,'sensi',.85);
    q(i1).elecNoStim2C     =roc(q(i1),q(1),mask2NS,S,  q(i1).gldStimCNoStim.s_th_ampli);
    q(i1).gldElecSNoStim2C  =extractDtc(q(i1),q(1),mask2NS.ind,S,q(i1).elecNoStimC.roc,'sensi',.85);
    q(i1).elecNoStim3C     =roc(q(i1),q(1),mask3NS,S,  q(i1).gldStimCNoStim.s_th_ampli);
    q(i1).gldElecSNoStim3C  =extractDtc(q(i1),q(1),mask3NS.ind,S,q(i1).elecNoStimC.roc,'sensi',.85);
    
    %     q(i1).elecNoStimD     =roc(q(i1),q(1),stimAll,S,  q(i1).gldStimDNoStim.s_th_ampli);
    %     q(i1).gldElecSNoStimD  =extractDtc(q(i1),q(1),stimAll.ind,S,q(i1).elecNoStimD.roc,'sensi',.85);
    
    q(i1).SCKSNoStim     =roc(q(1),q(i1),stim,S,q(1).gldStimB.s_th_ampli);
    q(i1).gldSCKSSNoStim  =extractDtc(q(1),q(i1),[],S,q(1).SCKSNoStim.roc,'sensi',.85);
    q(i1).SCKSNoStimC     =roc(q(1),q(i1),stimAll,S,  q(1).gldStimCNoStim.s_th_ampli);
    q(i1).gldSCKSSNoStimC  =extractDtc(q(1),q(i1),stimAll.ind,S,q(1).SCKSNoStimC.roc,'sensi',.85);
%     q(i1).SCKSNoStimD     =roc(q(1),q(i1),stimAll,S,  q(1).gldStimDNoStim.s_th_ampli);
%     q(i1).gldSCKSSNoStimD  =extractDtc(q(1),q(i1),stimAll.ind,S,q(1).SCKSNoStimD.roc,'sensi',.85);
%     
    
    q(i1).SCKSAll=roc(q(1),q(i1),noMask,S,q(1).gldStimB.s_th_ampli);
    q(i1).gldSCKSSAll  =extractDtc(q(1),q(i1),[],S,q(1).SCKSAll.roc,'sensi',.85);
    
end
try,disp( q(5).elecNoStimD) ,disp( q(end).elecNoStimD), end
% find the threshold on stim that is the best
if Options.thld && Options.Plot==0
    ThldSensi;
else
    thld5=[];
end




SCKSResult.roc=reduce1(q);
q(1).time0=time0;

%snr of the input that match the stim
ind1=50; % index to remove initial time that is not really good fit
index=q(1).gldStim.vct.dtc2;
index=q(1).indMax(logical(index));
index=index(index>ind1);
snr=r.vctB(index)./abs(r.vctN(index));
SCKSResult.snr.snrU=mean(snr);
SCKSResult.snr.meanU=mean(r.vctB(index));
SCKSResult.snr.stdB=std(r.vctB(ind1:end));
SCKSResult.snr.snrBU=mean(r.vctB(index))/std(r.vctB(ind1:end));
try
    SCKSResult.thld=thld5;
end
if Options.Plot==0
    return
end

%plot ROC where goldStandard is electric stim
figure(1)
trait={'B-','G.-','R:','K--','M*',':'};
index=[   4 5 8  1] ;
for i0=1:length(index)
    i1=index(i0);
    leg=[q(i1).name2 ', area = ' num2str(q(i1).Stim.ROCarea,'%1.2f')];
    aa=plot(q(i1).Stim.roc(:,2), q(i1).Stim.roc(:,1), trait{i0},'displayname', leg);
    leg=[q(i1).name2 ', area = ' num2str(q(i1).Stim1.ROCarea,'%1.2f')];
    aa=plot(q(i1).Stim1.roc(:,2), q(i1).Stim1.roc(:,1), trait{i0},'displayname', leg);
    
    hold on
end
hold off;
a=legend('show');
set(a,'location','southeast')
title('ROC, for electric stimulations'),set(gcf,'tag','ROC')
xlabel('1-Specificity'), ylabel('Sensitivity')
figure(1111)
trait={'B-','G.-','R:','K--','M*',':'};
index=[   4 5 8  1] ;
for i0=1:length(index)
    i1=index(i0);
    leg=[q(i1).name2 ', area = ' num2str(q(i1).Stim1.ROCarea,'%1.2f')];
    aa=plot(q(i1).Stim1.roc(:,2), q(i1).Stim1.roc(:,1), trait{i0},'displayname', leg);
    
    hold on
end
hold off;
a=legend('show');
set(a,'location','southeast')
title('ROC, for electric stimulations'),set(gcf,'tag','ROC')
xlabel('1-Specificity'), ylabel('Sensitivity')
%plot stim

figure(2);
index=[1 5 8];

for i0=1:length(index)
    i1=index(i0);
    subplot(3,1,i0);
    plotSensi(q(i1).gldStim,[],[],q(i1).vct,stim,S,{'Stim',q(i1).name2});
end

subplot(3,1,1);xlabel([]); title([]);a=legend('hide');set(a,'location','eastoutside'); set(gca,'xTickMode', 'manual','xtick',[],'yTickMode', 'manual','ytick',[])
set(gca,'position',[.13 .71 .6 .25])
subplot(3,1,2);xlabel([]); title([]);a=legend('hide');set(a,'location','eastoutside'); set(gca,'xTickMode', 'manual','xtick',[],'yTickMode', 'manual','ytick',[])
set(gca,'position',[.13 .43 .6 .25])
subplot(3,1,3); title([]);a=legend('show'); set(a,'location','eastoutside'); set(gca,'yTickMode', 'manual','ytick',[])
set(gca,'position',[.13 .15 .6 .25])

%Consider only periods without stimulations, and compare SCKS (series 2)
%against LFP (series 1)
figure(5);plotSensi( q(2).gldStimB , q(1).gldStimB ,q(2).vct,q(1).vct,stim,S,{'LFP','SCKS'});
title('SCKS against LFP outside stimulation period')
%ROC, also for periods without stimulations



coul={'B','G','R' 'K','Y'};
figure(6)
index=[  3 4 5 8] ;
for i0=1:length(index)
    i1=index(i0);
    plot(q(i1).elecNoStim.roc(:,2),q(i1).elecNoStim.roc(:,1),coul{i0},'displayname',['SCKS ' q(i1).name2 ', area = ' num2str(q(i1).elecNoStim.ROCarea,'%1.2f')])
    hold on
    plot(q(i1).elecNoStimB.roc(:,2),q(i1).elecNoStimB.roc(:,1),[coul{i0},'--'],'displayname',['SCKS ' q(i1).name2 'Out, area = ' num2str(q(i1).elecNoStimB.ROCarea,'%1.2f')])
    plot(q(i1).elecAll.roc(:,2),q(i1).elecAll.roc(:,1),[coul{i0} ':'],'displayname',['SCKS ' q(i1).name2 'All, area = ' num2str(q(i1).elecAll.ROCarea,'%1.2f')])
    plot(q(i1).elecNoStimC.roc(:,2),q(i1).elecNoStimC.roc(:,1),[coul{i0},'-.'],'displayname',['SCKS ' q(i1).name2 'STDout, area = ' num2str(q(i1).elecNoStimC.ROCarea,'%1.2f')])
    %       plot(q(i1).elecNoStimD.roc(:,2),q(i1).elecNoStimD.roc(:,1),[coul{i0},'-.'],'displayname',['SCKS ' q(i1).name2 'STD, area = ' num2str(q(i1).elecNoStimC.ROCarea,'%1.2f')])
    plot(q(i1).elecAllC.roc(:,2),q(i1).elecAllC.roc(:,1),[coul{i0} '-x'],'displayname',['SCKS ' q(i1).name2 'AllC, area = ' num2str(q(i1).elecAllC.ROCarea,'%1.2f')])
    
end
hold off
a=legend('show');
set(a,'location','southeast')
title('ROC, SCKS against LFP,MUA outside stimulation period'),set(gcf,'tag','ROC')
text(.5,.3,[ 'RSB reconstruction '  num2str(SCKSResult.EQM.ratiooz,'%2.2f  ')])
xlabel('1-Specificity'), ylabel('Sensitivity')




figure(66)
index=[  3 4 5 8] ;
for i0=1:length(index)
    i1=index(i0);
    plot(q(i1).SCKSNoStim.roc(:,2),q(i1).SCKSNoStim.roc(:,1),coul{i0},'displayname',['SCKSNoStim ' q(i1).name2 ', area = ' num2str(q(i1).SCKSNoStim.ROCarea,'%1.2f')])
    hold on
    plot(q(i1).SCKSAll.roc(:,2),q(i1).SCKSAll.roc(:,1),[coul{i0} ':'],'displayname',['SCKS ' q(i1).name2 'All, area = ' num2str(q(i1).SCKSAll.ROCarea,'%1.2f')])
end
hold off
a=legend('show');
set(a,'location','southeast')
title('ROC, SCKS against LFP,MUA outside stimulation period'),set(gcf,'tag','ROC')
text(.5,.3,[ 'RSB reconstruction '  num2str(SCKSResult.EQM.ratiooz,'%2.2f  ')])
xlabel('1-Specificity'), ylabel('Sensitivity')




figure(7)
index=[  7 8 9 10] ;
for i0=1:length(index)
    i1=index(i0);
    plot(q(i1).elecNoStim.roc(:,2),q(i1).elecNoStim.roc(:,1),coul{i0},'displayname',['SCKS ' q(i1).name2 ', area = ' num2str(q(i1).elecNoStim.ROCarea,'%1.2f')])
    hold on
    plot(q(i1).elecAll.roc(:,2),q(i1).elecAll.roc(:,1),[coul{i0} ':'],'displayname',['SCKS ' q(i1).name2 'All, area = ' num2str(q(i1).elecAll.ROCarea,'%1.2f')])
end
hold off
a=legend('show');
set(a,'location','southeast')
title('ROC, SCKS against LFP,MUA outside stimulation period'),set(gcf,'tag','ROC')
text(.5,.3,[ 'RSB reconstruction '  num2str(SCKSResult.EQM.ratiooz,'%2.2f  ')])
xlabel('1-Specificity'), ylabel('Sensitivity')

if 0
    % plot activité vs reconstruction
    subplot(2,1,1)
    plot(q(4).amp,(q(1).amp-q(4).amp).^2,'x')
    subplot(2,1,2)
    plot(q(4).amp,q(1).amp,'x')
end
%
% figure(8)
% index=[  5 11] ;
% for i0=1:length(index)
%     i1=index(i0);
%     plot(q(i1).SCKSNoStim(:,4),q(i1).SCKSNoStim(:,3),coul{i0},'displayname',['SCKS ' q(i1).name2 ', area = ' num2str(q(i1).infoB.ROCarea,'%1.2f')])
%     hold on
%     plot(q(i1).SCKSAll(:,4),q(i1).SCKSAll(:,3),[coul{i0} ':'],'displayname',['SCKS ' q(i1).name2 'All, area = ' num2str(q(i1).infoAll.ROCarea,'%1.2f')])
% end
% hold off
% a=legend('show');
% set(a,'location','southeast')
% title('ROC, Standard et Power'),set(gcf,'tag','ROC')
% text(.5,.3,[ 'RSB reconstruction '  num2str(SCKSResult.ratiooz,'%2.2f  ')])
% xlabel('1-Specificity'), ylabel('Sensitivity')

%  graph % compare SCKSdtc et electrogld
index=[  3 4 5 8] ;
for i0=1:length(index)
    i1=index(i0);
    figure(100+i1);
    %     plotSensi(q(i1).gldStimB,q(i1).gldElecSNoStim,q(i1).vct,q(1).vct,stimAll,S,{q(i1).name2,q(1).name2});
    plotSensi(q(i1).gldStimB,q(i1).gldElecSNoStimB,q(i1).vct,q(1).vct,stim,S,{q(i1).name2,q(1).name2});
end
1;
if 1   % seuil STD
    for i0=1:length(index)
        i1=index(i0);
        figure(400+i1);
        %     plotSensi(q(i1).gldStimB,q(i1).gldElecSNoStim,q(i1).vct,q(1).vct,stimAll,S,{q(i1).name2,q(1).name2});
        plotSensi(q(i1).gldStimB,q(i1).gldElecSNoStimC,q(i1).vct,q(1).vct,stim,S,{q(i1).name2,q(1).name2});
        %         plotSensi(q(i1).gldStimB,q(i1).gldElecSNoStimD,q(i1).vct,q(1).vct,stim,S,{q(i1).name2,q(1).name2});
    end
end

if 0  %  graph % compare SCKSgld et electrodtc
    for i0=1:length(index)
        i1=index(i0);
        figure(300+i1);
        plotSensi(q(1).gldStimB,q(1).gldSCKSSNoStim,q(1).vct,q(i1).vct,stimAll,S,{q(1).name2 ,q(i1).name2});
    end
end
if 0  %forward model
    
    for i0=1:length(index)
        i1=index(i0);
        figure(200+i1);
        plotSensi(q(i1).info2,q(i1).infoB,q(i1).vct,vctF,stim,S,{q(i1).name2,q(1).name2});
    end
end


if 0 % graph xy
    x=q(1).amp;
    stim=find(q(1).gldElecSNoStimB.vct.mask==1);
    noStim=find(q(1).gldElecSNoStimB.vct.mask==0);
    xstim=x(stim);
    xNS=x(noStim);
    y=q(5).amp;
    ystim=y(stim);
    yNS=y(noStim);
    subplot(2,1,1)
    
    plot(SCKS.Y(1,:),SCKS.qU.r{1}(1,:),'.') , xlabel('HbR, observed'), ylabel(('HbR, expected'))
    subplot(2,1,2)
    plot(x,y,'.'), xlabel('LFP, observed'), ylabel(('U(t), electro, expected'))
    printpdf(gcf,'b:\observedExpected.png','png',[],1)
end





function [score scoGld scoDtc ]=sensi(Gld, Dtc, mask,S)
%Gld: gold standard
%Dtc: data series to search
%dst (e.g. dt): window size for detection, in seconds
% mask : mask ==0 selected.   mask==1 stim not selected
%S: ROC structure
%find sensitivity matrix
if isfield(Gld,'thld')
    Gld=Gld.ind(Gld.amp>Gld.thld);
else
    Gld=Gld.ind;
end
if isfield(Dtc,'thld')
    Dtc=Dtc.ind(Dtc.amp>Dtc.thld);
else
    Dtc=Dtc.ind;
end

if isfield(mask,'ind')
    mask=mask.ind;
end
indexdtc=ceil(Dtc/S.freq/S.dtwindow);
indexgld=ceil(Gld/S.freq/S.dtwindow);
nBox=max([indexdtc,indexgld]);
if isempty(nBox); error('doonnées vide, vérifier intervalle de temps'); end
aa=zeros(nBox,1);
bb=aa;
aa(indexgld)=1;
bb(indexdtc)=1;

score.gldp2=length(Gld);
score.dtcp2=length(Dtc);
score.mask = mask;
if isempty(mask    )  % if mask empty select all data
    mask=zeros(size(aa));
else
    indexmask = ceil(mask/S.freq/S.dtwindow);
    mask = aa*0;
    indexmask=min(unique(max(indexmask,1)),nBox);
    mask(indexmask) = 1;
    
end



score.vct.mask = mask';
score.tp2=sum(aa&bb&~mask); %exclusive mask
score.fn2=sum(aa&~bb&~mask);
score.fp2=sum(~aa&bb&~mask);
score.tn2=sum(~aa&~bb&~mask);
score.fp2=sum(~aa&bb&~mask);
score.tpfp2=score.tp2+score.fp2;
score.gldpct=sum(aa&~mask); % percentage of time a true occurs un gld
score.dtcpct=sum(bb&~mask)/sum(~mask); % percentage of time a true occurs un dtc
% end
score.sensi2=score.tp2/(score.tp2+score.fn2);
score.speci2=score.tn2/(score.tn2+score.fp2);
score.gldpct2=score.tp2/length(Gld);
score.dtcPct2=score.tp2/length(Dtc);
score.ppv=score.tp2/(score.tp2+score.fp2);
score.vct.gld2=aa';
score.vct.dtc2=bb';
score.vct.temps2=(1:nBox)*S.dtwindow';







function [  out]=roc(gld,dtc,mask1,S,thldGld)
%vctDtc (e.g. vctB): time series to search
%gld (e.g. stim): time series considered as gold standard
%S: structure of ROC
%thldGld (e.g. S.s_th): threshold, used only after ROC . Find the
%threshold having a define sensibility
%rocIB :  [0 0 sensibility 1-specificity threshold];
%out : structure containing detection amplitude and position

if isempty(mask1)
    mask1.ind=[];
end
rocIB=[];
ii=1;
found_th = 0; %Boolean to search for point along ROC curve near thresholdSensi

if any( S.ROC_method==[2,3,4])
    %find if vctDtc exceeds threshold in each interval
    
    % fix acceptable limit  (remove extreme values)
    aa=sort(dtc.vct);
    max2 = aa(round(length(aa)*.99))*1.1; %
    min2 =  aa(round(length(aa)*.01))*1.1;
    bb=meanReg(aa,6,0);
    out.mean=mean(bb);
    out.std=std(bb);
    out.min2=min2;
    out.max2=max2;
    ratiobmax = max2;
    ratiobdelta = (max2-min2)/100;
    ratiob = min2-ratiobdelta;
end


while ratiob<=ratiobmax+ratiobdelta && ii<1010
    ratiob=ratiob+ratiobdelta;
    %select proper stim
    
    dtc.thld=ratiob;
    if nargin == 5  % select just the goldstandard that cross a certain threshold
        gld.thld= thldGld ;
        
    end
    stimB=sensi(gld,dtc, mask1.ind,S);
    % form ROC
    rocIB(ii,1)=stimB.sensi2;
    rocIB(ii,2)=1-stimB.speci2;
    rocIB(ii,3)=ratiob;
    rocIB(ii,4)=stimB.ppv;
    rocIB(ii,5)=stimB.tpfp2;
    rocIB(ii,6)=stimB.tp2;
    if any( S.ROC_method==[2,3,4])
        if rocIB(ii,3) == 0
            ratiob = ratiobmax+2*ratiobdelta; %to exit the while
        end
    end
    ii=ii+1;
end

%Compute area under ROC curve
if rocIB(1,1)~=1; rocIB=[ 1 1 -5 1 1 1; rocIB]; end
if rocIB(end,1)~=0; rocIB=[ rocIB;0 0 5 0 0 0]; end
rocIB( isnan(rocIB))=0;
meanROC=(rocIB(1:end-1,1)+rocIB(2:end,1))/2;
out.ROCarea = -sum(meanROC.*diff(rocIB(:,2)));
out.roc=rocIB;

function plotSensi(infoGld,infoDtc,vctGld,vctDtc,stim,S,names)
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

%set(gcf,'tag','doScale')

hold off
return
if 0  % try to find what correlate
    %
    bpm =extPhysio(physio,'bpm'); %Local field potential
    bpm=filtredecim(bpm,S.decim);
    bpm=filtfilt(b,a,bpm);
    bpm = bpm(SCKS.sv);
    pre =extPhysio(physio,'Press'); %Local field potential
    pre=filtredecim(pre,S.decim);
    pre=filtfilt(b,a,pre);
    pre = pre(SCKS.sv);
    
    aa=[ SCKS.pU.r' vctL*5 bpm/20  pre/20];
    
    bb=repmat((1:size(aa,2))*.15,size(aa,1),1);
    [b1,a1] = butter(5,0.005,'high');
    aa=filtfilt(b1,a1,aa);
    a= plot(time0 ,[aa+bb]);
    set(a(1),'color','B'), set(a(2),'color','G'), set(a(3),'color','K'), set(a(4),'color','M');
    legend({'HbR','HbT','Flow','LFP' 'bpm','Pression'})
end

function out=findEQM(Ele,vctB,out,time0)
%find difference between original and reconstruction
dt=[30 5 ];
[b,a] = butter(4,1./dt./(5/2),'bandpass');


vctB=vctB-mean(vctB);

% field=fieldnames(Ele);
for i1=1:length(Ele)
    ff=Ele(i1).name;
    electro=Ele(i1).vct-mean(Ele(i1).vct);
    out.EQM.(ff).scalesElectro = (electro\vctB); % scales
    scaleElectro=electro*(electro\vctB); % scales
    Ele1.(ff)=scaleElectro;
    out.EQM.( ff).elecIni=mean(electro.^2);
    out.EQM.( ff).Elec=mean(scaleElectro.^2);
    out.EQM.( ff).B=mean(vctB .^2);
    out.EQM.( ff).BminusElec=mean((vctB-scaleElectro).^2);
    out.EQM.( ff).BminusElecRatio=mean((vctB-scaleElectro).^2)/mean(vctB.^2);
    %     out.([ ff 'ratio'])=out.EQM.( ff).BminusElecRatio;
    electroFilt=filtfilt(b,a,electro);
    
    scaleElectroFilt=electroFilt*(electroFilt\vctB); % scales
    Ele2.(ff)=scaleElectroFilt;
    out.EQM.( ff).BminusElecFiltRatio=mean((vctB-scaleElectroFilt).^2)/mean(vctB.^2);
    %     out.([ ff 'Filtratio'])=out.EQM.( ff).BminusElecFiltRatio;
    
end
return
1;
plot(time0, [vctB Ele2.L Ele2.LFPL Ele2.LFPH Ele2.LFPLPw] )
legend({'b','LFPL','LFPH','LOW' 'LFPLPw'})
figure
plot(time0, [vctB Ele1.L Ele1.LFPL Ele1.LFPH] )
legend({'b','LFPL','LFPH','LOW'})

function varargout=formResample(S,vct,mask)
%form new dtc with different time intervall
S.wind2s = S.sk_npts2:S.winpts:S.n; %array of starting points
S.n2=floor(S.n/S.winpts);
mat=reshape(vct(1:(S.n2*S.winpts)),S.winpts,S.n2);
S.nwin2 = size(mat,1); %number of windows
t3.ind=[];
switch S.ROC_method
    case 2
        [t2amp t2ind]=max(mat,[],1);
    case 3
        %compute signed area
        t2amp= mean(mat,1);
        t2ind = t2amp*0+S.winpts/2; %just put index in center of window
    case 4
        [t2amp t3ind]=max(mat,[],1);
        t2ind = t2amp*0+S.winpts/2; %
end
t2ind=t2ind+S.wind2s(1:length(t2ind))-1;
t3ind=t3ind+S.wind2s(1:length(t3ind))-1;

if nargin==3
    ind=intersect(1:length(t2ind),find(mask));
    t2amp=t2amp(ind);
    t2ind=t2ind(ind);
    t3ind=t3ind(ind);
end

if nargout==1;
    varargout{1}.indMax=t3ind;
    varargout{1}.amp=t2amp;
    varargout{1}.ind=t2ind;
    varargout{1}.vct=vct;
else
    varargout={t2amp t2ind t3ind vct };
end


function stim=findStim(vctS,S)
[stim.amp stim.ind]=findpeaks(vctS,'minpeakdistance',S.dtFreq,'minpeakheight',.01);
if S.UseFullStims
    tmp1 = []; tmp2 = [];
    for i1=1:length(stim.ind)
        switch mod(i1,3)
            case 0
                ct1=15; %15 seconds stimulations
            case 1
                ct1=1; %1 second
            case 2
                ct1=5; %5 seconds
        end
        tmp3 = []; tmp4 = [];
        for i2=1:ct1
            tmp3 = [tmp3 stim.ind(i1)+(i2-1)*S.freq];
            tmp4 = [tmp4 stim.amp(i1)];
        end
        tmp1 = [tmp1 tmp3];
        tmp2 = [tmp2 tmp4];
    end
    stim.ind = tmp1;
    stim.amp = tmp2;
    %     vctS=vctS*0;
    %     vctS(stim.ind)=1;
    %     S.ROC_method=3;
    %     stim=formResample(S,vctS);
    
end

function vctDtc=shift(vctDtc,S)
% by phylippe not used

if S.ShiftedPower
    error
    vctDtc = vctDtc-min(vctDtc); %shift to positive values
    vctDtc = vctDtc.^2; %power
end
addn = round(S.shift*S.freq);
if S.shift < 0
    error
    vctDtc = vctDtc(abs(addn):end);
elseif addn>0
    error
    vctDtc = [zeros(addn,1); vctDtc];
end
function q=reduce1(q)
%reduce size for saving
for i1=2:length(q)
    
    field= fieldnames(q(1));
    for i2=1:length(field)
        data=q(i1).(field{i2});
        if strncmp('roc',field{i2},3)
        elseif isstruct(data)
            field2= fieldnames(data);
            for i3=1:length(field2)
                data2=data.(field2{i3});
                if strncmp('roc',field2{i3},3)
                elseif isstruct(data2)
                    q(i1).(field{i2}).(field2{i3})=[];
                elseif length(data2)>300
                    q(i1).(field{i2}).(field2{i3})=[];
                end
            end
            
        elseif length(data)>300
            q(i1).(field{i2}) =[];
            
        end
    end
end



function out=extractDtc(gld,dtc,mask,S,ROC,mode,thld)
% for a given threshold extract the TN FN TP FP
% mask : mask ==0 selected.   mask==1 stim not selected

switch(mode)
    case 'sensi'
        index=find(ROC(:,1)<=thld,1);
        dtc.thld=ROC(index,3);
    case 'speci'  %1-speci
        index=find(ROC(:,2)<=thld,1);
        dtc.thld=ROC(index,3);
    case 'thld'
        index=find(ROC(:,3)>=thld,1);
        dtc.thld=ROC(index,3);
    case 'std'
        out=sensi(gld,dtc, mask,S);
        aa=dtc.amp(out.vct.mask==0);
        dtc.thld= mean(aa)+thld*std(aa);
        %         disp(['pct detection outside stim =' num2str(mean(aa>dtc.thld)*100)]);
        index=1;
    otherwise
        error
end



out=sensi(gld,dtc, mask,S);
out.mode=mode;
out.s_th=thld;

out.s_th_ampli= dtc.thld;


function [mask1 mask2 mask3 mask1NS mask2NS mask3NS]=maskSensi(time0,info1,stim,stimAll)
% form stim 1 5 et 15 s
for i1=1:3
    %     mask.ind=stim.ind;
    %     mask.amp=stim.amp*0;
    mask.ind=[];
    mask.amp=[];
    maskNS=stimAll;
    index=find(info1.stim1(:,7)==i1);
    index2=min(index+1,length(info1.stim1));
    for i2=1:length(stim.ind)
        if any(time0(stim.ind(i2))>(info1.stim1(index,3)-3) & time0(stim.ind(i2))<(info1.stim1(index2,3)-3))
            %             mask.ind=stim.ind(i2);
            %             mask.amp(i2)=1;
        else
            mask.ind(end+1)=stim.ind(i2);
            mask.amp(end+1)=1;
            maskNS.ind(end+1)=stim.ind(i2);
            maskNS.amp(end+1)=1;
        end
    end
    [void ind]=unique(maskNS.ind);
    maskNS.ind= maskNS.ind(ind);
    maskNS.amp= maskNS.amp(ind);
    eval(['mask' num2str(i1) '=mask;'])
    eval(['mask' num2str(i1) 'NS=maskNS;'])
end
