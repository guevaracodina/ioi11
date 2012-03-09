function ROC=ioi_ROC(SCKS,Options)
% to execute : sensitivitySCKS(SCKS,physio,info1)
%S: ROC structureaf
%ROCmethod:
%2: Set threshold on amplitude
%3: Set threshold on area
S.ROC_method = 3; %this is set for the frequency of 0.2 Hz. So changing this mode has no effect
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
S.min_peak = 0.5;
%FIND STIM index
stim0 = SCKS.pU.v{1};
vctS = full(stim0);
time0 = (1:length(vctS))*SCKS.dt;
S.stimInd=find(vctS>0.05)';
vctB=SCKS.qU.v{1}';vctB(1:S.sk_npts)=0;
%FIND STIM index
S.n = size(SCKS.Y.y,1);
%Find stimulation times 
%[stim.ampGld stim.indGld]=findpeaks(vctS,'minpeakdistance',S.dtFreq,'minpeakheight',S.min_peak);
stim.indGld =S.stimInd; 
stim.ampGld = ones(1,length(S.stimInd));
if S.UseFullStims
    tmp1 = []; tmp2 = [];
    for i1=1:length(stim.indGld)
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
            tmp3 = [tmp3 stim.indGld(i1)+(i2-1)*S.freq];
            tmp4 = [tmp4 stim.ampGld(i1)];
        end
        tmp1 = [tmp1 tmp3];
        tmp2 = [tmp2 tmp4];
    end               
    stim.indGld = tmp1;
    stim.ampGld = tmp2;
end

% find ROC for SCKS
[rocSB infoSB vctB] = ioi_roc_detect(vctB,stim,S,S.s_th); 
disp(['SCKS: ' num2str(infoSB.ROCarea,'%1.3f')]);

%Find ROC for simple threshold on flow, HbR, HbT
off0 = 0+1;
vctF = SCKS.Y(3,:)';
vctF = [vctF(off0:end); vctF(end)*ones(off0-1,1)];
[rocF infoF vctF]=roc(vctF,stim,S,S.s_th); disp(['Flow: ' num2str(infoF.ROCarea,'%1.3f')]);
vctR = -SCKS.Y(1,:)';
vctR = [vctR(off0:end); vctR(end)*ones(off0-1,1)];
[rocR infoR vctR]=roc(vctR,stim,S,S.s_th); disp(['HbR: ' num2str(infoR.ROCarea,'%1.3f')]);
vctT = SCKS.Y(2,:)';
vctT = [vctT(off0:end); vctT(end)*ones(off0-1,1)];
[rocT infoT vctT]=roc(vctT,stim,S,S.s_th); disp(['HbT: ' num2str(infoT.ROCarea,'%1.3f')]);

%plot ROC
plot(rocSB(:,4),rocSB(:,3),'-');
leg={['SCKS, area = ' num2str(infoSB.ROCarea,'%1.3f')]};
legend(leg,'location','southeast')
title('ROC'),set(gcf,'tag','ROC')
xlabel('1-Specificity'), ylabel('Sensitivity')

%plot stim
ioi_plotSensi(infoSB,vctB,stim,S,{'Stim','SCKS'});

S.s_th = 0.7; %Sensitivity threshold
[rocSB2 infoSB2 vctB2]=roc(vctB,vctL,S,S.s_th,stim,infoSL);
plot(rocSB2(:,4),rocSB2(:,3));
leg={['SCKS, area = ' num2str(infoSB2.ROCarea,'%1.3f')],};
legend(leg,'location','southeast')
title('ROC'),set(gcf,'tag','ROC')
xlabel('1-Specificity'), ylabel('Sensitivity')
plotSensi(infoSL,infoSB2,vctL,vctB2,stim,S,{'LFP','SCKS'});
