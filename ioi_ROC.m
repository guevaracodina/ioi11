function ROC=ioi_ROC(SCKS,O)
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
S.set_limits = 0; %Boolean: 1: use S.xlim for x-axis limits, 0: display whole series
%S.xlim=[1 200]; %Display limits
S.xlim=[1 2000]; %Display limits
S.decim=3; %Low pass filter frequency, in Hz
S.freq=5; %Acquisition frequency in Hz
S.dtwindow=2; %detection window length in seconds  % need to be a multiple of the acquisition time 0.2s
S.s_th = 0.95; %0.95; %Sensitivity threshold
S.sk_npts2 = 1; %becomes an offset??? - point 1 is time zero
S.sk_npts = 0; %5; %number of points set to zero at beginning of backward pass SCKS deconvolved inputs
S.dtFreq=round(S.dtwindow*S.freq);
S.winpts = round(S.dtwindow*S.freq);
S.min_peak = 0.5;
%FIND STIM index
stim0 = SCKS.pU.v{1};
vctS = full(stim0);
%time0 = (1:length(vctS))*SCKS.dt;
S.stimInd=find(vctS>0.05)';
vctB=SCKS.qU.v{1}';vctB(1:S.sk_npts)=0;
vctB = ButterHPF(5,0.02,3,vctB);
%FIND STIM index
S.n = size(SCKS.Y.y,1);
stim.indGld =S.stimInd; 
stim.ampGld = ones(1,length(S.stimInd));

% find ROC for SCKS
[rocSB infoSB vctB] = ioi_roc_detect(vctB,stim,S,S.s_th); 
%disp(['SCKS: ' num2str(infoSB.ROCarea,'%1.3f')]);
%Store
ROC.rocSB = rocSB;
ROC.infoSB = infoSB;

h1 = figure;
plot(rocSB(:,4),rocSB(:,3));
leg={['SCKS, area = ' num2str(infoSB.ROCarea,'%1.3f')]};
legend(leg,'location','southeast')
title(O.tit) %,set(gcf,'tag','ROC')
xlabel('1-Specificity'), ylabel('Sensitivity')
ioi_save_figures(O.save_figures,O.generate_figures,h1,O.fname,O.dir);

%plot ROC
ioi_plotSensi(infoSB,vctB,stim,S,{'Stim','SCKS'},O);
