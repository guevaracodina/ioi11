%% Filtered & downsampled data Load ROI 13 time traces
ioiMat='E:\Edgar\Data\IOS_Results\12_09_28,EG09\Onsets\GLMfcIOS\corrMap\ROI\Series\FiltNDown\IOI.mat';
load(ioiMat)
% Define session and color to use later
s1 = 3;
% time limits (see script_seizure_duration.m)
tBegin = 329;
tEnd = 465;
c1 = 5;
% ROI
r1 = 13;
load(IOI.fcIOS.filtNdown.fname)
roi13 = filtNdownROI{r1};

%% Load ROI 8 time trace
ioiMat = 'E:\Edgar\Data\IOS_Results\12_09_28,EG09\Onsets\IOI.mat';
load(ioiMat)
load(IOI.fcIOS.filtNdown.fname)
filtNdownROI{r1} = roi13; clear roi13
% ROI time vector
tROI = (0:numel(filtNdownROI{r1}{s1, c1})-1)/IOI.fcIOS.filtNdown.downFreq;

%% Load LFP data
load(IOI.res.el2{s1})
% LFP time vector
tEl = (0:numel(el2)-1)'/IOI.res.sfel;
tBefore = 10;
tAfter = 10;
tLim = [tBegin-tBefore tEnd+tAfter];
idx1 = find(tROI>=(tBegin-tBefore),1,'first');
idx2 = find(tROI>=tBegin,1,'first');

%% Display filtered & downsampled data
h = figure;
set(gcf,'color','w')
axisFontSize = 18;
lineWidth = 3;

% ROI
subplot(211)

useBaseValue = false;

c1 = 5;
r1 = 13;
if useBaseValue
    baseValue = mean(filtNdownROI{r1}{s1, c1}(idx1:idx2));
    plot(tROI, (filtNdownROI{r1}{s1, c1}-baseValue)./baseValue,'r-','LineWidth',lineWidth)
else
    plot(tROI, filtNdownROI{r1}{s1, c1}/100,'r-','LineWidth',lineWidth)
end
hold on
r1 = 8;
if useBaseValue
    baseValue = mean(filtNdownROI{r1}{s1, c1}(idx1:idx2));
    plot(tROI, (filtNdownROI{r1}{s1, c1}-baseValue)./baseValue,'r--','LineWidth',lineWidth)
else
    plot(tROI, filtNdownROI{r1}{s1, c1}/100,'r--','LineWidth',lineWidth)
end

c1 = 6;
r1 = 13;
if useBaseValue
    baseValue = mean(filtNdownROI{r1}{s1, c1}(idx1:idx2));
    plot(tROI, (filtNdownROI{r1}{s1, c1}-baseValue)./baseValue,'b-','LineWidth',lineWidth)
else
    plot(tROI, filtNdownROI{r1}{s1, c1}/100,'b-','LineWidth',lineWidth)
end
r1 = 8;
if useBaseValue
    baseValue = mean(filtNdownROI{r1}{s1, c1}(idx1:idx2));
    plot(tROI, (filtNdownROI{r1}{s1, c1}-baseValue)./baseValue,'b--','LineWidth',lineWidth)
else
    plot(tROI, filtNdownROI{r1}{s1, c1}/100,'b--','LineWidth',lineWidth)
end


c1 = 7;
r1 = 13;
if useBaseValue
    baseValue = mean(filtNdownROI{r1}{s1, c1}(idx1:idx2));
    plot(tROI, (filtNdownROI{r1}{s1, c1}-baseValue)./baseValue,'k-','LineWidth',lineWidth)
else
    plot(tROI, filtNdownROI{r1}{s1, c1}/1e4,'k-','LineWidth',lineWidth)
end
r1 = 8;
if useBaseValue
    baseValue = mean(filtNdownROI{r1}{s1, c1}(idx1:idx2));
    plot(tROI, (filtNdownROI{r1}{s1, c1}-baseValue)./baseValue,'k--','LineWidth',lineWidth)
else
    plot(tROI, filtNdownROI{r1}{s1, c1}/1e4,'k--','LineWidth',lineWidth)
end

c1 = 8;
r1 = 13;
if useBaseValue
    baseValue = mean(filtNdownROI{r1}{s1, c1}(idx1:idx2));
    plot(tROI, (filtNdownROI{r1}{s1, c1}-baseValue)./baseValue,'c-','LineWidth',lineWidth)
else
    plot(tROI, filtNdownROI{r1}{s1, c1},'c-','LineWidth',lineWidth)
end
r1 = 8;
if useBaseValue
    baseValue = mean(filtNdownROI{r1}{s1, c1}(idx1:idx2));
    plot(tROI, (filtNdownROI{r1}{s1, c1}-baseValue)./baseValue,'c--','LineWidth',lineWidth)
else
    plot(tROI, filtNdownROI{r1}{s1, c1},'c--','LineWidth',lineWidth)
end

legend({'HbO_2 focus' 'HbO_2 surround' 'HbR focus' 'HbR surround' ...
    'CBF focus' 'CBF surround' 'CMRO_2 focus' 'CMRO_2 surround'},...
    'FontSize', axisFontSize)
xlim(tLim)
set(gca, 'xTick', [])
set(gca, 'FontSize', axisFontSize)
ylabel('\Delta [a.u.]', 'FontSize', axisFontSize)

% LFP
subplot(212)
plot(tEl, el2, 'k-')
axis tight
xlim(tLim)
xlabel('t [s]', 'FontSize', axisFontSize)
ylabel('LFP [mV]', 'FontSize', axisFontSize)
set(gca, 'FontSize', axisFontSize)

%% Print graphics
addpath(genpath('D:\Edgar\ssoct\Matlab'))
export_fig(fullfile('D:\Edgar\Documents\Dropbox\Docs\Epilepsy\figs','seizure_traces'),'-png',gcf)

% EOF
