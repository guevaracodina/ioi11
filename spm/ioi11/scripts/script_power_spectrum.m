%% Script that plots spectrum of raw, filtered, and regressed ROI time course
% Loading data
dirName = 'D:\Edgar\Data\IOS_Res\12_07_03,CO06';
load(fullfile(dirName,'IOI.mat'))
% Raw time course
load(IOI.ROI.ROIfname)
% Filtered and downsampled data
load(IOI.fcIOS.filtNdown.fname)
% Filtered, downsampled and regressed data
load(IOI.fcIOS.SPM.fnameROIregress)

%% Get power spectra
% Choose seed/ROI
r1 = 3;
% Choose session
s1 = 4;
% Choose color
c1 = 5;
% Color names
colorNames = fieldnames(IOI.color);
% Original sampling frequency
Fs = 1/IOI.dev.TR;
% Downsampling frequency
downFs = IOI.fcIOS.filtNdown.fs;
% Raw time course
[Xraw, freq] = ioi_positiveFFT(ROI{r1}{s1, c1}, Fs);
% Filtered/downsampled time course
[XfiltNdown, freqfiltNdown] = ioi_positiveFFT(filtNdownROI{r1}{s1, c1}, downFs);
% Raw time course
[Xregressed, freqregressed] = ioi_positiveFFT(ROIregress{r1}{s1, c1}, downFs);

% Choose line style
switch c1
    case 1
        plotSymbol = 'r:';
    case 2
        plotSymbol = 'g:';
    case 3
        plotSymbol = 'y:';
    case 4
        plotSymbol = 'k:';
    case 5
        plotSymbol = 'r-';
    case 6
        plotSymbol = 'b-';
    case 7
        plotSymbol = 'k-';
    otherwise
        plotSymbol = 'c-';
end

% Display power spectra
close all
figure;
set(gcf,'color','w')

subplot(311)
semilogx(freq, abs(Xraw),plotSymbol,'LineWidth',3)
xlabel('f [Hz]','FontSize',12)
ylabel('Power [a.u.]','FontSize',12)
set(gca,'FontSize',12)
title(sprintf('%s raw time-course for seed %d S%d(%s)',IOI.subj_name,r1,s1,colorNames{1+c1}),'interpreter', 'none')
xlim([min(freq) max(freq)])

subplot(312)
semilogx(freqfiltNdown, abs(XfiltNdown),plotSymbol,'LineWidth',3)
line([IOI.fcIOS.filtNdown.BPFfreq(1) IOI.fcIOS.filtNdown.BPFfreq(1)], ...
    [min(abs(XfiltNdown)) max(abs(XfiltNdown))],...
    'Color','k','LineWidth',2,'LineStyle','--');
line([IOI.fcIOS.filtNdown.BPFfreq(2) IOI.fcIOS.filtNdown.BPFfreq(2)], ...
    [min(abs(XfiltNdown)) max(abs(XfiltNdown))],...
    'Color','k','LineWidth',2,'LineStyle','--');
xlabel('f [Hz]','FontSize',12)
ylabel('Power [a.u.]','FontSize',12)
set(gca,'FontSize',12)
title(sprintf('Filtered/downsampled time-course for seed %d S%d(%s)',r1,s1,colorNames{1+c1}))
xlim([min(freq) max(freq)])

subplot(313)
semilogx(freqregressed, abs(Xregressed),plotSymbol,'LineWidth',3)
line([IOI.fcIOS.filtNdown.BPFfreq(1) IOI.fcIOS.filtNdown.BPFfreq(1)], ...
    [min(abs(XfiltNdown)) max(abs(XfiltNdown))],...
    'Color','k','LineWidth',2,'LineStyle','--');
line([IOI.fcIOS.filtNdown.BPFfreq(2) IOI.fcIOS.filtNdown.BPFfreq(2)], ...
    [min(abs(XfiltNdown)) max(abs(XfiltNdown))],...
    'Color','k','LineWidth',2,'LineStyle','--');
xlabel('f [Hz]','FontSize',12)
ylabel('Power [a.u.]','FontSize',12)
set(gca,'FontSize',12)
title(sprintf('Filtered/downsampled and regressed time-course for seed %d S%d(%s)',r1,s1,colorNames{1+c1}))
xlim([min(freq) max(freq)])
