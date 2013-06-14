%% Load data
clear all
load('E:\Edgar\Data\IOS_Results\12_08_24,EG04\Onsets\IOI.mat')
% third session (post 4-AP)
s1 = 3;
load(IOI.res.el2{s1})
load(IOI.res.elinfo{s1})

%% Filter data
% The signal was filtered between 10-5000 Hz, amplified 1000× with a
% microelectrode AC amplifier (model 1800, A-M systems, Sequim, WA), then it was
% digitized at 10 kHz. LFP data was further filtered between 0.2 and 130 Hz
% using an order 4 Butterworth digital filter in post-processing
fs = 1/ConvertedData.Data.MeasuredData(1,12).Property(1,3).Value;
cutoff = [0.2 130];
FilterOrder = 4;
% el2 in (mV), because ConvertedData.Data.MeasuredData(1,12).Data has a 1000× amplification 
% See ConvertedData.Data.MeasuredData(1,12).Property(1,6).Value for units
t = (0:numel(el2)-1)'.*ConvertedData.Data.MeasuredData(1,12).Property(1,3).Value;
el2filt = temporalBPF('butter', fs, cutoff, FilterOrder, el2);

%% Electrophysiology of 4-AP seizure.
% Limits of the ictal discharge
ictalDischargeLimits = [318 392];
% Time window for onset/ictus/offset
tWin = 4;
% Zoomed-in windows
onsetLimits = [ictalDischargeLimits(1) ictalDischargeLimits(1)+tWin];
ictusLimits = [mean(ictalDischargeLimits) mean(ictalDischargeLimits)+tWin];
offsetLimits = [ictalDischargeLimits(2)-tWin ictalDischargeLimits(2)];

% Font sizes
axisFont = 14;
axisLabelFont = 16;
dottedLineWidth = 2;

% Plot results
h = figure; set(gcf,'color','w')
% A, An example of LFP recording demonstrates several ictal discharges following
% a single 4-AP injection. The arrow highlights the 4-AP injection time.
subplot(311)
plot(t,el2,'k-'); hold on
axis tight; box off; set(gca,'FontSize',axisFont)
yAxisLim = get(gca,'yLim');
% ylabel('LFP [mV]','FontSize',axisLabelFont)
% Plot ictal discharge limits
plot([ictalDischargeLimits(1) ictalDischargeLimits(1)],...
    [yAxisLim(1) yAxisLim(2)],'r--','LineWidth',dottedLineWidth);
plot([ictalDischargeLimits(2) ictalDischargeLimits(2)],...
    [yAxisLim(1) yAxisLim(2)],'r--','LineWidth',dottedLineWidth);


% B, An expanded view of a single ictal discharge.
subplot(312)
plot(t,el2,'k-'); hold on
xlim([ictalDischargeLimits(1)-1 ictalDischargeLimits(2)+1]); 
ylim(yAxisLim); box off; set(gca,'FontSize',axisFont)
ylabel('LFP [mV]','FontSize',axisLabelFont)
% Plot onset limits
plot([onsetLimits(1) onsetLimits(1)],[yAxisLim(1) yAxisLim(2)],'r--','LineWidth',dottedLineWidth);
plot([onsetLimits(2) onsetLimits(2)],[yAxisLim(1) yAxisLim(2)],'r--','LineWidth',dottedLineWidth);
% Plot ictus limits
plot([ictusLimits(1) ictusLimits(1)],[yAxisLim(1) yAxisLim(2)],'r--','LineWidth',dottedLineWidth);
plot([ictusLimits(2) ictusLimits(2)],[yAxisLim(1) yAxisLim(2)],'r--','LineWidth',dottedLineWidth);
% Plot offset limits
plot([offsetLimits(1) offsetLimits(1)],[yAxisLim(1) yAxisLim(2)],'r--','LineWidth',dottedLineWidth);
plot([offsetLimits(2) offsetLimits(2)],[yAxisLim(1) yAxisLim(2)],'r--','LineWidth',dottedLineWidth);

% C, Further expended view shows the evolution of the ictal discharge starting
% with the onset (focal spike and recruiting rhythm), development of
% intermittent spike-and-wave activity, and gradual ictal offset. 
subplot(337)
plot(t,el2,'k-'); hold on
xlim(onsetLimits); ylim(yAxisLim); box off; set(gca,'FontSize',axisFont)
set(gca,'xTick',onsetLimits)
xlabel('Onset','FontSize',axisLabelFont)
% ylabel('LFP [mV]','FontSize',axisLabelFont)

subplot(338)
plot(t,el2,'k-'); hold on
xlim(ictusLimits); ylim(yAxisLim); box off; set(gca,'FontSize',axisFont)
set(gca,'xTick',ictusLimits)
set(gca,'yTick',[])
xlabel({'Ictus'; 'time [s]'}, 'FontSize', axisLabelFont)

subplot(339)
plot(t,el2,'k-'); hold on
xlim(offsetLimits); ylim(yAxisLim); box off; set(gca,'FontSize',axisFont)
set(gca,'xTick',offsetLimits)
set(gca,'yTick',[])
xlabel('Offset','FontSize',axisLabelFont)

%% Set windows real size
job.figSize = [6 2.5];
job.figRes = 300;
% Specify window units
set(h, 'units', 'inches')
% Change figure and paper size
set(h, 'Position', [0.1 0.1 job.figSize(1) job.figSize(2)])
set(h, 'PaperPosition', [0.1 0.1 job.figSize(1) job.figSize(2)])

% Save as PNG at the user-defined resolution
print(h, '-dpng', ...
    fullfile('D:\Edgar\Documents\Dropbox\Docs\Epilepsy\figs','LFP_4AP'),...
    sprintf('-r%d',job.figRes));
% .fig file huge (~195 Mb)
% saveas(h, fullfile('D:\Edgar\Documents\Dropbox\Docs\Epilepsy\figs','LFP_4AP.fig'), 'fig');

% Return the property to its default
set(h, 'units', 'pixels')

%%
addpath(genpath('D:\Edgar\ssoct\Matlab'))
export_fig(fullfile('D:\Edgar\Documents\Dropbox\Docs\Epilepsy\figs','LFP_4AP'),'-png',h)
% EOF
