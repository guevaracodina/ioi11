%% Load data
load('E:\Edgar\Data\IOS_Results\12_08_24,EG04\ROI_Fig1\ROI.mat');
% Region 1 (At the injection site)
r1 = 1;
% Session 3
s1 = 3;
% ROIs
a = ROI{r1}{s1,5};
b = ROI{r1}{s1,6};
% Filtering
af = ButterHPF(5,0.004,2,a);
bf = ButterHPF(5,0.004,2,b);
af = ButterLPF(5,0.08,2,af);
bf = ButterLPF(5,0.08,2,bf);

%% Plot resting state somatosensory seeds time traces.
% Font sizes
axisFont = 14;
axisLabelFont = 16;
dottedLineWidth = 2;
markSize = 12;
% xLimits = [200 500];
% Time window for onset/ictus/offset
tWin = 4;
ictalDischargeLimits = [200 500]; %[4-tWin 392+300+tWin];
% time vector idx2plot
idx2plot = (ictalDischargeLimits(1):ictalDischargeLimits(2));

h = figure; set(gcf,'color','w')
hold on
% HbO
plot(idx2plot, af(idx2plot), 'r-', 'LineWidth', dottedLineWidth, 'Color', [161 0 62]/255)
% HbR
plot(idx2plot, bf(idx2plot), 'b-', 'LineWidth', dottedLineWidth, 'Color', [1 139 255]/255)
axis tight
xlabel('time [s]','FontSize',axisLabelFont);
ylabel('\DeltaHb [\muM]','FontSize',axisLabelFont);
legend({'HbO_2'; 'HbR'}, 'Location', 'NorthWest')
set(gca,'FontSize',axisFont);

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
    fullfile('D:\Edgar\Documents\Dropbox\Docs\Epilepsy\figs','Fig1D'),...
    sprintf('-r%d',job.figRes));
% Return the property to its default
set(h, 'units', 'pixels')

%% Print graphics
% addpath(genpath('D:\Edgar\ssoct\Matlab'))
% export_fig(fullfile('D:\Edgar\Documents\Dropbox\Docs\Epilepsy\figs','Fig1D'),'-png',gcf)

% EOF
