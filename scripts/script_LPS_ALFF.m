%% Load data
clear; clc
dataFolder = 'D:\Edgar\OIS_Results\alff';

LPSIOImat{1} = 'D:\Edgar\OIS_Results\16_02_25,LP01a\ROI\LPF\FiltNDown\GLM\corrMap\IOI.mat';
LPSIOImat{2} = 'D:\Edgar\OIS_Results\16_02_25,LP01b\ROI\LPF\FiltNDown\GLM\corrMap\IOI.mat';
LPSIOImat{3} = 'D:\Edgar\OIS_Results\16_07_07,LP02a\ROI\LPF\FiltNDown\GLM\corrMap\IOI.mat';
LPSIOImat{4} = 'D:\Edgar\OIS_Results\16_07_07,LP02b\ROI\LPF\FiltNDown\GLM\corrMap\IOI.mat';
LPSIOImat{5} = 'D:\Edgar\OIS_Results\16_07_08,LP03\ROI\LPF\FiltNDown\GLM\corrMap\IOI.mat';
LPSIOImat{6} = 'D:\Edgar\OIS_Results\16_07_08,LP04\ROI\LPF\FiltNDown\GLM\corrMap\IOI.mat';
LPSIOImat{7} = 'D:\Edgar\OIS_Results\16_07_08,LP05a\ROI\LPF\FiltNDown\GLM\corrMap\IOI.mat';

NaClIOImat{1} = 'D:\Edgar\OIS_Results\16_02_25,NC01\ROI\LPF\FiltNDown\GLM\corrMap\IOI.mat';
NaClIOImat{2} = 'D:\Edgar\OIS_Results\16_02_25,NC02\ROI\LPF\FiltNDown\GLM\corrMap\IOI.mat';
NaClIOImat{3} = 'D:\Edgar\OIS_Results\16_02_25,NC03a\ROI\LPF\FiltNDown\GLM\corrMap\IOI.mat';
NaClIOImat{4} = 'D:\Edgar\OIS_Results\16_02_25,NC03b\ROI\LPF\FiltNDown\GLM\corrMap\IOI.mat';
NaClIOImat{5} = 'D:\Edgar\OIS_Results\16_02_26,NC05a\ROI\LPF\FiltNDown\GLM\corrMap\IOI.mat';
NaClIOImat{6} = 'D:\Edgar\OIS_Results\16_02_26,NC05b\ROI\LPF\FiltNDown\GLM\corrMap\IOI.mat';
NaClIOImat{7} = 'D:\Edgar\OIS_Results\16_02_26,NC06a\ROI\LPF\FiltNDown\GLM\corrMap\IOI.mat';
NaClIOImat{8} = 'D:\Edgar\OIS_Results\16_02_26,NC06b\ROI\LPF\FiltNDown\GLM\corrMap\IOI.mat';
NaClIOImat{9} = 'D:\Edgar\OIS_Results\16_07_07,NC07\ROI\LPF\FiltNDown\GLM\corrMap\IOI.mat';
NaClIOImat{10} = 'D:\Edgar\OIS_Results\16_07_08,NC08\ROI\LPF\FiltNDown\GLM\corrMap\IOI.mat';
NaClIOImat{11} = 'D:\Edgar\OIS_Results\16_07_08,NC09\ROI\LPF\FiltNDown\GLM\corrMap\IOI.mat';

nLPS = numel(LPSIOImat);
nNaCl = numel(NaClIOImat);
for iLPS = 1:nLPS,
    load (LPSIOImat{iLPS})
    load(IOI.fcIOS.mask.fnameSeries)
    yGlobalHbO = ioi_alff(brainMaskSeries{1}{5});
    yGlobalHbR = ioi_alff(brainMaskSeries{1}{6});
    load(IOI.ROI.ROIfname)
    for iR = 1:numel(ROI),
        HbOLPS(iLPS, iR) = ioi_alff(ROI{iR}{5}, yGlobalHbO);
        HbRLPS(iLPS, iR) = ioi_alff(ROI{iR}{6}, yGlobalHbR);
    end
end

for iNaCl = 1:nNaCl,
    load (NaClIOImat{iNaCl})
    load(IOI.fcIOS.mask.fnameSeries)
    yGlobalHbO = ioi_alff(brainMaskSeries{1}{5});
    yGlobalHbR = ioi_alff(brainMaskSeries{1}{6});
    load(IOI.ROI.ROIfname)
    for iR = 1:numel(ROI),
        HbONaCl(iNaCl, iR) = ioi_alff(ROI{iR}{5}, yGlobalHbO);
        HbRNaCl(iNaCl, iR) = ioi_alff(ROI{iR}{6}, yGlobalHbR);
    end
end

for iR = 1:numel(ROI),
    [p.HbO(iR), h.HbO(iR)] = ranksum(HbOLPS(:, iR), HbONaCl(:, iR));
    [p.HbR(iR), h.HbR(iR)] = ranksum(HbRLPS(:, iR), HbRNaCl(:, iR));
end

q.HbO = ioi_fdr(p.HbO);
q.HbR = ioi_fdr(p.HbR);

% Save data
save(fullfile(dataFolder,'alff_vals.mat'))

%% Load data
clear; clc
dataFolder = 'D:\Edgar\OIS_Results\alff';
load(fullfile(dataFolder,'alff_vals.mat'))

%% job options
% ------------------------------------------------------------------------------
% Define anonymous functions for affine transformations
% ------------------------------------------------------------------------------
rotx = @(theta) [1 0 0 0; 0 cos(theta) -sin(theta) 0; 0 sin(theta) cos(theta) 0; 0 0 0 1];
roty = @(theta) [cos(theta) 0 sin(theta) 0; 0 1 0 0; -sin(theta) 0 cos(theta) 0; 0 0 0 1];
rotz = @(theta) [cos(theta) -sin(theta) 0 0; sin(theta) cos(theta) 0 0; 0 0 1 0; 0 0 0 1];
translate = @(a,b) [1 0 a 0; 0 1 b 0; 0 0 1 0; 0 0 0 1];
% ------------------------------------------------------------------------------
nColorLevels = 256;
% ------------------------------------------------------------------------------
% Define matlab batch job with the required fields
% ------------------------------------------------------------------------------
job(1).figCmap                                  = ioi_get_colormap('bipolar');     % colormap
job(1).figIntensity                             = 1;            % [0 - 1]
job(1).transM                                   = rotz(pi/2)*roty(pi/2);     % affine transform
job(1).figSize                                  = [6 3];    % inches
job(1).figRes                                   = 300;          % in dpi
job(1).drawCircle(1).drawCircle_On(1).circleLW  = 0.8;          % line width
job(1).drawCircle(1).drawCircle_On(1).circleLS  = '-';          % line style
job(1).drawCircle(1).drawCircle_On(1).circleEC  = 'k';          % line color
job.generate_figures                            = true;         % display figure
job.save_figures                                = false;        % save figure
% ------------------------------------------------------------------------------

%% Distribution plots (HbO)
addpath('D:\Edgar\distributionPlot')
colorContrast = {[] [] [] [] 'r' 'b' 'g' 'm'};
% Load cell with data points
for iR = 1:numel(IOI.ROIname),
    dataPoints{2*iR-1} = HbONaCl(:, iR);
    dataPoints{2*iR} = HbOLPS(:, iR);
end
% Create figure    
hFig = figure; set(gcf,'color','w');
distributionPlot( dataPoints, ...
'color',...
repmat({colorContrast{5} 0.75*[1 1 1]}, [1, 10]),...
    'addSpread',true,'showMM',6,'variableWidth',true);
set(gca, 'XTick',1.5:2:20)       
set(gca, 'XTickLabel',IOI.ROIname,...
    'FontSize',14)
legend({'NaCl' 'LPS'},'FontSize',14,'Location','NorthEast')
ylabel('Standardized ALFF','FontSize',14)
ylim([-1.5 7])
% Specify window units
set(hFig, 'units', 'inches')
% Change figure and paper size
set(hFig, 'Position', [0.1 0.1 job.figSize(1) job.figSize(2)])
set(hFig, 'PaperPosition', [0.1 0.1 job.figSize(1) job.figSize(2)])
if job.save_figures
    % Save as PNG at the user-defined resolution
    print(hFig, '-dpng', ...
        fullfile(dataFolder, sprintf('ALFF_HbO')),...
        sprintf('-r%d',job.figRes));
    % Return the property to its default
    set(hFig, 'units', 'pixels')
    close(hFig)
end

%% Distribution plots (HbR)
addpath('D:\Edgar\distributionPlot')
colorContrast = {[] [] [] [] 'r' 'b' 'g' 'm'};
% Load cell with data points
for iR = 1:numel(IOI.ROIname),
    dataPoints{2*iR-1} = HbRNaCl(:, iR);
    dataPoints{2*iR} = HbRLPS(:, iR);
end
% Create figure    
hFig = figure; set(gcf,'color','w');
distributionPlot( dataPoints, ...
'color',...
repmat({colorContrast{6} 0.75*[1 1 1]}, [1, 10]),...
    'addSpread',true,'showMM',6,'variableWidth',true);
set(gca, 'XTick',1.5:2:20)       
set(gca, 'XTickLabel',IOI.ROIname,...
    'FontSize',14)
legend({'NaCl' 'LPS'},'FontSize',14,'Location','NorthEast')
ylabel('Standardized ALFF','FontSize',14)
ylim([-1.5 7])
% Specify window units
set(hFig, 'units', 'inches')
% Change figure and paper size
set(hFig, 'Position', [0.1 0.1 job.figSize(1) job.figSize(2)])
set(hFig, 'PaperPosition', [0.1 0.1 job.figSize(1) job.figSize(2)])
if job.save_figures
    % Save as PNG at the user-defined resolution
    print(hFig, '-dpng', ...
        fullfile(dataFolder, sprintf('ALFF_HbR')),...
        sprintf('-r%d',job.figRes));
    % Return the property to its default
    set(hFig, 'units', 'pixels')
    close(hFig)
end


% EOF


