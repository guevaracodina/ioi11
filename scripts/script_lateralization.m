% script_lateralization
%% Load data
clear; clc
dataFolder = 'D:\Edgar\OIS_Results\lateralization';

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

%% Compute lateralization index
nLPS = numel(LPSIOImat);
nNaCl = numel(NaClIOImat);
for iLPS = 1:nLPS,
    load (LPSIOImat{iLPS})
    load(IOI.fcIOS.corr.corrMatrixFname)
    HbOMat = seed2seedCorrMat{1}{5};
    HbRMat = seed2seedCorrMat{1}{6};
    for iR = 1:2:numel(IOI.ROIname),
        LPS.HbO.LL(fix(iR/2)+1) = mean(HbOMat(iR+1,4:2:10));
        LPS.HbO.LR(fix(iR/2)+1) = mean(HbOMat(iR+1,3:2:9));
        LPS.HbO.RL(fix(iR/2)+1) = mean(HbOMat(iR,4:2:10));
        LPS.HbO.RR(fix(iR/2)+1) = mean(HbOMat(iR,3:2:9));
        
        LPS.HbR.LL(fix(iR/2)+1) = mean(HbRMat(iR+1,4:2:10));
        LPS.HbR.LR(fix(iR/2)+1) = mean(HbRMat(iR+1,3:2:9));
        LPS.HbR.RL(fix(iR/2)+1) = mean(HbRMat(iR,4:2:10));
        LPS.HbR.RR(fix(iR/2)+1) = mean(HbRMat(iR,3:2:9));
        
        fcliLPS.HbO(fix(iR/2)+1, iLPS) = ioi_lateralization_idx (LPS.HbO.LL(fix(iR/2)+1),...
            LPS.HbO.LR(fix(iR/2)+1), LPS.HbO.RL(fix(iR/2)+1), LPS.HbO.RR(fix(iR/2)+1));
        fcliLPS.HbR(fix(iR/2)+1, iLPS) = ioi_lateralization_idx (LPS.HbR.LL(fix(iR/2)+1),...
            LPS.HbR.LR(fix(iR/2)+1), LPS.HbR.RL(fix(iR/2)+1), LPS.HbR.RR(fix(iR/2)+1));
    end
end

for iNaCl = 1:nNaCl,
    load (NaClIOImat{iNaCl})
    load(IOI.fcIOS.corr.corrMatrixFname)
    HbOMat = seed2seedCorrMat{1}{5};
    HbRMat = seed2seedCorrMat{1}{6};
    for iR = 1:2:numel(IOI.ROIname),
        NaCl.HbO.LL(fix(iR/2)+1) = mean(HbOMat(iR+1,4:2:10));
        NaCl.HbO.LR(fix(iR/2)+1) = mean(HbOMat(iR+1,3:2:9));
        NaCl.HbO.RL(fix(iR/2)+1) = mean(HbOMat(iR,4:2:10));
        NaCl.HbO.RR(fix(iR/2)+1) = mean(HbOMat(iR,3:2:9));
        
        NaCl.HbR.LL(fix(iR/2)+1) = mean(HbRMat(iR+1,4:2:10));
        NaCl.HbR.LR(fix(iR/2)+1) = mean(HbRMat(iR+1,3:2:9));
        NaCl.HbR.RL(fix(iR/2)+1) = mean(HbRMat(iR,4:2:10));
        NaCl.HbR.RR(fix(iR/2)+1) = mean(HbRMat(iR,3:2:9));
        
        fcliNaCl.HbO(fix(iR/2)+1, iNaCl) = ioi_lateralization_idx (NaCl.HbO.LL(fix(iR/2)+1),...
            NaCl.HbO.LR(fix(iR/2)+1), NaCl.HbO.RL(fix(iR/2)+1), NaCl.HbO.RR(fix(iR/2)+1));
        fcliNaCl.HbR(fix(iR/2)+1, iNaCl) = ioi_lateralization_idx (NaCl.HbR.LL(fix(iR/2)+1),...
            NaCl.HbR.LR(fix(iR/2)+1), NaCl.HbR.RL(fix(iR/2)+1), NaCl.HbR.RR(fix(iR/2)+1));
    end
end

%% Perform stat test and save data
for iR = 1:2:numel(IOI.ROIname),
    [p.HbO(fix(iR/2)+1), h.HbO(fix(iR/2)+1)] = ranksum(fcliNaCl.HbO(fix(iR/2)+1, :), ...
        fcliLPS.HbO(fix(iR/2)+1, :));
    [p.HbR(fix(iR/2)+1), h.HbR(fix(iR/2)+1)] = ranksum(fcliNaCl.HbR(fix(iR/2)+1, :), ...
        fcliLPS.HbR(fix(iR/2)+1, :));
end

q.HbO = ioi_fdr(p.HbO);
q.HbR = ioi_fdr(p.HbR);

save(fullfile(dataFolder,'lateral_idx.mat'))

%% Create distribution plots
job.save_figures = true;
job(1).figSize                                  = [6 3];    % inches
job(1).figRes                                   = 300;      % in dpi

addpath('D:\Edgar\distributionPlot')
colorContrast = {[] [] [] [] 'r' 'b' 'g' 'm'};

%% Load cell with HBO data points
for iR = 1:2:numel(IOI.ROIname),
    dataPoints{iR} = fcliNaCl.HbO(fix(iR/2)+1, :);
    dataPoints{iR+1} = fcliLPS.HbO(fix(iR/2)+1, :);
end

%% Create HbO figure
hFig = figure; set(gcf,'color','w');
distributionPlot( dataPoints, ...
'color',...
repmat({colorContrast{5} 0.75*[1 1 1]}, [1, 10]),...
    'addSpread',true,'showMM',6,'variableWidth',true);
set(gca, 'XTick',1.5:2:20)       
set(gca, 'XTickLabel',{'V' 'C' 'S' 'R' 'M'},...
    'FontSize',14)
legend({'NaCl' 'LPS'},'FontSize',14,'Location','NorthEast')
ylabel('L.I.(\DeltaHbO_2) [a.u.]','FontSize',14)
ylim([-1 1])
% Specify window units
set(hFig, 'units', 'inches')
% Change figure and paper size
set(hFig, 'Position', [0.1 0.1 job.figSize(1) job.figSize(2)])
set(hFig, 'PaperPosition', [0.1 0.1 job.figSize(1) job.figSize(2)])
if job.save_figures
    % Save as PNG at the user-defined resolution
    print(hFig, '-dpng', ...
        fullfile(dataFolder, sprintf('LI_HbO')),...
        sprintf('-r%d',job.figRes));
    % Return the property to its default
    set(hFig, 'units', 'pixels')
    close(hFig)
end

%% Load cell with HBR data points
for iR = 1:2:numel(IOI.ROIname),
    dataPoints{iR} = fcliNaCl.HbR(fix(iR/2)+1, :);
    dataPoints{iR+1} = fcliLPS.HbR(fix(iR/2)+1, :);
end

%% Create HbR figure
hFig = figure; set(gcf,'color','w');
distributionPlot( dataPoints, ...
'color',...
repmat({colorContrast{6} 0.75*[1 1 1]}, [1, 10]),...
    'addSpread',true,'showMM',6,'variableWidth',true);
set(gca, 'XTick',1.5:2:20)       
set(gca, 'XTickLabel',{'V' 'C' 'S' 'R' 'M'},...
    'FontSize',14)
legend({'NaCl' 'LPS'},'FontSize',14,'Location','NorthEast')
ylabel('L.I.(\DeltaHbR) [a.u.]','FontSize',14)
ylim([-1 1])
% Specify window units
set(hFig, 'units', 'inches')
% Change figure and paper size
set(hFig, 'Position', [0.1 0.1 job.figSize(1) job.figSize(2)])
set(hFig, 'PaperPosition', [0.1 0.1 job.figSize(1) job.figSize(2)])
if job.save_figures
    % Save as PNG at the user-defined resolution
    print(hFig, '-dpng', ...
        fullfile(dataFolder, sprintf('LI_HbR')),...
        sprintf('-r%d',job.figRes));
    % Return the property to its default
    set(hFig, 'units', 'pixels')
    close(hFig)
end

% EOF