%% Load data
clear; clc
dataFolder = 'D:\Edgar\OIS_Results\extent';

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
alphaVal = 0.000001;

for c1 = 5:6,                       % Contrast Loop
    %% LPS
    for iLPS = 1:nLPS,              % Subjects loop
        load(LPSIOImat{iLPS})
        load(IOI.fcIOS.corr.fname)
        % Brain Mask
        vol = spm_vol(IOI.fcIOS.mask.fname);
        brainMask = logical(fix(ioi_MYimresize(spm_read_vols(vol), [512, 512])));
        myIdx = 1;
        for iR = 3:numel(seed_based_fcIOS_map), % ROI by ROI 
            pMap = seed_based_fcIOS_map{iR}{c1}.pValue;
            
            % FDR-correction
            pMask = ~isnan(pMap) & brainMask;
            pMapFDRtmp = ioi_fdr(pMap(pMask));
            pMapFDR = nan(size(pMap));
            pMapFDR(pMask) = pMapFDRtmp;
            
            % Apply threshold
            % Pmap_N_S: 'Binary map indicating significance at P<0.05 (fdr corrected)'
            Pmap_N_S = (pMapFDR <= alphaVal)  & brainMask;
            nSignifPixels = nnz(Pmap_N_S);
            ratSignifPixels = nSignifPixels / nnz(brainMask);
            switch(c1)
                case 5
                    LPSextent.HbO(myIdx, iLPS) = ratSignifPixels;
                case 6
                    LPSextent.HbR(myIdx, iLPS) = ratSignifPixels;
            end
            myIdx = myIdx + 1;
        end
    end
    
    %% NaCl
    for iNaCl = 1:nNaCl,
        load (NaClIOImat{iNaCl})
        load(IOI.fcIOS.corr.fname)
        % Brain Mask
        vol = spm_vol(IOI.fcIOS.mask.fname);
        brainMask = logical(fix(ioi_MYimresize(spm_read_vols(vol), [512, 512])));
        myIdx = 1;
        for iR = 3:numel(seed_based_fcIOS_map),
            pMap = seed_based_fcIOS_map{iR}{c1}.pValue;
            
            % FDR-correction
            pMask = ~isnan(pMap) & brainMask;
            pMapFDRtmp = ioi_fdr(pMap(pMask));
            pMapFDR = nan(size(pMap));
            pMapFDR(pMask) = pMapFDRtmp;
            
            % Apply threshold
            % Pmap_N_S: 'Binary map indicating significance at P<0.05 (fdr corrected)'
            Pmap_N_S = (pMapFDR <= alphaVal)  & brainMask;
            nSignifPixels = nnz(Pmap_N_S);
            ratSignifPixels = nSignifPixels / nnz(brainMask);
            switch(c1)
                case 5
                    NaClextent.HbO(myIdx, iNaCl) = ratSignifPixels;
                case 6
                    NaClextent.HbR(myIdx, iNaCl) = ratSignifPixels;
            end
            myIdx = myIdx + 1;
        end
    end
end

%% ROI by ROI comparison
for iR = 1:8,
    [p.HbO(iR), h.HbO(iR)] = ranksum(LPSextent.HbO(iR, :), NaClextent.HbO(iR, :));
    [p.HbR(iR), h.HbR(iR)] = ranksum(LPSextent.HbR(iR, :), NaClextent.HbR(iR, :));
end

q.HbO = ioi_fdr(p.HbO);
q.HbR = ioi_fdr(p.HbR);

% Save data
save(fullfile(dataFolder,'spatial_extent_thresh001.mat'), 'p', 'h', 'q', 'LPSextent', 'NaClextent')

%% job options
% ------------------------------------------------------------------------------
% Define matlab batch job with the required fields
% ------------------------------------------------------------------------------
job(1).figSize                                  = [3.5 3.5];    % inches
job(1).figRes                                   = 300;          % in dpi
job.generate_figures                            = true;         % display figure
job.save_figures                                = true;        % save figure
% ------------------------------------------------------------------------------

%% ROI by ROI  Distribution plots (HbO)
addpath('D:\Edgar\distributionPlot')
colorContrast = {[] [] [] [] 'r' 'b' 'g' 'm'};
% Load cell with data points
for iR = 1:8,
    dataPoints{2*iR-1} = 100*NaClextent.HbO(iR, :);
    dataPoints{2*iR} = 100*LPSextent.HbO(iR, :);
end
% Create figure    
hFig = figure; set(gcf,'color','w');
distributionPlot( dataPoints, ...
'color',...
repmat({colorContrast{5} 0.75*[1 1 1]}, [1, 8]),...
    'addSpread',true,'showMM',6,'variableWidth',true);
set(gca, 'XTick',1.5:2:16)       
set(gca, 'XTickLabel',IOI.ROIname(3:10),...
    'FontSize',14)
legend({'NaCl' 'LPS'},'FontSize',14,'Location','NorthEast')
ylabel('Significant pixels [%]','FontSize',14)
ylim([60 120])
% Specify window units
set(hFig, 'units', 'inches')
% Change figure and paper size
set(hFig, 'Position', [0.1 0.1 job.figSize(1) job.figSize(2)])
set(hFig, 'PaperPosition', [0.1 0.1 job.figSize(1) job.figSize(2)])
if job.save_figures
    % Save as PNG at the user-defined resolution
    print(hFig, '-dpng', ...
        fullfile(dataFolder, sprintf('spatial_extent_HbO')),...
        sprintf('-r%d',job.figRes));
    % Return the property to its default
    set(hFig, 'units', 'pixels')
    close(hFig)
end

%% ROI by ROI Distribution plots (HbR)
% Load cell with data points
for iR = 1:8,
    dataPoints{2*iR-1} = 100*NaClextent.HbR(iR, :);
    dataPoints{2*iR} = 100*LPSextent.HbR(iR, :);
end
% Create figure    
hFig = figure; set(gcf,'color','w');
distributionPlot( dataPoints, ...
'color',...
repmat({colorContrast{6} 0.75*[1 1 1]}, [1, 8]),...
    'addSpread',true,'showMM',6,'variableWidth',true);
set(gca, 'XTick',1.5:2:16)       
set(gca, 'XTickLabel',IOI.ROIname(3:10),...
    'FontSize',14)
legend({'NaCl' 'LPS'},'FontSize',14,'Location','NorthEast')
ylabel('Significant pixels [%]','FontSize',14)
ylim([60 120])
% Specify window units
set(hFig, 'units', 'inches')
% Change figure and paper size
set(hFig, 'Position', [0.1 0.1 job.figSize(1) job.figSize(2)])
set(hFig, 'PaperPosition', [0.1 0.1 job.figSize(1) job.figSize(2)])
if job.save_figures
    % Save as PNG at the user-defined resolution
    print(hFig, '-dpng', ...
        fullfile(dataFolder, sprintf('spatial_extent_HbR')),...
        sprintf('-r%d',job.figRes));
    % Return the property to its default
    set(hFig, 'units', 'pixels')
    close(hFig)
end

% -------------------------------------------------------------------------
%% Global spatial extent
clear; close all; clc
% addpath('D:\Edgar\distributionPlot')
addpath('C:\Edgar\Dropbox\Matlab\distributionPlot')
load('C:\Edgar\Dropbox\PostDoc\Newborn\OIS_Results\averaged_maps\spatial_extent_LPS.mat',...
    'nSignifPixelsLPS', 'nPixelsLPS')
load('C:\Edgar\Dropbox\PostDoc\Newborn\OIS_Results\averaged_maps\spatial_extent_NaCl.mat',...
    'nSignifPixelsNaCl', 'nPixelsNaCl')
nSignifPixelsLPS = nSignifPixelsLPS(3:end, 5:6);
nPixelsLPS = nPixelsLPS(3:end, 5:6);
nSignifPixelsNaCl = nSignifPixelsNaCl(3:end, 5:6);
nPixelsNaCl = nPixelsNaCl(3:end, 5:6);
extLPS.HbO = nSignifPixelsLPS(:,1) ./ nPixelsLPS(:,1);
extLPS.HbR = nSignifPixelsLPS(:,2) ./ nPixelsLPS(:,2);
extNaCl.HbO = nSignifPixelsNaCl(:,1) ./ nPixelsNaCl(:,1);
extNaCl.HbR = nSignifPixelsNaCl(:,2) ./ nPixelsNaCl(:,2);

%% Stat test
[p.HbO, h.HbO] = ranksum(extLPS.HbO, extNaCl.HbO);
[p.HbR, h.HbR] = ranksum(extLPS.HbR, extNaCl.HbR);

%% Prepare data points
dataPoints{1} = 100*extNaCl.HbO;
dataPoints{2} = 100*extLPS.HbO;
dataPoints{3} = 100*extNaCl.HbR;
dataPoints{4} = 100*extLPS.HbR;

%% Create figure
hFig = figure; set(gcf,'color','w');
distributionPlot( dataPoints, ...
'color',...
{[1 0 0] 0.75*[1 1 1] [0 0 1] 0.75*[1 1 1]},...
    'addSpread',true,'showMM',6,'variableWidth',true);
set(gca, 'XTick',1:4)       
set(gca, 'XTickLabel',{'NaCl' 'LPS' 'NaCl' 'LPS'},...
    'FontSize',14)
% legend({'NaCl' 'LPS'},'FontSize',14,'Location','NorthEast')
ylabel('Significant pixels [%]','FontSize',14)
ylim([-20 100])
% Specify window units
set(hFig, 'units', 'inches')
% Change figure and paper size
set(hFig, 'Position', [0.1 0.1 job.figSize(1) job.figSize(2)])
set(hFig, 'PaperPosition', [0.1 0.1 job.figSize(1) job.figSize(2)])
if job.save_figures
    % Save as PNG at the user-defined resolution
    print(hFig, '-dpng', ...
        fullfile('C:\Edgar\Dropbox\PostDoc\Newborn\OIS_Results\averaged_maps',...
        sprintf('spatial_extension_tmp.png')),...
        sprintf('-r%d',job.figRes));
    % Return the property to its default
    set(hFig, 'units', 'pixels')
    close(hFig)
end

%% Create figure to show differences
clc
% 'C_R' 'C_L' 'S_R' 'S_L' 'R_R' 'R_L' 'M_R' 'M_L'
ROInames = {'C_R' 'C_L' 'S_R' 'S_L' 'R_R' 'R_L' 'M_R' 'M_L'};
hFig = figure; set(hFig,'color','w')
NaCl_coordinates = 0.5:2:16;
LPS_coordinates = 1.5:2:16.5;
MrkSz = 8;
LinWdt = 2;

% ---------------------------- HbO2 ---------------------------------------
subplot(211)
hold on
% plot([0.75:8, 1.25:8.25],[dataPoints{1}; dataPoints{2}]','--')
plot(NaCl_coordinates, dataPoints{1},'ro', 'MarkerSize',MrkSz, 'LineWidth',LinWdt,...
    'MarkerFaceColor','r')
plot(LPS_coordinates, dataPoints{2},'rx', 'MarkerSize',MrkSz, 'LineWidth',LinWdt)
j = 1;
for i=fix((LPS_coordinates + NaCl_coordinates) / 2),
    plot([NaCl_coordinates(1) LPS_coordinates(1)]'+i-1,[dataPoints{1}(j); dataPoints{2}(j)]','r--')
    j = j + 1;
end
plot(NaCl_coordinates, dataPoints{1},'ro', 'MarkerSize',MrkSz, 'LineWidth',LinWdt,...
    'MarkerFaceColor','r')
plot(LPS_coordinates, dataPoints{2},'rx', 'MarkerSize',MrkSz, 'LineWidth',LinWdt)
set(gca, 'xTick', fix((LPS_coordinates + NaCl_coordinates) / 2), 'xTickLabel', ROInames,'FontSize',14)
xlabel('ROI','FontSize',14)
% ylabel('Significant pixels [%]','FontSize',14)
% legend({'NaCl_{HbO_2}' 'LPS_{HbO_2}'})
set(gca, 'FontSize', 14)
ylim([0 70])
set(gca, 'yTick',[0 35 70]);

% ---------------------------- HbR ---------------------------------------
subplot(212)
hold on
plot(NaCl_coordinates, dataPoints{3},'bo', 'MarkerSize',MrkSz, 'LineWidth',LinWdt,...
    'MarkerFaceColor','b')
plot(LPS_coordinates, dataPoints{4},'bx', 'MarkerSize',MrkSz, 'LineWidth',LinWdt)
j = 1;
for i=fix((LPS_coordinates + NaCl_coordinates) / 2),
    plot([NaCl_coordinates(1) LPS_coordinates(1)]'+i-1,[dataPoints{3}(j); dataPoints{4}(j)]','b--')
    j = j + 1;
end
plot(NaCl_coordinates, dataPoints{3},'bo', 'MarkerSize',MrkSz, 'LineWidth',LinWdt,...
    'MarkerFaceColor','b')
plot(LPS_coordinates, dataPoints{4},'bx', 'MarkerSize',MrkSz, 'LineWidth',LinWdt)
set(gca, 'xTick', fix((LPS_coordinates + NaCl_coordinates) / 2), 'xTickLabel', ROInames,'FontSize',14)
xlabel('ROI','FontSize',14)
ylabel('Significant pixels [%]','FontSize',14)
% legend({'NaCl_{HbR}' 'LPS_{HbR}'})
set(gca, 'FontSize', 14)
ylim([0 70])
set(gca, 'yTick',[0 35 70]);

% Specify window units
set(hFig, 'units', 'inches')
% Change figure and paper size
set(hFig, 'Position', [0.1 0.1 job.figSize(1) job.figSize(2)])
set(hFig, 'PaperPosition', [0.1 0.1 job.figSize(1) job.figSize(2)])

if job.save_figures
    % Save as PNG at the user-defined resolution
    print(hFig, '-dpng', ...
        fullfile('C:\Edgar\Dropbox\PostDoc\Newborn\OIS_Results\averaged_maps',...
        sprintf('extent_comparison.png')),...
        sprintf('-r%d',job.figRes));
    % Return the property to its default
    set(hFig, 'units', 'pixels')
    close(hFig)
end
% EOF
% Comment as test