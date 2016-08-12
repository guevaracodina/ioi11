
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
alphaVal = 0.05;

for c1 = 5:6,                       % Contrast Loop
    %% LPS
    for iLPS = 1:nLPS,
        load(LPSIOImat{iLPS})
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

%% ROI
for iR = 1:8,
    [p.HbO(iR), h.HbO(iR)] = ranksum(LPSextent.HbO(iR, :), NaClextent.HbO(iR, :));
    [p.HbR(iR), h.HbR(iR)] = ranksum(LPSextent.HbR(iR, :), NaClextent.HbR(iR, :));
end

q.HbO = ioi_fdr(p.HbO);
q.HbR = ioi_fdr(p.HbR);

% Save data
save(fullfile(dataFolder,'spatial_extent.mat'), 'p', 'h', 'q', 'LPSextent', 'NaClextent')

%% job options
% ------------------------------------------------------------------------------
% Define matlab batch job with the required fields
% ------------------------------------------------------------------------------
job(1).figSize                                  = [6 3];    % inches
job(1).figRes                                   = 300;          % in dpi
job.generate_figures                            = true;         % display figure
job.save_figures                                = false;        % save figure
% ------------------------------------------------------------------------------

%% Distribution plots (HbO)
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

%% Distribution plots (HbR)
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

% EOF
