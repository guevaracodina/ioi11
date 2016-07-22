%% script_overlay_maps_LPS
clear
c1 = 5;             % Color index
r1 = 1:10;          % ROI index

%% Load anatomical & brainmask
clc; close all;
% Load IOI matrix of the source image
IOImat = 'D:\Edgar\OIS_Results\16_02_25,NC01\ROI\LPF\FiltNDown\GLM\corrMap\IOI.mat';
load(IOImat)
vol = spm_vol('D:\Edgar\OIS_Results\averaged_maps\AVG_Atlas.img');
% Must flip average anatomical image, createed with Fiji
Underlay = flipud(rot90(spm_read_vols(vol), 1));
Underlay = ioi_MYimresize(Underlay, [512 512]);
%     figure; imagesc(Underlay); axis image; colormap(ioi_get_colormap('wob'));
vol = spm_vol(IOI.fcIOS.mask.fname);
brainMask = rot90(spm_read_vols(vol), 1);
brainMask = logical(fix(ioi_MYimresize(brainMask, [512 512])));
%     figure; imagesc(brainMask);  axis image; colormap(ioi_get_colormap('wob'));
% Range of values to map to the full range of colormap: [minVal maxVal]
fcMapRange = [-1 1];
% Range of values to map to display non-transparent pixels: [minVal maxVal]
alphaRange = [-1 1];
% Spatial extension threshold
spatialThreshold = 0.5;

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
job(1).figSize                                  = [3.5 3.5];    % inches
job(1).figRes                                   = 300;          % in dpi
job(1).drawCircle(1).drawCircle_On(1).circleLW  = 0.8;          % line width
job(1).drawCircle(1).drawCircle_On(1).circleLS  = '-';          % line style
job(1).drawCircle(1).drawCircle_On(1).circleEC  = 'k';          % line color
job.generate_figures                            = true;         % display figure
job.save_figures                                = false;        % save figure
% ------------------------------------------------------------------------------

%% ROI loop
clear pSpatial hSpatial
for iR = r1
    % Folder where to save the images (Change accordingly to contrast c1)
    switch(c1)
        case 5
            figFolder = fullfile('D:\Edgar\OIS_Results\averaged_maps\HbO\',num2str(r1(iR)));
        case 6
            figFolder = fullfile('D:\Edgar\OIS_Results\averaged_maps\HbR\',num2str(r1(iR)));
    end
    job.parent_results_dir{1}                       = fullfile(figFolder,'overlay');
    % ------------------------------------------------------------------------------
    % Load LPS images
    % String to identify the group
    groupID = 'LP';
    LPSListNIfTI = dir(fullfile(figFolder, ['*' groupID '*.nii']));
    LPSimages = struct2cell(LPSListNIfTI);
    LPSimages = LPSimages(1,:)';
    LPSimages = cellfun(@(x) fullfile(figFolder, [x ',1']), LPSimages, 'UniformOutput', false);
    [LPSimages, sts] = cfg_getfile([1 Inf],'image','Select images',LPSimages, figFolder, '.*');
    clear LPS %LPS_spatial_extension
    % Reading loop
    for iFiles = 1:numel(LPSimages)
        vol = spm_vol(LPSimages{iFiles});
        LPS(:,:,iFiles) = spm_read_vols(vol) .* brainMask;
        pixelMask = false(size(brainMask));
        pixelMask(squeeze(LPS(:,:,iFiles)) > spatialThreshold) = true;
        % Spatial extension % defined as displayed to brain pixels ratio.
        LPS_spatial_extension(iR) = nnz(pixelMask) / nnz(brainMask);
    end
    % ------------------------------------------------------------------------------
    % Load NaCl images
    % String to identify the group
    groupID2 = 'NC';
    NaClListNIfTI = dir(fullfile(figFolder, ['*' groupID2 '*.nii']));
    NaClimages = struct2cell(NaClListNIfTI);
    NaClimages = NaClimages(1,:)';
    NaClimages = cellfun(@(x) fullfile(figFolder, [x ',1']), NaClimages, 'UniformOutput', false);
    [NaClimages, sts] = cfg_getfile([1 Inf],'image','Select images',NaClimages, figFolder, '.*');
    clear NaCl %NaCl_spatial_extension
    % Reading loop
    for iFiles = 1:numel(NaClimages)
        vol = spm_vol(NaClimages{iFiles});
        NaCl(:,:,iFiles) = spm_read_vols(vol) .* brainMask;
        pixelMask = false(size(brainMask));
        pixelMask(squeeze(NaCl(:,:,iFiles)) > spatialThreshold) = true;
        % Spatial extension % defined as displayed to brain pixels ratio.
        NaCl_spatial_extension(iR) = nnz(pixelMask) / nnz(brainMask);
        
    end

    % Compare spatial extension
    [pSpatial(iR), hSpatial(iR)] = ranksum(LPS_spatial_extension, NaCl_spatial_extension);
%     avg_I = nanmean(I,3);

% save(fullfile(figFolder, sprintf('stats_R%d_C%d.mat', r1(iR), c1)), 'pSpatial', 'hSpatial', ...
%     'NaCl_spatial_extension', 'LPS_spatial_extension', 'LPSimages', 'NaClimages',...
%     'LPS', 'NaCl', 'brainMask', 'Underlay')

end % ROI loop

%% Stat test
[pSpatial, hSpatial] = ranksum(LPS_spatial_extension, NaCl_spatial_extension);
qSpatial = ioi_fdr(pSpatial);
% save(fullfile('D:\Edgar\OIS_Results\averaged_maps\', sprintf('stats_C%d.mat', c1)), ...
%     'LPS_spatial_extension', 'NaCl_spatial_extension', 'pSpatial', 'hSpatial', 'qSpatial')

%% Distribution plots
load('D:\Edgar\OIS_Results\averaged_maps\stats_C5.mat')
addpath('D:\Edgar\distributionPlot')
colorContrast = {[] [] [] [] 'r' 'b'};
dataPoints{1} = NaCl_spatial_extension;
dataPoints{2} = LPS_spatial_extension;
load('D:\Edgar\OIS_Results\averaged_maps\stats_C6.mat')
dataPoints{3} = NaCl_spatial_extension;
dataPoints{4} = LPS_spatial_extension;
hFig = figure; set(gcf,'color','w');
distributionPlot( dataPoints, ...
'color',...
    {colorContrast{5} 0.75*[1 1 1] colorContrast{6} 0.75*[1 1 1]},...
    'addSpread',true,'showMM',6,'variableWidth',true);
set(gca, 'XTick',[1 2 3 4])       
set(gca, 'XTickLabel',{'NaCl_{HbO_2}' 'LPS_{HbO_2}' 'NaCl_{HbR}' 'LPS_{HbR}'},...
    'FontSize',14)
% xlabel({'NaCl_{HbO_2}'; 'LPS_{HbO_2}'; 'NaCl_{HbR}'; 'LPS_{HbR}'},...
%     'FontSize',14,'Interpreter', 'tex')
ylabel('Spatial extension','FontSize',14)
ylim([-0.05 0.3])
% Specify window units
set(hFig, 'units', 'inches')
% Change figure and paper size
set(hFig, 'Position', [0.1 0.1 job.figSize(1) job.figSize(2)])
set(hFig, 'PaperPosition', [0.1 0.1 job.figSize(1) job.figSize(2)])
if job.save_figures
    % Save as PNG at the user-defined resolution
    print(hFig, '-dpng', ...
        fullfile('D:\Edgar\OIS_Results\averaged_maps\', sprintf('spatial_extension')),...
        sprintf('-r%d',job.figRes));
    % Return the property to its default
    set(hFig, 'units', 'pixels')
    close(hFig)
end

%% Create statistical maps from saved data
clear; clc;
c1 = 6;
r1 = 1:10;
%% ROI loop
for iR = r1
    % Folder where to save the images (Change accordingly to contrast c1)
    switch(c1)
        case 5
            figFolder = fullfile('D:\Edgar\OIS_Results\averaged_maps\HbO\',num2str(iR));
        case 6
            figFolder = fullfile('D:\Edgar\OIS_Results\averaged_maps\HbR\',num2str(iR));
    end
    job.parent_results_dir{1}                       = fullfile(figFolder,'overlay');
    load(fullfile(figFolder, sprintf('stats_R%d_C%d.mat', iR, c1)), ...
    'LPS', 'NaCl', 'brainMask', 'Underlay')
    nLPS = size(LPS,3);
    nNaCl = size(NaCl,3);
    nPixels = nnz(brainMask);
    
    brainMask3DLPS = repmat(brainMask,1,1,nLPS);
    LPSresh = LPS(brainMask3DLPS);
    LPSresh = reshape(LPSresh, [nLPS nPixels]);
    
    brainMask3D = repmat(brainMask,1,1,nNaCl);
    NaClresh = NaCl(brainMask3D);
    NaClresh = reshape(NaClresh, [nNaCl nPixels]);
    
    p = zeros([nPixels 1]);
    h = zeros([nPixels 1]);
    zVal = zeros([nPixels 1]);
    rankVal = zeros([nPixels 1]);
    
    Pmap_N_S = nan([size(LPS,1) size(LPS,2)]);
    Tmap_N_S = nan([size(LPS,1) size(LPS,2)]);
    Bmap_N_S = nan([size(LPS,1) size(LPS,2)]);
    Hmap_N_S = nan([size(LPS,1) size(LPS,2)]);
    tic
    for iPixels = 1:nPixels
        if all(isnan(LPSresh(:,iPixels))) || all(isnan(NaClresh(:,iPixels)))
            p(iPixels) = nan;
            h(iPixels) = nan;
%             zVal(iPixels) = nan;
            rankVal(iPixels) = nan;
        else
%         [p(iPixels), h(iPixels), ST] = ranksum(LPSresh(:,iPixels), NaClresh(:,iPixels));
            [h(iPixels),p(iPixels),CI,ST] = ttest2(LPSresh(:,iPixels), NaClresh(:,iPixels));
            rankVal(iPixels) = ST.tstat;
        end
    end
    
    Pmap_N_S(brainMask) = p;
    Tmap_N_S(brainMask) = rankVal;
%     Bmap_N_S(brainMask) = rankVal;
    Bdiff = nanmean(NaCl,3) - nanmean(LPS,3);
    Bmap_N_S(brainMask) = Bdiff(brainMask);
    Hmap_N_S(brainMask) = h;
    
    disp(['Statistical comparison in: ' datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS')]);
    save(fullfile(figFolder, sprintf('statMap_R%d_C%d.mat', iR, c1)), 'Pmap_N_S',...
        'Tmap_N_S', 'Bmap_N_S', 'Hmap_N_S', ...
        'LPS', 'NaCl', 'brainMask', 'Underlay')
%     dualcodeImage
end

%% Load maps from saved data
clear; clc;
c1 = 5;
r1 = 1;
%% ROI loop
for iR = r1
    % Folder where to save the images (Change accordingly to contrast c1)
    switch(c1)
        case 5
            figFolder = fullfile('D:\Edgar\OIS_Results\averaged_maps\HbO\',num2str(iR));
        case 6
            figFolder = fullfile('D:\Edgar\OIS_Results\averaged_maps\HbR\',num2str(iR));
    end
    job.parent_results_dir{1} = fullfile(figFolder,'overlay');
    load(fullfile(figFolder, sprintf('statMap_R%d_C%d.mat', iR, c1)))
    nLPS = size(LPS,3);
    nNaCl = size(NaCl,3);
    nPixels = nnz(brainMask);
    Pmapresh = Pmap_N_S(brainMask);
    Pmapresh = ioi_fdr(Pmapresh);
    Pthresh = 0.001;
    PmapBW = Pmapresh < Pthresh;
%     PmapBW_img = zeros(size(brainMask));
    Pmap_N_S(brainMask) = PmapBW;
    figure; imagesc(Pmap_N_S); axis image; colormap(ioi_get_colormap('wob'))
%     dualcodeImage
end
% EOF
