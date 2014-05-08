%% script_overlay_blend for 3MT contest
clc
% String to identify the group
groupID = 'NC';
% Folder where to save the images
figFolder = 'D:\Edgar\Data\IOS_Carotid_Res\alignment\aligned_HbO\';
% Range of values to map to the full range of colormap: [minVal maxVal]
% fcMapRange = [-0.2733 0.9052];
fcMapRange = [0 1.6];
% Range of values to map to display non-transparent pixels: [minVal maxVal]
alphaRange = [0.5 1.6+ eps] ;
% Print anatomy
PRINT_ANATOMICAL = false;

%% Code
% dirListNIfTI = dir(fullfile(figFolder, [groupID '*.nii']));
dirListNIfTI = dir(fullfile(figFolder, [groupID '*.img']));
% Load IOI matrix of the source image
if strcmp(groupID, 'CC')
    IOImat = 'D:\Edgar\Data\IOS_Carotid_Res\12_10_19,CC10\GLMfcIOS\corrMap\IOI.mat';
elseif strcmp(groupID, 'NC')
    IOImat = 'D:\Edgar\Data\IOS_Carotid_Res\12_10_18,NC09\GLMfcIOS\corrMap\IOI.mat';
else
    fprintf('No IOI matrix found for %s.\n', groupID);
    return
end
images2overlay = struct2cell(dirListNIfTI);
images2overlay = images2overlay(1,:)';
images2overlay = cellfun(@(x) fullfile(figFolder, [x ',1']), images2overlay, 'UniformOutput', false);
[images2overlay, sts] = cfg_getfile([1 Inf],'image','Select images',images2overlay, figFolder, '.*');
% ------------------------------------------------------------------------------
% Define anonymous functions for affine transformations
% ------------------------------------------------------------------------------
rotx = @(theta) [1 0 0 0; 0 cos(theta) -sin(theta) 0; 0 sin(theta) cos(theta) 0; 0 0 0 1];
roty = @(theta) [cos(theta) 0 sin(theta) 0; 0 1 0 0; -sin(theta) 0 cos(theta) 0; 0 0 0 1];
rotz = @(theta) [cos(theta) -sin(theta) 0 0; sin(theta) cos(theta) 0 0; 0 0 1 0; 0 0 0 1];
translate = @(a,b) [1 0 a 0; 0 1 b 0; 0 0 1 0; 0 0 0 1];
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
% Define matlab batch job with the required fields
% ------------------------------------------------------------------------------
job(1).figCmap                                  = pat_get_colormap('linlhot');     % colormap
job(1).figIntensity                             = 1;            % [0 - 1]
job(1).transM                                   = rotz(pi);     % affine transform
job(1).figSize                                  = [3 3];    % inches
job(1).figRes                                   = 300;          % in dpi
job(1).drawCircle                               = [];
% job(1).drawCircle(1).drawCircle_On(1).circleLW  = 0.8;          % line width
% job(1).drawCircle(1).drawCircle_On(1).circleLS  = '-';          % line style
% job(1).drawCircle(1).drawCircle_On(1).circleEC  = 'w';          % line color
job.parent_results_dir{1}                       = fullfile(figFolder,'overlay');
job.generate_figures                            = true;         % display figure
job.save_figures                                = true;        % save figure
% ------------------------------------------------------------------------------

%% Anatomical image
if PRINT_ANATOMICAL
    v           = spm_vol(fullfile(figFolder,'12_10_18,NC09_anat_S01.nii'));
    anatomical  = spm_read_vols(v);
    v           = spm_vol(fullfile(figFolder,'12_10_18,NC09_anat_brainmask.nii'));
    brainmask   = spm_read_vols(v);
    % orient images
    anatomical  = fliplr(anatomical');
    brainmask  = fliplr(brainmask');
    hFig = figure;
    imagesc(anatomical.*brainmask); axis image; colormap (gray(255))
    set(hFig, 'color', 'k')
    % Allow printing of black background
    set(hFig, 'InvertHardcopy', 'off');
    % Specify window units
    set(hFig, 'units', 'inches')
    % Change figure and paper size
    set(hFig, 'Position', [0.1 0.1 job.figSize(1) job.figSize(2)])
    set(hFig, 'PaperPosition', [0.1 0.1 job.figSize(1) job.figSize(2)])
    if job.save_figures
        % Save as PNG at the user-defined resolution
        print(hFig, '-dpng', ...
            fullfile(figFolder,'anatomical.png'),...
            sprintf('-r%d',job.figRes));
        % Return the property to its default
        set(hFig, 'units', 'pixels')
        close(hFig)
    end
end
%% Main loop
for iFiles = 1:numel(images2overlay)
    [displayed_pixels{iFiles}, total_pixels(iFiles)] = ioi_overlay_blend(IOImat, job, images2overlay{iFiles}, fcMapRange, alphaRange);
end

% EOF
