%% script_overlay_blend
% String to identify the group
groupID = 'NC';
% Folder where to save the images
% figFolder = 'D:\Edgar\Documents\Dropbox\Docs\Carotid\Figures\aligned\';
figFolder = 'F:\Edgar\Data\IOS_Carotid_Res2\alignment\figures\aligned_CBF\';
% Range of values to map to the full range of colormap: [minVal maxVal]
fcMapRange = [-0.2733 0.9052];
% Range of values to map to display non-transparent pixels: [minVal maxVal]
alphaRange = [-1 1];

%% Code
dirListNIfTI = dir(fullfile(figFolder, [groupID '*.nii']));
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
job(1).figCmap                                  = jet(256);     % colormap
job(1).figIntensity                             = 1;            % [0 - 1]
job(1).transM                                   = rotz(pi);     % affine transform
job(1).figSize                                  = [1.5 1.5];    % inches
job(1).figRes                                   = 300;          % in dpi
job(1).drawCircle(1).drawCircle_On(1).circleLW  = 0.8;          % line width
job(1).drawCircle(1).drawCircle_On(1).circleLS  = '-';          % line style
job(1).drawCircle(1).drawCircle_On(1).circleEC  = 'w';          % line color
job.parent_results_dir{1}                       = fullfile(figFolder,'overlay');
job.generate_figures                            = true;         % display figure
job.save_figures                                = false;        % save figure
% ------------------------------------------------------------------------------

% Main loop
for iFiles = 1:numel(images2overlay)
    [displayed_pixels{iFiles}, total_pixels(iFiles)] = ioi_overlay_blend(IOImat, job, images2overlay{iFiles}, fcMapRange, alphaRange);
end

% EOF
