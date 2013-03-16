%% script overlay blend
groupID = 'NC';
figFolder = 'D:\Edgar\Documents\Dropbox\Docs\Carotid\Figures\aligned';

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
job(1).figCmap                                  = jet(256);
job(1).figRange                                 = [-1 1]; % Fixed for correlation values
job(1).figIntensity                             = 0.8;
job(1).figAlpha                                 = 1;
job(1).transM                                   = rotz(pi);
job(1).figSize                                  = [1.5 1.5];
job(1).figRes                                   = 300;
job(1).drawCircle(1).drawCircle_On(1).circleLW  = 0.8;
job(1).drawCircle(1).drawCircle_On(1).circleLS  = '-';
job.parent_results_dir{1}                       = fullfile(figFolder,'overlay');

% Main loop
for iFiles = 1:numel(images2overlay)
    ioi_overlay_blend(IOImat, job, images2overlay{iFiles});
end

% EOF
