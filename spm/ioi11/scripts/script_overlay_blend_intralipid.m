%% script_overlay_blend
% String to identify the group
groupID = 'WL';
% Folder where to save the images
figFolder = 'F:\Edgar\Data\IOS_Resolution\Results\averaged_maps\CBF\R04\';
% Range of values to map to the full range of colormap: [minVal maxVal]
fcMapRange = [-1 1];
% Range of values to map to display non-transparent pixels: [minVal maxVal]
alphaRange = [-1 1];

%% Code
dirListNIfTI = dir(fullfile(figFolder, [groupID '*.nii']));
% Load IOI matrix of the source image
if strcmp(groupID, 'WL')
    IOImat = 'F:\Edgar\Data\IOS_Resolution\Results\13_07_30,WL01\IOI.mat';
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
job(1).transM                                   = rotz(pi/2)*roty(pi/2);     % affine transform
job(1).figSize                                  = [1.5 1.5];    % inches
job(1).figRes                                   = 300;          % in dpi
job(1).drawCircle(1).drawCircle_On(1).circleLW  = 0.8;          % line width
job(1).drawCircle(1).drawCircle_On(1).circleLS  = '-';          % line style
job(1).drawCircle(1).drawCircle_On(1).circleEC  = 'w';          % line color
job.parent_results_dir{1}                       = fullfile(figFolder,'overlay');
job.generate_figures                            = true;         % display figure
job.save_figures                                = true;        % save figure
% ------------------------------------------------------------------------------
nColorLevels = 256;
r1 = 4;
% Main loop
for iFiles = 1:numel(images2overlay)
    [displayed_pixels{iFiles}, total_pixels(iFiles)] = ioi_overlay_blend_intralipid(...
        IOImat, job, images2overlay{iFiles}, fcMapRange, alphaRange, nColorLevels, r1);
end

%% Average images
% for i=1:10
%     vol(i) = spm_vol(fullfile('F:\Edgar\Data\IOS_Resolution\Results\averaged_maps\CBF\R04',['13_07_30,WL' sprintf('%02d',i) '_R04_S01_C7_fcIOS_map.nii']));
%     I(:,:,i) = spm_read_vols(vol(i));
% end
% avg_I = mean(I,3);
% ioi_save_nifti(avg_I, 'F:\Edgar\Data\IOS_Resolution\Results\averaged_maps\CBF\R04\AVG_R04_C7.nii', [1 1 1]);
% EOF
