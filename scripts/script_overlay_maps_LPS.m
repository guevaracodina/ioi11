%% script_overlay_maps_LPS
r1 = 1:10;             % ROI index
for iR = 1:numel(r1)
    c1 = 6;             % Color index
    % Folder where to save the images (Change accordingly to contrast c1)
    figFolder = fullfile('D:\Edgar\OIS_Results\averaged_maps\HbR\',num2str(r1(iR)));
    % String to identify the group
    groupID = 'LP';
    % Range of values to map to the full range of colormap: [minVal maxVal]
    fcMapRange = [-1 1];
    % Range of values to map to display non-transparent pixels: [minVal maxVal]
    alphaRange = [-1 1];
    
    %% Code
    dirListNIfTI = dir(fullfile(figFolder, ['*' groupID '*.nii']));
    % Load IOI matrix of the source image
    IOImat = 'D:\Edgar\OIS_Results\16_02_25,NC01\ROI\LPF\FiltNDown\GLM\corrMap\IOI.mat';
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
    nColorLevels = 256;
    % ------------------------------------------------------------------------------
    % Define matlab batch job with the required fields
    % ------------------------------------------------------------------------------
    job(1).figCmap                                  = ioi_get_colormap('bipolar');     % colormap
    job(1).figIntensity                             = 1;            % [0 - 1]
    job(1).transM                                   = rotz(pi/2)*roty(pi/2);     % affine transform
    job(1).figSize                                  = [1.5 1.5];    % inches
    job(1).figRes                                   = 300;          % in dpi
    job(1).drawCircle(1).drawCircle_On(1).circleLW  = 0.8;          % line width
    job(1).drawCircle(1).drawCircle_On(1).circleLS  = '-';          % line style
    job(1).drawCircle(1).drawCircle_On(1).circleEC  = 'k';          % line color
    job.parent_results_dir{1}                       = fullfile(figFolder,'overlay');
    job.generate_figures                            = true;         % display figure
    job.save_figures                                = false;        % save figure
    % ------------------------------------------------------------------------------
    
    % Main loop
    for iFiles = 1:numel(images2overlay)
        [displayed_pixels{iFiles}, total_pixels(iFiles)] = ioi_overlay_blend_intralipid(...
            IOImat, job, images2overlay{iFiles}, fcMapRange, alphaRange, nColorLevels, r1(iR), c1);
    end
    
    %% Average images
    % reload atlas
    IOImat = 'D:\Edgar\OIS_Results\16_02_25,NC01\ROI\LPF\FiltNDown\GLM\corrMap\IOI.mat';
    for i=1:size(images2overlay,1)
        vol = spm_vol(images2overlay{i});
        I(:,:,i) = spm_read_vols(vol);
    end
    avg_I = nanmean(I,3);
    ioi_save_nifti(avg_I, fullfile(figFolder, sprintf('%s_AVG_R%d_C%d.nii',groupID, r1(iR), c1)), [1 1 1]);
    [displayed_pixels{iFiles}, total_pixels(iFiles)] = ioi_overlay_blend_intralipid(...
        IOImat, job, fullfile(figFolder,sprintf('%s_AVG_R%d_C%d.nii',groupID, r1(iR), c1)), fcMapRange, alphaRange, nColorLevels, r1(iR), c1);
    
    %% Reprocess for controls NC
    groupID = 'NC';
    % Folder where to save the images
    % figFolder = 'D:\Edgar\OIS_Results\averaged_maps\HbO\1';
    % Range of values to map to the full range of colormap: [minVal maxVal]
    fcMapRange = [-1 1];
    % Range of values to map to display non-transparent pixels: [minVal maxVal]
    alphaRange = [-1 1];
    
    %% Code
    dirListNIfTI = dir(fullfile(figFolder, ['*' groupID '*.nii']));
    % Load IOI matrix of the source image
    % if strcmp(groupID, 'LP')
    IOImat = 'D:\Edgar\OIS_Results\16_02_25,NC01\ROI\LPF\FiltNDown\GLM\corrMap\IOI.mat';
    % else
    %     fprintf('No IOI matrix found for %s.\n', groupID);
    %     return
    % end
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
    nColorLevels = 256;
    % ------------------------------------------------------------------------------
    % Define matlab batch job with the required fields
    % ------------------------------------------------------------------------------
    job(1).figCmap                                  = ioi_get_colormap('bipolar');     % colormap
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
    % r1 = 1;             % ROI index
    % c1 = 5;             % Color index
    % Main loop
    for iFiles = 1:numel(images2overlay)
        [displayed_pixels{iFiles}, total_pixels(iFiles)] = ioi_overlay_blend_intralipid(...
            IOImat, job, images2overlay{iFiles}, fcMapRange, alphaRange, nColorLevels, r1(iR), c1);
    end
    
    %% Average images
    % reload atlas
    IOImat = 'D:\Edgar\OIS_Results\16_02_25,NC01\ROI\LPF\FiltNDown\GLM\corrMap\IOI.mat';
    for i=1:size(images2overlay,1)
        vol = spm_vol(images2overlay{i});
        I(:,:,i) = spm_read_vols(vol);
    end
    avg_I = nanmean(I,3);
    ioi_save_nifti(avg_I, fullfile(figFolder, sprintf('%s_AVG_R%d_C%d.nii',groupID, r1(iR), c1)), [1 1 1]);
    [displayed_pixels{iFiles}, total_pixels(iFiles)] = ioi_overlay_blend_intralipid(...
        IOImat, job, fullfile(figFolder,sprintf('%s_AVG_R%d_C%d.nii',groupID, r1(iR), c1)), fcMapRange, alphaRange, nColorLevels, r1(iR), c1);
end % ROI loop
% EOF
