%% script_fcOIS_movie
% script demonstrating the principle of seed-based functional connectivity

% pathName    = 'F:\Edgar\Dropbox\PhD\Thesis\Figures\fcOIS_movie';
pathName    = 'D:\Edgar\OIS_Results\movie';
fileName    = 'rs_OIS_LPS_movie';
fileNameROI = 'rs_OIS_ROI_LPS_movie';
vidWidth    = 1920/2;
vidHeight   = 1080/2;
frameRate   = 60;
saveVideo   = true;

%% Load data
IOImat = 'D:\Edgar\OIS_Results\16_07_08,NC09\ROI\LPF\FiltNDown\GLM\corrMap\IOI.mat';
load(IOImat);
% Session (s1 = 1)
s1 = 1;
% ROI/seed (S_R = 5, S_L = 6)
r1 = 5;
% Contrast (HbO2 = 5, HbR = 6)
c1 = 6;
% Labels font size
labelFont = 12;
% Functional map
% images2overlay = IOI.fcIOS.corr.corrMapName{r1}{s1, c1};
% Folder where to save the images
% figFolder = 'D:\Edgar\Documents\Dropbox\Docs\Thesis\Figures\seed_based_principle';
figFolder = pathName;

% Read maps
vol = spm_vol(IOI.res.file_anat);
anatomical = rot90(spm_read_vols(vol), 1);

% vol = spm_vol(IOI.fcIOS.corr.corrMapName{r1}{c1});
vol = spm_vol(IOI.fcIOS.SPM.fname{s1, c1});
fcMap = spm_read_vols(vol);
ioi_text_waitbar(0, 'Please wait...');
for iFrames=1:size(fcMap,3)
    fcMap(:,:,iFrames) = rot90(squeeze(fcMap(:,:,iFrsprintfames)), 1);
    ioi_text_waitbar(iFrames/size(fcMap,3), ('Processing frame %d from %d',...
        iFrames, size(fcMap,3)));
end
ioi_text_waitbar('Clear');

vol = spm_vol(IOI.fcIOS.mask.fname);
brainMask = rot90(spm_read_vols(vol), 1);
% brainMask = spm_read_vols(vol);

% Range of values to map to the full range of colormap: [minVal maxVal]
% fcMapRange = [-1 1];
fcMapRange = [min(fcMap(:)) max(fcMap(:))];
% Range of values to map to display non-transparent pixels: [minVal maxVal]
% alphaRange = [0 fcMapRange(2) fcMapRange(1) 0];
alphaRange = [];
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
job(1).figSize                                  = [3 3];        % inches
job(1).figRes                                   = 300;          % in dpi
job(1).drawCircle(1).drawCircle_On(1).circleLW  = 03;          % line width
job(1).drawCircle(1).drawCircle_On(1).circleLS  = '-';          % line style
job(1).drawCircle(1).drawCircle_On(1).circleEC  = 'r';          % line color
job.parent_results_dir{1}                       = fullfile(figFolder,'overlay');
job.generate_figures                            = true;         % display figure
job.save_figures                                = false;         % save figure
% ------------------------------------------------------------------------------

%% Display overlay map
% Frames to show
idx = 1:size(fcMap, 3);

textOffset = [40 30];
% h = figure;
% Preallocate structure movIn
% videoOut(1:numel(idx)) = ...
%     struct('cdata', zeros(vidHeight, vidWidth, 3, 'uint8'),...
%     'colormap', []);

%  Begin frame grabbing
ioi_text_waitbar(0, 'Please wait...');
if saveVideo
    writerObj = VideoWriter(fullfile(pathName,[fileName '.mp4']), 'MPEG-4');
    writerObj.FrameRate = frameRate;
    open(writerObj);
end

for iFrames = idx
    %     images2overlay = sprintf('%s,%d',IOI.fcIOS.SPM.fname{s1, c1},iFrames);
    %     ioi_overlay_blend(IOImat, job, images2overlay, fcMapRange, alphaRange, 256, r1, c1);
    [fcMapBlend, h] = ioi_overlay_blend(anatomical, squeeze(fcMap(:,:,iFrames)), ...
        brainMask, fcMapRange, alphaRange, ioi_get_colormap('bipolar'), 1);
    %     hold on
    
    % S_R
    % Add text S_R
    text(IOI.res.ROI{1, 5}.center(1) + textOffset(1), ...
        size(fcMap,2)-IOI.res.ROI{1, 5}.center(2)/4 - textOffset(2),...
        'S_R', 'FontWeight', 'bold', 'FontSize', labelFont, 'color', 'w')
    % Draw Right Somatosensory seed
    seedDims = [IOI.res.ROI{1, 5}.center(1)  ...
        size(fcMap,2)-IOI.res.ROI{1, 5}.center(2)/4   ...
        IOI.res.ROI{1, 5}.radius   IOI.res.ROI{1, 5}.radius];
    rectangle('Position',seedDims,...
        'Curvature',[1,1],...
        'LineWidth',job.drawCircle.drawCircle_On.circleLW,...
        'LineStyle',job.drawCircle.drawCircle_On.circleLS,...
        'EdgeColor',[235 101 35]/256);
    
    % S_L
    % Add text S_L
    text(IOI.res.ROI{1, 6}.center(1) - textOffset(1), ...
        size(fcMap,2)-IOI.res.ROI{1, 6}.center(2)/4 - textOffset(2),...
        'S_L', 'FontWeight', 'bold', 'FontSize', labelFont, 'color', 'w')
    % Draw Left Somatosensory seed
    seedDims = [IOI.res.ROI{1, 6}.center(1)  ...
        size(fcMap,2)-IOI.res.ROI{1, 6}.center(2)/4   ...
        IOI.res.ROI{1, 6}.radius   IOI.res.ROI{1, 6}.radius];
    rectangle('Position',seedDims,...
        'Curvature',[1,1],...
        'LineWidth',job.drawCircle.drawCircle_On.circleLW,...
        'LineStyle',job.drawCircle.drawCircle_On.circleLS,...
        'EdgeColor','b');
    % Specify window units
    set(h, 'units', 'pixels')
    set(h, 'position', [100 150 vidWidth vidHeight])
    
    if saveVideo
        set(gca,'nextplot','replacechildren');
        set(gcf,'Renderer','zbuffer');
    end
    
    % Grab frames
    %     videoOut(iFrames-idx(1)+1) = getframe(h);
    if saveVideo
        frame = getframe(h);
        writeVideo(writerObj,frame);
    end
    %     hold off
%     print(h, '-dpng', ...
%             fullfile(pathName,[fileName sprintf('_%04d.png', iFrames)]));
    saveas(gcf, fullfile(pathName,[fileName sprintf('_%04d.png', iFrames)]), 'bmp')
    if iFrames == idx(end),
        print(h, '-dpng', ...
            fullfile(pathName,[fileName '_last_frame.png']),...
            sprintf('-r%d',job.figRes));
    end
    ioi_text_waitbar(iFrames/numel(idx), sprintf('Processing frame %d from %d',...
        iFrames, numel(idx)));
end
if saveVideo
    close(writerObj);
end
ioi_text_waitbar('Clear');
fprintf('Processing video done!...\n');

%% Load data for ROI time course
% IOImat = 'D:\Edgar\Data\IOS_Carotid_Res\12_10_18,NC09\GLMfcIOS\corrMap\ROI_test\IOI.mat';
load(IOImat)
load(IOI.fcIOS.SPM.fnameROIregress)

%% Plot ROI time courses as a video
close all; clear videoOut

t = linspace(0,numel(idx)-1,numel(idx));
% Axis limits
yLim = [-1.195 1.1];
xLim = [t(1) t(end)];
% Font sizes
axisFont = 14;
axisLabelFont = 16;
dottedLineWidth = 2;
figName = 'fcOIS_ROI_timecourse';
h = figure; set(h, 'color', 'w')
vidWidth = 640;
vidHeight = 480/2;
fprintf('corr(%s, %s) = %0.4f \n',IOI.ROIname{1}, IOI.ROIname{2}, corr(ROIregress{1}{s1,c1}(idx)', ROIregress{2}{s1,c1}(idx)'))
fprintf('corr(%s, %s) = %0.4f \n',IOI.ROIname{1}, IOI.ROIname{3}, corr(ROIregress{1}{s1,c1}(idx)', ROIregress{3}{s1,c1}(idx)'))
% Folder where to save the images
figFolder = pathName;
job.figSize = [6.5 3];
job.figRes = 300;
% Specify window units
set(h, 'units', 'pixels')
% Change figure and paper size
set(h, 'Name', figName)
set(h, 'position', [100 150 vidWidth vidHeight])
% Preallocate structure movIn
videoOut(1:numel(t)) = ...
    struct('cdata', zeros(vidHeight, vidWidth, 3, 'uint8'),...
    'colormap', []);
% Begin frame grabbing
for iFrames = idx
    r1 = 1;
    plot(t(1:iFrames-idx(1)+1), 10*ROIregress{r1}{s1,c1}(idx(1):iFrames),'-','color',[1 0 0], 'LineWidth', dottedLineWidth)
    hold on
    r1 = 2;
    plot(t(1:iFrames-idx(1)+1), 10*ROIregress{r1}{s1,c1}(idx(1):iFrames),'--','color',[235 101 35]/256, 'LineWidth', dottedLineWidth)
    r1 = 3;
    plot(t(1:iFrames-idx(1)+1), 10*ROIregress{r1}{s1,c1}(idx(1):iFrames),'-.','color',[0 0 1], 'LineWidth', dottedLineWidth)
    % axis tight
    ylim(yLim)
    xlim(xLim)
    box off
    set(gca,'FontSize',axisFont)
    xlabel('Time [s]', 'FontSize', axisLabelFont)
    ylabel('\DeltaCMRO_2 [%]', 'FontSize', axisLabelFont)
    hl = legend({'S_R' 'S_L' 'V_R'}, 'Location', 'EastOutside', 'FontSize', axisFont);
    set(hl, 'Box', 'off')
    set(hl, 'Color', 'none')
    videoOut(iFrames-idx(1)+1) = getframe(h);
    if iFrames == idx(end)
        set(h, 'PaperPosition', [100 150 vidWidth vidHeight])
        print(h, '-dpng', ...
            fullfile(pathName,[figName '_last_frame.png']),...
            sprintf('-r%d',job.figRes));
    end
end
fprintf('Processing video done!...\n');

%% Display movie
% h = figure;
% movie(h, videoOut, 1, frameRate);

% EOF
