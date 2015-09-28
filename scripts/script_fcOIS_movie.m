%% script_fcOIS_movie
% script demonstrating the principle of seed-based functional connectivity
% pathName    = 'F:\Edgar\Dropbox\PhD\Thesis\Figures\fcOIS_movie';
pathName    = 'C:\Edgar\OneDrive for Business\Poly\PhD\Thesis\Figures';
fileName    = 'fcOIS_resting_state.avi';
vidWidth    = 640/2;
vidHeight   = 480/2;
frameRate   = 30;
saveVideo   = false;

%% Load data
IOImat = 'D:\Edgar\Data\IOS_Carotid_Res\12_10_18,NC09\GLMfcIOS\corrMap\IOI.mat';
load(IOImat);
% Session (s1 = 1)
s1 = 1;
% ROI/seed (S_R = 7, S_L = 8, V_R = 11)
r1 = 7;
% Contrast (CMRO2 = 8)
c1 = 8;
% Labels font size
labelFont = 12;
% Functional map
images2overlay = IOI.fcIOS.corr.corrMapName{r1}{s1, c1};
% Folder where to save the images
figFolder = 'D:\Edgar\Documents\Dropbox\Docs\Thesis\Figures\seed_based_principle';
% Range of values to map to the full range of colormap: [minVal maxVal]
% fcMapRange = [-1 1];
fcMapRange = [-0.08 0.08];
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
idx = 114:414;
textOffset = [20 15];
h = figure;
% Preallocate structure movIn
videoOut(1:numel(idx)) = ...
    struct('cdata', zeros(vidHeight, vidWidth, 3, 'uint8'),...
           'colormap', []);
       
%  Begin frame grabbing      
% for iFrames = idx
%     images2overlay = sprintf('%s,%d',IOI.fcIOS.SPM.fname{s1, c1},iFrames);
%     ioi_overlay_blend(IOImat, job, images2overlay, fcMapRange, alphaRange, 256, r1, c1);
% %     hold on
%     % Specify window units
%     set(h, 'units', 'pixels')
%     set(h, 'position', [100 150 vidWidth vidHeight])
%     % Add text S_R
%     text(233.8680 - textOffset(1), 195.1770 - textOffset(2), 'S_R', 'FontWeight', 'bold', 'FontSize', labelFont)
%     % Draw Left Somatosensory seed
%     seedDims = [63.7500  195.7781   14.0000   14.0000];
%     rectangle('Position',seedDims,...
%             'Curvature',[1,1],...
%             'LineWidth',job.drawCircle.drawCircle_On.circleLW,...
%             'LineStyle',job.drawCircle.drawCircle_On.circleLS,...
%             'EdgeColor',[235 101 35]/256);
%     % Add text S_L
%     text(seedDims(1) - textOffset(1), seedDims(2) - textOffset(2), 'S_L', 'FontWeight', 'bold', 'FontSize', labelFont)
% 	% Draw Right Visual seed
%     seedDims = [218.8399  - 50  278.1320 + 20  14.0000   14.0000];
%     rectangle('Position',seedDims,...
%             'Curvature',[1,1],...
%             'LineWidth',job.drawCircle.drawCircle_On.circleLW,...
%             'LineStyle',job.drawCircle.drawCircle_On.circleLS,...
%             'EdgeColor','b');
%     % Add text V_R
% 	text(seedDims(1) - textOffset(1), seedDims(2) - textOffset(2), 'V_R', 'FontWeight', 'bold', 'FontSize', labelFont)
%     videoOut(iFrames-idx(1)+1) = getframe(h);
% %     hold off
% if iFrames == idx(end)
%     print(h, '-dpng', ...
%         fullfile(pathName,[fileName '_last_frame.png']),...
%         sprintf('-r%d',job.figRes));
% end
% end
% fprintf('Processing video done!...\n');

%% Save processed video file
% if saveVideo
%     fprintf('Saving video as %s...\n',fullfile(pathName,fileName));
%     movie2avi(videoOut, fullfile(pathName,fileName), 'fps', frameRate);
%     fprintf('Saving video done!...\n');
% end

%% Load data for ROI time course
IOImat = 'D:\Edgar\Data\IOS_Carotid_Res\12_10_18,NC09\GLMfcIOS\corrMap\ROI_test\IOI.mat';
load(IOImat)
load(IOI.fcIOS.SPM.fnameROIregress)

%% Plot ROI time courses as a video
close all; clear videoOut
% Session
s1 = 1;
% Contrast
c1 = 8;
idx = 114:414;
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

%% Save processed video file
if saveVideo
    fileName        = fullfile(pathName,figName);
    fprintf('Saving video as %s...\n',fileName);
    movie2avi(videoOut, fileName, 'fps', frameRate);
    fprintf('Saving video done!...\n');
end

EOF
