% script demonstrating the principle of seed-based functional connectivity
%% Load data
IOImat = 'D:\Edgar\Data\IOS_Carotid_Res\12_10_18,NC09\GLMfcIOS\corrMap\IOI.mat';
load(IOImat);
% Session
s1 = 1;
% ROI/seed
r1 = 7;
% Contrast
c1 = 8;
% Functional map
images2overlay = IOI.fcIOS.corr.corrMapName{r1}{s1, c1};
% Folder where to save the images
figFolder = 'D:\Edgar\Documents\Dropbox\Docs\Thesis\Figures\seed_based_principle';
% Range of values to map to the full range of colormap: [minVal maxVal]
fcMapRange = [-1 1];
% Range of values to map to display non-transparent pixels: [minVal maxVal]
alphaRange = [0.35 1 -1 -0.35];
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
job(1).drawCircle(1).drawCircle_On(1).circleLW  = 0.8;          % line width
job(1).drawCircle(1).drawCircle_On(1).circleLS  = '-';          % line style
job(1).drawCircle(1).drawCircle_On(1).circleEC  = 'r';          % line color
job.parent_results_dir{1}                       = fullfile(figFolder,'overlay');
job.generate_figures                            = true;         % display figure
job.save_figures                                = true;         % save figure
% ------------------------------------------------------------------------------


%% Display overlay map
ioi_overlay_blend(IOImat, job, images2overlay, fcMapRange, alphaRange, 256, r1, c1);

%% Load data for ROI time course
IOImat = 'D:\Edgar\Data\IOS_Carotid_Res\12_10_18,NC09\GLMfcIOS\corrMap\ROI_test\IOI.mat';
load(IOImat)
load(IOI.fcIOS.SPM.fnameROIregress)

%% Plot ROI time courses
close all
% Session
s1 = 1;
% Contrast
c1 = 8;
idx = 114:414;
t = linspace(0,300,numel(idx));
% Font sizes
axisFont = 14;
axisLabelFont = 16;
dottedLineWidth = 2;
h = figure; set(h, 'color', 'w')
r1 = 1;
plot(t, 10*ROIregress{r1}{s1,c1}(idx),'-','color',[1 0 0], 'LineWidth', dottedLineWidth)
hold on
r1 = 2;
plot(t, 10*ROIregress{r1}{s1,c1}(idx),'--','color',[235 101 35]/256, 'LineWidth', dottedLineWidth)
r1 = 3;
plot(t, 10*ROIregress{r1}{s1,c1}(idx),'-.','color',[0 0 1], 'LineWidth', dottedLineWidth)
axis tight
box off
set(gca,'FontSize',axisFont)
xlabel('Time [s]', 'FontSize', axisLabelFont)
ylabel('\DeltaCMRO_2 [%]', 'FontSize', axisLabelFont)
% legend({'M_R' 'M_L' 'V_R'}, 'Location', 'NorthWest', 'FontSize', axisFont)
fprintf('corr(%s, %s) = %0.4f \n',IOI.ROIname{1}, IOI.ROIname{2}, corr(ROIregress{1}{s1,c1}(idx)', ROIregress{2}{s1,c1}(idx)'))
fprintf('corr(%s, %s) = %0.4f \n',IOI.ROIname{1}, IOI.ROIname{3}, corr(ROIregress{1}{s1,c1}(idx)', ROIregress{3}{s1,c1}(idx)'))
% Folder where to save the images
figFolder = 'D:\Edgar\Documents\Dropbox\Docs\Thesis\Figures\seed_based_principle';
job.figSize = [6.5 3];
job.figRes = 300;
% Specify window units
set(h, 'units', 'inches')
% Change figure and paper size
set(h, 'Position', [0.1 0.1 job.figSize(1) job.figSize(2)])
set(h, 'PaperPosition', [0.1 0.1 job.figSize(1) job.figSize(2)])

%% Print seeds time course
% Save as PNG at the user-defined resolution
print(h, '-dpng', ...
    fullfile(figFolder,'seeds_time_course_01'),...
    sprintf('-r%d',job.figRes));
% Return the property to its default
set(h, 'units', 'pixels')

% EOF
