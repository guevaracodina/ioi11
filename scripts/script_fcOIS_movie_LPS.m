%% script_fcOIS_movie
% script demonstrating the principle of seed-based functional connectivity
clear; close all; clc;
% % pathName    = 'F:\Edgar\Dropbox\PhD\Thesis\Figures\fcOIS_movie';
% pathName    = 'C:\Users\Ramón\Desktop\Edgar\movie';
% fileName    = 'rs_OIS_LPS_movie';
% fileNameROI = 'rs_OIS_ROI_LPS_movie';
% vidWidth    = 1920/2;
% vidHeight   = 1080/2;
% frameRate   = 60;
% saveVideo   = true;
% 
% %% Load data
% IOImat = 'D:\Edgar\\OIS_Results\16_07_08,NC09\ROI\LPF\FiltNDown\GLM\corrMap\IOI.mat';
% load(IOImat);
% % Session (s1 = 1)
% s1 = 1;
% % ROI/seed (S_R = 5, S_L = 6)
% r1 = 5;
% % Contrast (HbO2 = 5, HbR = 6)
% c1 = 6;
% % Labels font size
% labelFont = 12;
% % Functional map
% % images2overlay = IOI.fcIOS.corr.corrMapName{r1}{s1, c1};
% % Folder where to save the images
% % figFolder = 'D:\Edgar\Documents\Dropbox\Docs\Thesis\Figures\seed_based_principle';
% figFolder = pathName;
% 
% % Read maps
% % vol = spm_vol(IOI.res.file_anat);
% vol = spm_vol('D:\Edgar\OIS_Results\averaged_maps\AVG_Atlas.hdr');
% anatomical = rot90(spm_read_vols(vol), 1);
% 
% % vol = spm_vol(IOI.fcIOS.corr.corrMapName{r1}{c1});
% vol = spm_vol(IOI.fcIOS.SPM.fname{s1, c1});
% fcMap = spm_read_vols(vol);
% ioi_text_waitbar(0, 'Please wait...');
% for iFrames=1:size(fcMap,3)
%     fcMap(:,:,iFrames) = rot90(squeeze(fcMap(:,:,iFrsprintfames)), 1);
%     ioi_text_waitbar(iFrames/size(fcMap,3), sprintf('Processing frame %d from %d',iFrames, size(fcMap,3)));
% end
% ioi_text_waitbar('Clear');
% 
% % vol = spm_vol(IOI.fcIOS.mask.fname);
% vol = spm_vol('D:\Edgar\OIS_Results\averaged_maps\16_02_25,NC01_anat_brainmask.nii');
% brainMask = rot90(spm_read_vols(vol), 1);
% % brainMask = spm_read_vols(vol);
% 
% % Range of values to map to the full range of colormap: [minVal maxVal]
% % fcMapRange = [-1 1];
% fcMapRange = [min(fcMap(:)) max(fcMap(:))];
% % Range of values to map to display non-transparent pixels: [minVal maxVal]
% % alphaRange = [0 fcMapRange(2) fcMapRange(1) 0];
% alphaRange = [];
% % ------------------------------------------------------------------------------
% % Define anonymous functions for affine transformations
% % ------------------------------------------------------------------------------
% rotx = @(theta) [1 0 0 0; 0 cos(theta) -sin(theta) 0; 0 sin(theta) cos(theta) 0; 0 0 0 1];
% roty = @(theta) [cos(theta) 0 sin(theta) 0; 0 1 0 0; -sin(theta) 0 cos(theta) 0; 0 0 0 1];
% rotz = @(theta) [cos(theta) -sin(theta) 0 0; sin(theta) cos(theta) 0 0; 0 0 1 0; 0 0 0 1];
% translate = @(a,b) [1 0 a 0; 0 1 b 0; 0 0 1 0; 0 0 0 1];
% % ------------------------------------------------------------------------------
% 
% % ------------------------------------------------------------------------------
% % Define matlab batch job with the required fields
% % ------------------------------------------------------------------------------
% job(1).figCmap                                  = jet(256);     % colormap
% job(1).figIntensity                             = 1;            % [0 - 1]
% job(1).transM                                   = rotz(pi);     % affine transform
% job(1).figSize                                  = [3 3];        % inches
% job(1).figRes                                   = 300;          % in dpi
% job(1).drawCircle(1).drawCircle_On(1).circleLW  = 03;          % line width
% job(1).drawCircle(1).drawCircle_On(1).circleLS  = '-';          % line style
% job(1).drawCircle(1).drawCircle_On(1).circleEC  = 'r';          % line color
% job.parent_results_dir{1}                       = fullfile(figFolder,'overlay');
% job.generate_figures                            = true;         % display figure
% job.save_figures                                = false;         % save figure
% % ------------------------------------------------------------------------------

%% Load data from mat file
load('C:\Users\Ramón\Desktop\Edgar\movie_data\rs_OIS_LPS_movie_data.mat');
load('C:\Users\Ramón\Desktop\Edgar\movie_data\rs_OIS_LPS_movie_ROI.mat');
pathName = 'C:\Users\Ramón\Desktop\Edgar\movie';

%% Display overlay map
% Frames to show
idx = 1:size(fcMap, 3);

textOffset = [40 30];
labelFont = 6;
job.drawCircle.drawCircle_On.circleLW = 1;
scaleFactor=2;
fudgeFactor = 2.5352; % Necessary to scale panel A to [960 540]
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
        scaleFactor*IOI.res.ROI{1, 5}.radius   scaleFactor*IOI.res.ROI{1, 5}.radius];
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
    seedDims = [IOI.res.ROI{1, 6}.center(1)-80  ...
        size(fcMap,2)-IOI.res.ROI{1, 6}.center(2)/4   ...
        scaleFactor*IOI.res.ROI{1, 6}.radius   scaleFactor*IOI.res.ROI{1, 6}.radius];
    rectangle('Position',seedDims,...
        'Curvature',[1,1],...
        'LineWidth',job.drawCircle.drawCircle_On.circleLW,...
        'LineStyle',job.drawCircle.drawCircle_On.circleLS,...
        'EdgeColor','b');
    % Specify window units
    set(h,'inverthardcopy','off')
    set(h, 'units', 'inches')
%     set(h, 'position', [1 1 vidWidth vidHeight])
%     set(h, 'paperPosition', [1 1 vidWidth vidHeight])
    set(h, 'position', [1 1 fudgeFactor*vidWidth/job.figRes fudgeFactor*vidHeight/job.figRes])
    set(h, 'paperPosition', [1 1 fudgeFactor*vidWidth/job.figRes fudgeFactor*vidHeight/job.figRes])
    
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
    %     saveas(h, fullfile(pathName,[fileName sprintf('_%04d.png', iFrames)]), 'png')
    print(h, '-dpng', ...
        fullfile(pathName,[fileName sprintf('_%04d.png', iFrames)]),...
        sprintf('-r%d',job.figRes));
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

t = linspace(0,(numel(idx)-1)/IOI.fcIOS.filtNdown.fs,numel(idx));
% Axis limits
yLim = [0.9*min([ROIregress{7}{s1,c1}(:); ROIregress{7}{s1,c1}(:)]) 1.1*max([ROIregress{7}{s1,c1}(:); ROIregress{8}{s1,c1}(:)])];
xLim = [t(1) t(end)];
% Font sizes
axisFont = 14;
axisLabelFont = 16;
dottedLineWidth = 2;
figName = 'fcOIS_ROI_timecourse';
h = figure; set(h, 'color', 'w')
vidWidth = 1920;
vidHeight = 1080/2;
fprintf('corr(%s, %s) = %0.4f \n',IOI.ROIname{5}, IOI.ROIname{6}, corr(ROIregress{7}{s1,c1}(idx)', ROIregress{8}{s1,c1}(idx)'))
% fprintf('corr(%s, %s) = %0.4f \n',IOI.ROIname{1}, IOI.ROIname{6}, corr(ROIregress{1}{s1,c1}(idx)', ROIregress{3}{s1,c1}(idx)'))
% Folder where to save the images
figFolder = pathName;
job.figSize = [6.5 3];
job.figRes = 300;
% Specify window units
set(h, 'units', 'pixels')
% Change figure and paper size
set(h, 'Name', figName)
set(h, 'position', [100 150 vidWidth vidHeight])
%  Begin frame grabbing
ioi_text_waitbar(0, 'Please wait...');
if saveVideo
    writerObj = VideoWriter(fullfile(pathName,[figName '.mp4']), 'MPEG-4');
    writerObj.FrameRate = frameRate;
    open(writerObj);
end
for iFrames = idx
    % S_R, , r1 = 5
    r1 = 7;
    plot(t(1:iFrames-idx(1)+1), ROIregress{r1}{s1,c1}(idx(1):iFrames),'-','color',[235 101 35]/256, 'LineWidth', dottedLineWidth)
    hold on
    % S_L, r1 = 6
    r1 = 8;
    plot(t(1:iFrames-idx(1)+1), ROIregress{r1}{s1,c1}(idx(1):iFrames),'--','color','b', 'LineWidth', dottedLineWidth)
%     r1 = 3;
%     plot(t(1:iFrames-idx(1)+1), 10*ROIregress{r1}{s1,c1}(idx(1):iFrames),'-.','color',[0 0 1], 'LineWidth', dottedLineWidth)
    hold off
    axis tight
    ylim(yLim)
    xlim(xLim)
    box off
    set(gca,'FontSize',axisFont)
    xlabel('Time [s]', 'FontSize', axisLabelFont)
    ylabel('\DeltaHbR [%]', 'FontSize', axisLabelFont)
    hl = legend({'S_R' 'S_L'}, 'Location', 'EastOutside', 'FontSize', axisFont);
    set(hl, 'Box', 'off')
    set(hl, 'Color', 'none')
    if saveVideo
%         videoOut(iFrames-idx(1)+1) = getframe(h);
        frame = getframe(h);
        writeVideo(writerObj,frame);
    end
    saveas(gcf, fullfile(pathName,[figName sprintf('_%04d.png', iFrames)]), 'bmp')
    if iFrames == idx(end)
        set(h, 'PaperPosition', [100 150 vidWidth vidHeight])
        print(h, '-dpng', ...
            fullfile(pathName,[figName '_last_frame.png']),...
            sprintf('-r%d',job.figRes));
    end
%     drawnow
%     pause(0.1);
ioi_text_waitbar(iFrames/numel(idx), sprintf('Processing frame %d from %d',...
        iFrames, numel(idx)));
end
if saveVideo
    close(writerObj);
end
ioi_text_waitbar('Clear');
fprintf('Processing video done!...\n');

%% Display movie
% h = figure;
% movie(h, videoOut, 1, frameRate);


%% Create full-HD frame
clear; close all; clc
% figFolder = 'C:\Users\Ramón\Desktop\Edgar\movie';
% dataFolder = 'C:\Users\Ramón\Desktop\Edgar\movie_data';
figFolder = 'D:\Edgar\OIS_Results\movie';
dataFolder = 'D:\Edgar\OIS_Results\movie_data';
figA_name = 'rs_OIS_LPS_movie_';
figB_name = 'fcOIS_ROI_timecourse_';
figC_name = 'rs_OIS_LPS_movie_caption_';
fileExt = '.png';
% Number of frames
idx = 1:2000;
% Preallocate
fullFrameHD = uint8(zeros([1080 1920 3]));
fileName = 'fcOIS_fullHD';
vidWidth = 1920;
vidHeight = 1080;
frameRate = 60;
saveVideo = true;
%  Begin frame grabbing
ioi_text_waitbar(0, 'Please wait...');
if saveVideo
    writerObj = VideoWriter(fullfile(figFolder,[fileName '.mp4']), 'MPEG-4');
    writerObj.FrameRate = frameRate;
    open(writerObj);
end
for iFrames=idx,
    % Read HbR maps
    currFig = imread(fullfile(figFolder,sprintf('%s%04d%s', figA_name, iFrames, fileExt)));
    % Fill panel A (Top Left)
    fullFrameHD(1:vidHeight/2, 1:vidWidth/2, :) = currFig;
    
    % Read ROI time course
    [currFig, cMap] = imread(fullfile(figFolder,sprintf('%s%04d%s', figB_name, iFrames, fileExt)));
    currFig = im2uint8(ind2rgb(currFig, cMap));
    % Orange mask (S_R)
    imMaskOrange = false(size(currFig));
    imMaskOrange(:,:,1) = currFig(:,:,1)==224;
    imMaskOrange(:,:,2) = currFig(:,:,2)==96;
    imMaskOrange(:,:,3) = currFig(:,:,3)==64;
    imMaskANDOrange = ~(squeeze(imMaskOrange(:,:,1)) & squeeze(imMaskOrange(:,:,2)) & squeeze(imMaskOrange(:,:,3)));
    % Blue mask (S_L)
    imMaskBlue = false(size(currFig));
    imMaskBlue(:,:,1) = currFig(:,:,1)==0;
    imMaskBlue(:,:,2) = currFig(:,:,2)==0;
    imMaskBlue(:,:,3) = currFig(:,:,3)==255;
    imMaskANDBlue = ~(squeeze(imMaskBlue(:,:,1)) & squeeze(imMaskBlue(:,:,2)) & squeeze(imMaskBlue(:,:,3)));
    % Full mask
    imMask = repmat(imMaskANDBlue & imMaskANDOrange, [1 1 3]);
    % Invert all other colors different than orange or blue
    currFig(imMask) = imcomplement(currFig(imMask));
    % Fill panel B (Bottom)
    fullFrameHD(vidHeight/2+1:end, :, :) = currFig;
    
    if iFrames > 0 && iFrames <= 333
        % Read captions
        [currFig, cMap] = imread(fullfile(dataFolder,sprintf('%sA%s', figC_name, fileExt)));
        currFig = im2uint8(ind2rgb(currFig, cMap));
    end
    if iFrames > 333 && iFrames <= 666
        % Read captions
        currFig = imread(fullfile(dataFolder,sprintf('%sB%s', figC_name, fileExt)));
    end
    if iFrames > 666 && iFrames <= 1000
        % Read captions
        currFig = imread(fullfile(dataFolder,sprintf('%sC_L%s', figC_name, fileExt)));
    end
    if iFrames > 1000 && iFrames <= 1333
        % Read captions
        currFig = imread(fullfile(dataFolder,sprintf('S_L_map%s', fileExt)));
    end
    if iFrames > 1333 && iFrames <= 1666
        % Read captions
        currFig = imread(fullfile(dataFolder,sprintf('%sC_R%s', figC_name, fileExt)));
    end
    if iFrames > 1666 && iFrames <= 2000
        % Read captions
        currFig = imread(fullfile(dataFolder,sprintf('S_R_map%s', fileExt)));
    end
    
    % Fill panel C (Top Right)
    fullFrameHD(1:vidHeight/2, vidWidth/2+1:end, :) = currFig;
    
    % Create video or capture frames
    if saveVideo
        % Show HD movie frame (not necessary)
%         h = imshow(fullFrameHD);
%         drawnow;
%         set(gca,'nextplot','replacechildren');
%         set(h,'Renderer','zbuffer');
%         frame = getframe(h);
        writeVideo(writerObj,fullFrameHD);
    else
        imwrite(fullFrameHD, ...
            fullfile(figFolder,[fileName sprintf('_%04d.png', iFrames)]));
    end
    ioi_text_waitbar(iFrames/numel(idx), sprintf('Processing full-HD frame %d from %d',...
        iFrames, numel(idx)));
end
if saveVideo
    close(writerObj);
end
ioi_text_waitbar('Clear');
fprintf('Processing full-HD video done!...\n');
% EOF
