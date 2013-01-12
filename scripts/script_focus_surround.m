%% script_focus_surround
% Recordings in the focus and surround
% IOI matrix
ioiMat = 'E:\Edgar\Data\IOS_Results\12_08_24,EG03\FiltNDown\GLMfcIOS\corrMap\ROI\Series\IOI.mat';
load(ioiMat)

% Read anatomical NIFTI file
vol = spm_vol(IOI.res.file_anat);
im_anat = spm_read_vols(vol);
% Read brainmask file
vol = spm_vol(IOI.fcIOS.mask.fname);
im_mask = spm_read_vols(vol);

%% Display masked image
h = figure; set(gcf,'color','w')
imagesc(im_anat.*im_mask);
axis image;
colormap(gray(256));
set(gca, 'xTick',[]); set(gca, 'yTick',[])
% imcontrast; % Adjust manually and the print

%% Print anatomical image
job.figSize = [2 2];
job.figRes = 1200;
% Specify window units
set(h, 'units', 'inches')
% Change figure and paper size
set(h, 'Position', [0.1 0.1 job.figSize(1) job.figSize(2)])
set(h, 'PaperPosition', [0.1 0.1 job.figSize(1) job.figSize(2)])

% Save as PNG at the user-defined resolution
print(h, '-dpng', ...
    fullfile('D:\Edgar\Documents\Dropbox\Docs\Epilepsy\figs', 'focus_surround_anat'),...
    sprintf('-r%d',job.figRes));
% Return the property to its default
set(h, 'units', 'pixels')

%% Load ROIs images
% Read focus ROI file
vol = spm_vol(IOI.res.ROI{1,13}.fname);
im_ROI_focus = spm_read_vols(vol);
% Read focus ROI file
vol = spm_vol(IOI.res.ROI{1,8}.fname);
im_ROI_surround = spm_read_vols(vol);

%% Display ROIs
h2 = figure; set(gcf,'color','w')
imagesc(im_ROI_focus);
axis image;
colormap(gray(256));
set(gca, 'xTick',[]); set(gca, 'yTick',[])

h3 = figure; set(gcf,'color','w')
imagesc(im_ROI_surround);
axis image;
colormap(gray(256));
set(gca, 'xTick',[]); set(gca, 'yTick',[])

%% Print ROIs
job.figSize = [2 2];
job.figRes = 1200;
% Specify window units
set(h2, 'units', 'inches')
% Change figure and paper size
set(h2, 'Position', [0.1 0.1 job.figSize(1) job.figSize(2)])
set(h2, 'PaperPosition', [0.1 0.1 job.figSize(1) job.figSize(2)])

% Save as PNG at the user-defined resolution
print(h2, '-dpng', ...
    fullfile('D:\Edgar\Documents\Dropbox\Docs\Epilepsy\figs', 'focus_ROI'),...
    sprintf('-r%d',job.figRes));
% Return the property to its default
set(h2, 'units', 'pixels')

% Specify window units
set(h3, 'units', 'inches')
% Change figure and paper size
set(h3, 'Position', [0.1 0.1 job.figSize(1) job.figSize(2)])
set(h3, 'PaperPosition', [0.1 0.1 job.figSize(1) job.figSize(2)])

% Save as PNG at the user-defined resolution
print(h3, '-dpng', ...
    fullfile('D:\Edgar\Documents\Dropbox\Docs\Epilepsy\figs', 'surround_ROI'),...
    sprintf('-r%d',job.figRes));
% Return the property to its default
set(h3, 'units', 'pixels')


%% Load ROI 13 time traces
% ioiMat = 'E:\Edgar\Data\IOS_Results\12_08_24,EG03\FiltNDown\GLMfcIOS\corrMap\ROI\Series\IOI.mat';
ioiMat = 'E:\Edgar\Data\IOS_Results\12_09_28,EG09\Onsets\GLMfcIOS\corrMap\ROI\Series\IOI.mat';
load(ioiMat)
% Define session and color to use later
s1 = 4;
% time limits (see script_seizure_duration.m)
tBegin = 271;
tEnd = 432;
c1 = 5;
% ROI
r1 = 13;
load(IOI.ROI.ROIfname)
roi13=ROI{r1};

%% Load ROI 8 time trace
% load('E:\Edgar\Data\IOS_Results\12_08_24,EG03\FiltNDown\IOI.mat')
ioiMat = 'E:\Edgar\Data\IOS_Results\12_09_28,EG09\Onsets\IOI.mat';
load(ioiMat)
load(IOI.ROI.ROIfname)
ROI{r1} = roi13; clear roi13
% ROI time vector
tROI = (0:numel(ROI{r1}{s1, c1})-1)*IOI.dev.TR;

%% Load LFP data
load(IOI.res.el2{s1})
% LFP time vector
tEl = (0:numel(el2)-1)'/IOI.res.sfel;
tBefore = 5;
tAfter = 5;
tLim = [tBegin-tBefore tEnd+tAfter];
idx1 = find(tROI>(tBegin-tBefore),1,'first');
idx2 = find(tROI>tBegin,1,'first');

%% Display 
h = figure;
set(gcf,'color','w')
% ROI
subplot(211)

c1 = 5;
r1 = 13;
baseValue = mean(ROI{r1}{s1, c1}(idx1:idx2));
plot(tROI, (ROI{r1}{s1, c1}-baseValue)./baseValue,'r-')
hold on
r1 = 8;
baseValue = mean(ROI{r1}{s1, c1}(idx1:idx2));
plot(tROI, (ROI{r1}{s1, c1}-baseValue)./baseValue,'r--')

c1 = 6;
r1 = 13;
baseValue = mean(ROI{r1}{s1, c1}(idx1:idx2));
plot(tROI, (ROI{r1}{s1, c1}-baseValue)./baseValue,'b-')
r1 = 8;
baseValue = mean(ROI{r1}{s1, c1}(idx1:idx2));
plot(tROI, (ROI{r1}{s1, c1}-baseValue)./baseValue,'b--')

c1 = 7;
r1 = 13;
baseValue = mean(ROI{r1}{s1, c1}(idx1:idx2));
plot(tROI, (ROI{r1}{s1, c1}-baseValue)./baseValue,'k-')
r1 = 8;
baseValue = mean(ROI{r1}{s1, c1}(idx1:idx2));
plot(tROI, (ROI{r1}{s1, c1}-baseValue)./baseValue,'k--')


c1 = 8;
r1 = 13;
baseValue = mean(ROI{r1}{s1, c1}(idx1:idx2));
plot(tROI, (ROI{r1}{s1, c1}-baseValue)./baseValue,'c-')
r1 = 8;
baseValue = mean(ROI{r1}{s1, c1}(idx1:idx2));
plot(tROI, (ROI{r1}{s1, c1}-baseValue)./baseValue,'c--')

legend({'HbO_2 focus' 'HbO_2 surround' 'HbR focus' 'HbR surround' ...
    'CBF focus' 'CBF surround' 'CMRO_2 focus' 'CMRO_2 surround'})
xlim(tLim)

% LFP
subplot(212)
plot(tEl, el2, 'k-')
xlim(tLim)


