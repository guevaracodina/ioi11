%% Script to display correlation maps as overlays
% Loading data
dirName = 'D:\Edgar\Data\IOS_Results\12_06_20,CO04\GLMfcIOS\corrMap';
load(fullfile(dirName,'IOI.mat'))

% Load p-Values
load(IOI.fcIOS.corr.fname)

% Choose seed/ROI
r1                  = 1;
% Choose session
s1                  = 2;
% Choose color
c1                  = 5;

% ------------------------------------------------------------------------------
% Example
% ------------------------------------------------------------------------------
anatomical      = IOI.res.file_anat;
positiveMap     = IOI.fcIOS.corr.corrMapName{r1}{s1, c1};
negativeMap     = IOI.fcIOS.corr.corrMapName{r1}{s1, c1};
colorNames      = fieldnames(IOI.color);
mapRange        = [0.3 1];
titleString     = sprintf('%s seed %d S%d(%s)',IOI.subj_name,r1,s1,colorNames{1+c1});
% Display plots on SPM graphics window
% spm_figure('GetWin', 'Graphics');
% spm_figure('Clear', 'Graphics');
% ------------------------------------------------------------------------------

% Function
h = ioi_overlay_map(anatomical, positiveMap, negativeMap, mapRange, titleString);

% EOF
