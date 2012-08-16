%% Script to display correlation maps as overlays
% Loading data
dirName = 'D:\Edgar\Data\IOS_Results\12_06_20,CO04\GLMfcIOS\corrMap\';
load(fullfile(dirName,'IOI.mat'))

% Choose seed/ROI
r1 = 4;
% Choose session
s1 = 1;
% Choose color
c1 = 5;
% Necessary to get rid of old overlay objects
clear imagesOverlay
% Make cell with file names
imagesOverlay{1,1} = IOI.res.file_anat;                     % Anatomical image
imagesOverlay{2,1} = IOI.fcIOS.corr.corrMapName{r1}{s1, c1};% Pos. correlation
imagesOverlay{3,1} = IOI.fcIOS.corr.corrMapName{r1}{s1, c1};% Neg. correlation
figure(999);
% Create overlay object
slObj = slover(char(imagesOverlay));

%% Set overlay parameters
slObj.slices = 1;                           % For IOI images only 1 slice (2D)
                                            % Automatic range for image 1 (anatomy)
slObj.img(2).range = [0.2 1]';              % Range for positive map
slObj.img(3).range = -slObj.img(2).range;   % Same range for negative map
slObj.img(1).type = 'truecolor';            % Anatomical image
slObj.img(2).type = 'split';                % Pos. map
slObj.img(3).type = 'split';                % Neg. map
slObj.img(1).cmap = gray(256);              % Colormap for anatomy
slObj.img(2).cmap = hot(256);               % Colormap for positive map
slObj.img(3).cmap = winter(256);            % Colormap for negative map
slObj.cbar = [2 3];                         % Plot colorbars for images 2 & 3
slObj.area.valign = 'top';                  % Vertical alignment
slObj.area.halign = 'center';               % Horizontal alignment
slObj.img(1).prop = 1;                      % Proportion of intensity for anatomy
slObj.img(2).prop = 1;                      % Proportion of intensity for positive map
slObj.img(3).prop = 1;                      % Proportion of intensity for negative map
slObj.img(1).outofrange =  {0 255};         % Behavior for image values out of range 
slObj.img(2).outofrange =  {0 255};
slObj.img(3).outofrange =  {0 255};
slObj.labels = 'none';                      % No labels on this slice
% Color names
colorNames = fieldnames(IOI.color);
set(slObj.figure,'Name',sprintf('%s seed %d S%d(%s)',IOI.subj_name,r1,s1,colorNames{1+c1}));

% Redraw the object (e.g. after window maximization)
slObj = paint(slObj);

