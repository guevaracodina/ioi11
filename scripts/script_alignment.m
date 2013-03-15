function script_alignment
%% script_alignment
clc;
% Always 1 session
s1 = 1;
% Color index (5=HbO, 6=HbR, 7=CBF, 8=CMRO2)
c1 = 5;
% ROI index
r1 = 6;
% group ID string
groupID = 'CC';

%% Load manually aligned images with ImageJ plugin TurboReg
clear imDataArray
topDir = 'D:\Edgar\Data\IOS_Carotid_Res\alignment';
currentDir = sprintf('%s_R%02dC%02d', groupID, r1, c1);
[dirAlignment, sts] = cfg_getfile([1 1],'dir','Select folder',{fullfile(topDir,currentDir)}, topDir, '.*');
dirListNIfTI = dir(fullfile(topDir,[currentDir filesep '*.nii']));
dirListAnalyze = dir(fullfile(topDir,[currentDir filesep '*.img']));
images2align = [struct2cell(dirListAnalyze) struct2cell(dirListNIfTI)];
images2align = images2align(1,:)';
images2align = cellfun(@(x) fullfile(dirAlignment{1}, [x ',1']), images2align, 'UniformOutput', false);
[images2align, sts] = cfg_getfile([1 Inf],'image','Select images',images2align, dirAlignment{1}, '.*');
V = spm_vol(images2align);
for iVols = 1:numel(V),
    [imData, imXYZ] = spm_read_vols(V{iVols});
    % Set aligned IMG/HDR images to the orientation of source image
    if ndims(imData) ~= 2
    	imData = squeeze(imData(:,:,1));
        imData = fliplr(imData);
    end
    imDataArray(:,:,iVols) = imData;
end
% Average seed-based correlation maps.
imAvg = median(imDataArray, 3);

%% Display average image
figFolder = 'D:\Edgar\Documents\Dropbox\Docs\Carotid\Figures\aligned';
% h = figure; imagesc(imAvg, [-1 1]); axis image; colorbar
figName = fullfile(figFolder,[currentDir '_avg']);

%% Prepare job
% ------------------------------------------------------------------------------
% Define anonymous functions for affine transformations
% ------------------------------------------------------------------------------
rotx = @(theta) [1 0 0 0; 0 cos(theta) -sin(theta) 0; 0 sin(theta) cos(theta) 0; 0 0 0 1];
roty = @(theta) [cos(theta) 0 sin(theta) 0; 0 1 0 0; -sin(theta) 0 cos(theta) 0; 0 0 0 1];
rotz = @(theta) [cos(theta) -sin(theta) 0 0; sin(theta) cos(theta) 0 0; 0 0 1 0; 0 0 0 1];
translate = @(a,b) [1 0 a 0; 0 1 b 0; 0 0 1 0; 0 0 0 1];
% ------------------------------------------------------------------------------

job(1).figCmap                                  = jet(256);
job(1).figRange                                 = [-1 1];
job(1).figIntensity                             = 0.8;
job(1).figAlpha                                 = 1;
job(1).transM                                   = rotz(pi);
job(1).figSize                                  = [1.5 1.5];
job(1).figRes                                   = 300;
job(1).drawCircle(1).drawCircle_On(1).circleLW  = 0.8;
job(1).drawCircle(1).drawCircle_On(1).circleLS  = '-';
job.parent_results_dir{1}                       = figFolder;

%% Overlay slover
% Load IOI matrix of the source image
if strcmp(groupID, 'CC')
    load('D:\Edgar\Data\IOS_Carotid_Res\12_10_19,CC10\GLMfcIOS\corrMap\IOI.mat');
elseif strcmp(groupID, 'NC')
    load('D:\Edgar\Data\IOS_Carotid_Res\12_10_18,NC09\GLMfcIOS\corrMap\IOI.mat');
else
    fprintf('No IOI matrix found for %s.\n', groupID);
    return
end
if IOI.res.shrinkageOn
    vx = [IOI.res.shrink_x IOI.res.shrink_y 1];
else
    vx = [1 1 1];
end
% Read brain mask
maskVol = spm_vol(IOI.fcIOS.mask.fname);
brainMask = spm_read_vols(maskVol);
if size(brainMask,1)~= size(imAvg,1)|| size(brainMask,2)~= size(imAvg,2)
    brainMask = ioi_MYimresize(brainMask, [size(imAvg,1) size(imAvg,2)]);
end
% Make NaN all the transparent pixels
imAvg(~brainMask) = NaN;

% Save NIfTI to use slover
ioi_save_nifti(imAvg, [figName '.nii'], vx);

% Get parameters for overlay
anatomical      = IOI.res.file_anat;
corrMap         = [figName '.nii'];
% Seed annotation dimensions the lower left corner of the bounding rectangle
% at the point seedX, seedY (the + sign here is due to image rotation)
seedX = IOI.res.ROI{r1}.center(2) + IOI.res.ROI{r1}.radius;
seedY = IOI.res.ROI{r1}.center(1) + IOI.res.ROI{r1}.radius;
% Seed width
seedW = 2*IOI.res.ROI{r1}.radius;
% Seed height
seedH = 2*IOI.res.ROI{r1}.radius;
% Change seed circle size if shrunk
if isfield(IOI.res,'shrinkageOn')
    if IOI.res.shrinkageOn == 1
        seedW = seedW * IOI.res.shrink_x;
        seedH = seedH * IOI.res.shrink_y;
    end
end
internal_overlay_map(anatomical, corrMap,  job, [currentDir '_avg'], [seedX seedY seedW seedH]);
colorNames = fieldnames(IOI.color);
fprintf('Average seed-based correlation map done! Group: %s, R%02d, (%s)\n', groupID, r1, colorNames{c1+1});
end

function [h, varargout] = internal_overlay_map(anatomical, positiveMap,  job, titleString, seedDims)
% Creates an overlay image composed of positive & negative functional maps (any
% contrast, correlation, z, p-value, etc.) onto an anatomical image.
% The images must be saved as NIfTI (.nii) files and they can have different
% sizes. This function makes use of slover and paint functions of SPM8.
% The positive map is plotted in hot colormap, the negative map is cold and the
% anatomical is grayscale by default. The image has the following orientation:
% 
%         Rostral
%       |________
%       |
%       |
%  Left | 
%       |        
% 
% SYNTAX
% h = ioi_overlay_map(anatomical,positiveMap,negativeMap,mapRange,titleString)
% INPUTS
% anatomical    NIfTI (.nii) filename with anatomical image for background.
% positiveMap   NIfTI (.nii) filename with positive functional map on the
%               foreground.
% job           Matlab batch job 
% titleString   String with the title to be displayed.
% seedDims      Seed dimensions [seedX seedY seedW seedH]
% OUTPUT
% h             Handle to the figure
% slObj         [OPTIONAL] slover object
%_______________________________________________________________________________
% Copyright (C) 2012 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________


% Necessary to get rid of old overlay objects
clear imagesOverlay
% Make cell with file names
imagesOverlay{1,1}  = anatomical;           % Anatomical image
imagesOverlay{2,1}  = positiveMap;          % Pos. correlation
% imagesOverlay{3,1}  = negativeMap;          % Neg. correlation

% Handle of current figure;
h = figure(999);
set(h,'color','k')

% Get range
mapRange        = {job.figRange};

% Make mapRange a column vector
if size(mapRange, 1) == 1
    mapRange = mapRange';
    mapRange = cellfun(@(x) x', mapRange, 'UniformOutput', false);
end

% ------------------------------------------------------------------------------
% Define anonymous functions for affine transformations
% ------------------------------------------------------------------------------
rotx = @(theta) [1 0 0 0; 0 cos(theta) -sin(theta) 0; 0 sin(theta) cos(theta) 0; 0 0 0 1];
roty = @(theta) [cos(theta) 0 sin(theta) 0; 0 1 0 0; -sin(theta) 0 cos(theta) 0; 0 0 0 1];
rotz = @(theta) [cos(theta) -sin(theta) 0 0; sin(theta) cos(theta) 0 0; 0 0 1 0; 0 0 0 1];
translate = @(a,b) [1 0 a 0; 0 1 b 0; 0 0 1 0; 0 0 0 1];
% ------------------------------------------------------------------------------

% Create overlay object
slObj = slover(char(imagesOverlay));

slObj.slices = 1;                           % For IOI images only 1 slice (2D)
                                            % Automatic range for image 1 (anatomy)
slObj.img(2).range = mapRange{1};           % Range for positive map

slObj.img(1).type = 'truecolour';           % Anatomical image
slObj.img(2).type = 'truecolour';           % Functional map

slObj.img(1).cmap = gray(256);              % Colormap for anatomy
slObj.img(2).cmap = job.figCmap;            % Colormap for functional map

% slObj.cbar = [2 3];                         % Plot colorbars for images 2 & 3
slObj.area.valign = 'middle';               % Vertical alignment
slObj.area.halign = 'center';               % Horizontal alignment

slObj.img(1).prop = job.figIntensity;       % Proportion of intensity for anatomy
slObj.img(2).prop = job.figAlpha;           % Proportion of intensity for positive map

slObj.img(1).outofrange =  {0 255};         % Behavior for image values out of range 
slObj.img(2).outofrange =  {0 255};

slObj.labels = 'none';                      % No labels on this slice

% Apply affine transformation
slObj.transform = job.transM*translate(-slObj.img(1).vol.dim(1),-slObj.img(1).vol.dim(2));               % Apply affine transformation

% Change figure name
set(slObj.figure,'Name',titleString);

% Pass the slover object as output
varargout{1} = slObj;

% Specify window units
set(h, 'units', 'inches')
% Change figure and paper size
set(h, 'Position', [0.1 0.1 job.figSize(1) job.figSize(2)])
set(h, 'PaperPosition', [0.1 0.1 job.figSize(1) job.figSize(2)])
% Refresh figure
slObj = paint(slObj);

% Seed positions and sizes will be shown with black circles
if isfield(job.drawCircle,'drawCircle_On')
    figure(h);
    % New seeds coordinates (after affine transformation)
    newSeedCoord = job.transM*[seedDims(1) seedDims(2) 0 1]' + [slObj.img(1).vol.dim(2) slObj.img(1).vol.dim(1) 0 0]';
    seedDims(1:2) = newSeedCoord(2:-1:1)';
    % Display ROI
    rectangle('Position',seedDims,...
        'Curvature',[1,1],...
        'LineWidth',job.drawCircle.drawCircle_On.circleLW,...
        'LineStyle',job.drawCircle.drawCircle_On.circleLS);
end

% Save as PNG at the user-defined resolution
print(h, '-dpng', ...
    fullfile(job.parent_results_dir{1}, titleString),...
    sprintf('-r%d',job.figRes));

% Return the property to its default
set(h, 'units', 'pixels')
close(h)
end

% EOF
