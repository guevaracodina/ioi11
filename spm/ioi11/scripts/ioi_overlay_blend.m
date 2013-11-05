function [displayed_pixels, total_pixels, h] = ioi_overlay_blend(IOImat, job, fcMapFile, varargin)
% Overlay/blend a functional image from a NIfTI file onto an anatomical image.
% SYNTAX
% ioi_overlay_blend(IOImat, job, fcMapFile, fcMapRange, alphaRange, nColorLevels)
% INPUTS
% IOImat        .mat file with IOI structure
% job           Matlab batch job structure with the required fields.
%               For instance:
%               job(1).figCmap                                  = jet(256);     % colormap
%               job(1).figIntensity                             = 1;            % [0 - 1]
%               job(1).transM                                   = rotz(pi);     % affine transform
%               job(1).figSize                                  = [1.5 1.5];    % inches
%               job(1).figRes                                   = 300;          % in dpi
%               job(1).drawCircle(1).drawCircle_On(1).circleLW  = 0.8;          % line width
%               job(1).drawCircle(1).drawCircle_On(1).circleLS  = '-';          % line style
%               job(1).drawCircle(1).drawCircle_On(1).circleEC  = 'w';          % line color
%               job.parent_results_dir{1}                       = fullfile(figFolder,'overlay');
%               job.generate_figures                            = true;         % display figure
%               job.save_figures                                = true;         % save figure
% fcMapFile     NIfTI file name of the functional map (foreground). If you add a
%               coma, you can choose specific frame from 4-D NIfTI data.
% [fcMapRange]  Range of values to map to the full range of colormap
%               if empty, functional map is automatically displayed from [minVal maxVal]
% [alphaRange]  Range of values to map to display non-transparent pixels.
%               Useful to set a threshold of displayed pixels.
%               if empty, alpha values are automatically taken from [minVal maxVal]
%               Can specify up to two ranges to display [minVal1 maxVal1 minVal2 maxVal2]
%               Useful for symmetric data ranges, such as in correlation maps.
% [nColorLevels]Number of gray levels for the anatomical image (background)
%               Default: 256 gray levels.
% [r1]          ROI/seed index
% [c1]          Color/contrast index
% OUTPUTS
% none
%_______________________________________________________________________________
% Copyright (C) 2012 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________

%% -----------------------------------------------------------------------------
% Optional inputs handling
% ------------------------------------------------------------------------------
% only want 1 optional input at most
numvarargs                  = length(varargin);
if numvarargs > 5
    error('ioi_overlay_blend:TooManyInputs', ...
        'Requires at most 5 optional inputs');
end
% set defaults for optional inputs
[~, fName]                  = fileparts(fcMapFile);
optargs                     = { [] [] 256 ...
    str2double(regexp(fName, '(?<=(_R))(\d+)(?=(C))','match'))...
    str2double(regexp(fcMapFile, '(?<=(C))(\d+)(?=(_))','match')) + 1};
% now put these defaults into the optargs cell array, and overwrite the ones
% specified in varargin.
optargs(1:numvarargs)       = varargin;
% Place optional args in memorable variable names
[fcMapRange alphaRange ...
    nColorLevels r1 c1]= optargs{:};
% ------------------------------------------------------------------------------

%% Overlay blend
fcColorMap      = job.figCmap;
load(IOImat);
if ~exist(job.parent_results_dir{1},'dir'),
    mkdir(job.parent_results_dir{1})
end

% if isempty(r1)
% 	r1 = str2double(regexp(fName, '(?<=(_R))(\d+)(?=(C))','match'));
% end
currentName = regexp(fName, '.*(?=(_avg))', 'match');
if isempty(currentName)
    currentName = {fName};
end
figName = fullfile(job.parent_results_dir{1} ,[currentName{1} '_avg_overlay']);

%% Seed positions
% Seed annotation dimensions the lower left corner of the bounding rectangle
% at the point seedX, seedY (the + sign here is due to image rotation)
seedX = IOI.res.ROI{r1}.center(2) - IOI.res.ROI{r1}.radius;
seedY = IOI.res.ROI{r1}.center(1) + IOI.res.ROI{r1}.radius;
if r1 == 12 && ~isempty(regexp(IOImat, 'NC', 'once')), % Manual correction of seed 12 (NC)
    seedX = seedX - 15;
    seedY = seedY + 10;
end
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

%% Read files
% Get anatomical image
anatomicalFile      = IOI.res.file_anat;
anatVol             = spm_vol(anatomicalFile);
anatomical          = spm_read_vols(anatVol);

if exist(IOI.fcIOS.mask.fname, 'file')
    % Read brain mask
    maskVol             = spm_vol(IOI.fcIOS.mask.fname);
    brainMask           = logical(spm_read_vols(maskVol));
    if size(brainMask,1)~= size(anatomical,1)|| size(brainMask,2)~= size(anatomical,2)
        brainMask       = ioi_MYimresize(brainMask, [size(anatomical,1) size(anatomical,2)]);
    end
else
    brainMask = ones(size(anatomical));
end

% Get functional image
fcMapVol            = spm_vol(fcMapFile);
fcMap               = spm_read_vols(fcMapVol);

if size(fcMap,1)~= size(anatomical,1)|| size(fcMap,2)~= size(anatomical,2)
    fcMap           = ioi_MYimresize(fcMap, [size(anatomical,1) size(anatomical,2)]);
end

% Orient images
anatomical          = fliplr(anatomical');
fcMap               = fliplr(fcMap');
brainMask           = fliplr(brainMask');
seedDims = [size(fcMap,2) - seedY seedX seedW seedH];

%% If values are empty display min/max
if isempty(fcMapRange)
    fcMapRange = [min(fcMap(:)) max(fcMap(:))];
    fprintf('fcMapRange = [%0.4f %0.4f]\n', fcMapRange(1), fcMapRange(2));
end
if isempty(alphaRange) || (numel(alphaRange) ~= 2 && numel(alphaRange) ~= 4)
    alphaRange = [min(fcMap(:)) max(fcMap(:))];
    fprintf('alphaRange = [%0.4f %0.4f]\n', alphaRange(1), alphaRange(2));
end

%% Convert anatomical image to grayscale (weighted by job.figIntensity)
anatomicalGray      = job.figIntensity .* mat2gray(anatomical);
anatomicalGray      = repmat(anatomicalGray,[1 1 3]);
% Convert functional image to RGB
fcMapGray           = mat2gray(fcMap, fcMapRange); % Fix range for correlation maps
fcMapIdx            = gray2ind(fcMapGray, nColorLevels);
fcMapRGB            = ind2rgb(fcMapIdx, fcColorMap);
% Set transparency according to mask and pixels range
pixelMask = false(size(brainMask));
if numel(alphaRange) == 2,
    pixelMask(fcMap > alphaRange(1) & fcMap < alphaRange(2)) = true;
elseif numel(alphaRange) == 4,
    pixelMask(fcMap > alphaRange(1) & fcMap < alphaRange(2)) = true;
    pixelMask(fcMap > alphaRange(3) & fcMap < alphaRange(4)) = true;
end
fcMapRGB(repmat(~brainMask | ~pixelMask,[1 1 3])) = 0.5;
% Spatial extension % defined as displayed to brain pixels ratio.
spatial_extension = nnz(pixelMask) / nnz(brainMask);
displayed_pixels = fcMap(pixelMask);
total_pixels = nnz(brainMask);

%% Apply overlay blend algorithm
fcMapBlend = 1 - 2.*(1 - anatomicalGray).*(1 - fcMapRGB);
fcMapBlend(anatomicalGray<0.5) = 2.*anatomicalGray(anatomicalGray<0.5).*fcMapRGB(anatomicalGray<0.5);

%% Generate/Print figures
if job.generate_figures
    % h = figure;
    h = gcf;
    h = imshow(fcMapBlend, 'InitialMagnification', 'fit', 'border', 'tight');
    hFig = gcf;
    set(hFig, 'color', 'k')
    % Allow printing of black background
    set(hFig, 'InvertHardcopy', 'off');
    % Specify window units
    set(hFig, 'units', 'inches')
    % Change figure and paper size
    set(hFig, 'Position', [0.1 0.1 job.figSize(1) job.figSize(2)])
    set(hFig, 'PaperPosition', [0.1 0.1 job.figSize(1) job.figSize(2)])
    
    % Seed positions and sizes will be shown with black circles
    if isfield(job.drawCircle,'drawCircle_On')
        figure(hFig);
        % Display ROI
        rectangle('Position',seedDims,...
            'Curvature',[1,1],...
            'LineWidth',job.drawCircle.drawCircle_On.circleLW,...
            'LineStyle',job.drawCircle.drawCircle_On.circleLS,...
            'EdgeColor',job.drawCircle.drawCircle_On.circleEC);
    end
    if job.save_figures
        % Save as PNG at the user-defined resolution
        print(hFig, '-dpng', ...
            figName,...
            sprintf('-r%d',job.figRes));
        % Return the property to its default
        set(hFig, 'units', 'pixels')
        close(hFig)
    end
    colorNames = fieldnames(IOI.color);
%     if isempty(c1)
%         c1 = str2double(regexp(fcMapFile, '(?<=(C))(\d+)(?=(_))','match'));
%     else
%         c1 = c1-1;
%     end
    fprintf('Overlay blend done! File: %s, R%02d, (%s) %0.2f%% brain pixels displayed.\n',...
        fName, r1, colorNames{c1+1},100*spatial_extension);
end
end % ioi_overlay_blend
% EOF
