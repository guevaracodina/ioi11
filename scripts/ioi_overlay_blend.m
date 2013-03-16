function ioi_overlay_blend(IOImat, job, fcMapFile)
%% Overlay blend
fcMapRange = [-1 1]; % Fix range for correlation maps
nColorLevels = 256;
load(IOImat);
if ~exist(job.parent_results_dir{1},'dir'),
    mkdir(job.parent_results_dir{1})
end
[~, fName] = fileparts(fcMapFile);
r1 = str2double(regexp(fName, '(?<=(_R))(\d+)(?=(C))','match'));
currentName = regexp(fName, '.*(?=(_avg))', 'match');
figName = fullfile(job.parent_results_dir{1} ,[currentName{1} '_avg_overlay']);

%%
% Seed annotation dimensions the lower left corner of the bounding rectangle
% at the point seedX, seedY (the + sign here is due to image rotation)
seedX = IOI.res.ROI{r1}.center(2) - IOI.res.ROI{r1}.radius;
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

%% Read files

% Get anatomical image
anatomicalFile      = IOI.res.file_anat;
anatVol             = spm_vol(anatomicalFile);
anatomical          = spm_read_vols(anatVol);

% Read brain mask
maskVol             = spm_vol(IOI.fcIOS.mask.fname);
brainMask           = logical(spm_read_vols(maskVol));
if size(brainMask,1)~= size(anatomical,1)|| size(brainMask,2)~= size(anatomical,2)
    brainMask       = ioi_MYimresize(brainMask, [size(anatomical,1) size(anatomical,2)]);
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

%% Convert to grayscale
anatomicalGray      = mat2gray(anatomical);
anatomicalGray      = repmat(anatomicalGray,[1 1 3]);
% Convert to RGB
fcMapGray           = mat2gray(fcMap, fcMapRange); % Fix range for correlation maps
fcMapX              = gray2ind(fcMapGray, nColorLevels);
fcMapRGB            = ind2rgb(fcMapX, jet(nColorLevels));
% Set transparency according to mask
fcMapRGB(repmat(~brainMask,[1 1 3])) = 0.5;

%% Apply overlay blend
fcMapBlend = 1 - 2.*(1 - anatomicalGray).*(1 - fcMapRGB);
fcMapBlend(anatomicalGray<0.5) = 2.*anatomicalGray(anatomicalGray<0.5).*fcMapRGB(anatomicalGray<0.5);

%% Print
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
% Save as PNG at the user-defined resolution
print(hFig, '-dpng', ...
    figName,...
    sprintf('-r%d',job.figRes));
% Return the property to its default
set(hFig, 'units', 'pixels')
close(hFig)
colorNames = fieldnames(IOI.color);
c1 = str2double(regexp(fcMapFile, '(?<=(C))(\d+)(?=(_))','match'));
fprintf('Overlay blend done! File: %s, R%02d, (%s)\n', fName, r1, colorNames{c1+1});
end % script_overlay_blend
% EOF
