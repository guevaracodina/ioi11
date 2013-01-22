function out = ioi_fcIOS_maps_run(job)
% Prints multiple correlation maps at the same scale. 
% Ideal to make a mosaique figure.
%_______________________________________________________________________________
% Copyright (C) 2012 LIOM Laboratoire d'Imagerie Optique et Moleculaire
%                    Ecole Polytechnique de Montreal
%_______________________________________________________________________________

% ------------------------------------------------------------------------------
% REMOVE AFTER FINISHING THE FUNCTION //EGC
% ------------------------------------------------------------------------------
% fprintf('Work in progress...\nEGC\n')
% out.IOImat = job.IOImat;
% return
% ------------------------------------------------------------------------------

% Select a subset of sessions
[all_sessions selected_sessions] = ioi_get_sessions(job);
% Colors to include
IC = job.IC;

% Loop over subjects
for SubjIdx = 1:length(job.IOImat)
    try
        tic
        % Load IOI.mat information
        [IOI IOImat dir_ioimat]= ioi_get_IOI(job,SubjIdx);
        % Load brain mask
        brainMaskVol = spm_vol(IOI.fcIOS.mask.fname);
        brainMask = logical(spm_read_vols(brainMaskVol));
        % Read anatomical NIFTI file
        % imAnatVol = spm_vol(IOI.res.file_anat);
        % imAnat = spm_read_vols(imAnatVol);
        
        if IOI.res.shrinkageOn
            vx = [IOI.res.shrink_x IOI.res.shrink_y 1];
        else
            vx = [1 1 1];
        end
        
        % Loop over sessions
        for s1 = 1:length(IOI.sess_res)
            if all_sessions || sum(s1==selected_sessions)
                % Loop over available colors
                for c1=1:length(IOI.sess_res{s1}.fname)
                    doColor = ioi_doColor(IOI,c1,IC);
                    if doColor
                        [all_ROIs selected_ROIs] = ioi_get_ROIs(job);
                        nROI = 1:length(IOI.res.ROI); % All the ROIs
                        % Loop over ROI/seeds
                        for r1 = nROI,
                            if (all_ROIs || sum(r1==selected_ROIs)) && IOI.fcIOS.corr.corrMapOK{r1}{s1, c1}
                                %% Load correlation data
                                tempCorrMapVol = spm_vol(IOI.fcIOS.corr.corrMapName{r1}{s1, c1});
                                tempCorrMap = spm_read_vols(tempCorrMapVol);
                                [~, ~, oldExt] = fileparts(IOI.fcIOS.SPM.fnameROInifti{r1}{s1, c1});
                                newName = [sprintf('%s_R%02d_S%02d_C%d',IOI.subj_name,r1,s1,c1) '_fcIOS_corrMapMasked'];
                                dir_corrMapMaskedfig = fullfile(dir_ioimat,'fig_corrMapMasked');
                                if ~exist(dir_corrMapMaskedfig,'dir'), mkdir(dir_corrMapMaskedfig); end
                                
                                %% Mask out non-brain voxels
                                if size(brainMask,1)~= size(tempCorrMap,1)|| size(brainMask,2)~= size(tempCorrMap,2)
                                    brainMask = ioi_MYimresize(brainMask, [size(tempCorrMap,1) size(tempCorrMap,2)]);
                                end
                                % Make NaN all the transparent pixels
                                tempCorrMap(~brainMask) = NaN;
                                y = tempCorrMap;
                                IOI.fcIOS.corr.corrMapNameMask{r1}{s1, c1} = fullfile(dir_corrMapMaskedfig,[newName oldExt]);
                                % Save NIfTI to use slover
                                ioi_save_nifti(y, IOI.fcIOS.corr.corrMapNameMask{r1}{s1, c1}, vx);
                                
                                %% Display overlay data and print to png
                                % Get parameters for overlay
                                anatomical      = IOI.res.file_anat;
                                corrMap         = IOI.fcIOS.corr.corrMapNameMask{r1}{s1, c1};
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
                                internal_overlay_map(anatomical, corrMap,  job, newName, [seedX seedY seedW seedH]);
                            end
                        end % ROIs loop
                    end
                end % colors loop
            end
        end % sessions loop
        out.IOImat{SubjIdx} = IOImat;
        disp(['Elapsed time: ' datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS')]);
        disp(['Subject ' int2str(SubjIdx) ' (' IOI.subj_name ')' ' complete']);
    catch exception
        disp(exception.identifier)
        disp(exception.stack(1))
        out.IOImat{SubjIdx} = job.IOImat{SubjIdx};
    end
end % loop over subjects
end % Main function

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
% job           Matlab batch jobjob 
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
% rotx = @(theta) [1 0 0 0; 0 cos(theta) -sin(theta) 0; 0 sin(theta) cos(theta) 0; 0 0 0 1];
% roty = @(theta) [cos(theta) 0 sin(theta) 0; 0 1 0 0; -sin(theta) 0 cos(theta) 0; 0 0 0 1];
% rotz = @(theta) [cos(theta) -sin(theta) 0 0; sin(theta) cos(theta) 0 0; 0 0 1 0; 0 0 0 1];
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
