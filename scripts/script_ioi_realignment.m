%% script_ioi_realignment
% Run first this script, then script_sam_new2ioi11
%% Path (after running spm8)
addpath(genpath('D:\spm8\toolbox\ioi'))
doShrinkage = false;
shrink_factor = 1;
doFlip = true;

%% dir field
IOI.dir.dir_group_all = 'D:\Edgar\';
IOI.dir.dir_group_raw = 'D:\Edgar\OIS_Data\';
IOI.dir.dir_group_res = 'D:\Edgar\OIS_Results\';
atlasDir = fullfile(IOI.dir.dir_group_raw, '16_02_25,NC01');
h =  open(fullfile(atlasDir,['16_02_25,NC01' '.fig']));
atlas_fixed = mat2gray(getimage);
% Load moving points as a starting point
load('D:\Edgar\OIS_Data\16_02_25,NC02\16_02_25,NC02_coregistration.mat')
if ~doShrinkage
    movingPoints = movingPoints/2;
end

%% Subjects missing anatomical figure
% subjectList{1} = '16_02_25,NC03b';

%% Create anatomical figure
% for iSubjects=1:numel(subjectList)
%     tic
%     IOI.subj_name = subjectList{iSubjects};
%     IOI.dir.dir_subj_raw = fullfile(IOI.dir.dir_group_raw, IOI.subj_name);
%     load (fullfile(IOI.dir.dir_subj_raw,'Data_green.mat'));
%     im_anat = squeeze(Frames(:,:,1));
%     h = figure;
%     imagesc(im_anat); axis image; colormap('gray')
%     % Save figure
%     saveas(h, fullfile(IOI.dir.dir_subj_raw,[IOI.subj_name '.fig']), 'fig');
%     disp(['Save anatomical image in: ' datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS')]);
% end

%% Subject list for realignment
clear subjectList
subjectList{1} = '16_02_25,LP01a';
subjectList{2} = '16_02_25,LP01b';
subjectList{3} = '16_02_25,NC02';
subjectList{4} = '16_02_25,NC03a';
subjectList{5} = '16_02_25,NC03b';
subjectList{6} = '16_02_26,NC05a';
subjectList{7} = '16_02_26,NC05b';
subjectList{8} = '16_02_26,NC06a';
subjectList{9} = '16_02_26,NC06b';
subjectList{10} = '16_07_07,NC07';

%% Create alignment points
subjects2Run = 10;      % List of subject numbers to run
for iSubjects = subjects2Run
    IOI.subj_name = subjectList{iSubjects};
    IOI.dir.dir_subj_raw = fullfile(IOI.dir.dir_group_raw, IOI.subj_name);
    h = open(fullfile(IOI.dir.dir_subj_raw,[IOI.subj_name '.fig']));
    set(h,'units','inch')
    % Choose moving image
    anat_moving = mat2gray(getimage);
    if doFlip
        anat_moving = rot90(anat_moving, 2);
    end
    close(h)
    % Export control points (n x 2 matrix) fixedPoints movingPoints
    [movingPoints, fixedPoints] = cpselect(anat_moving, atlas_fixed, movingPoints, fixedPoints, 'Wait', true);
    % After exporting, compute transformation
    transformationType = 'projective';
    tform = fitgeotrans(movingPoints,fixedPoints,transformationType);
    % Apply transformation & display registered images
    anat_registered = imwarp(anat_moving,tform,'OutputView',imref2d(size(atlas_fixed)));
    falsecolorOverlay = imfuse(atlas_fixed, anat_registered);
    figure
    imshow(falsecolorOverlay,'InitialMagnification','fit');
    % Save registered image
    save(fullfile(IOI.dir.dir_subj_raw,[IOI.subj_name '_coregistration.mat']),'fixedPoints', 'movingPoints', 'anat_registered','anat_moving','atlas_fixed')
    h = figure;
    imagesc(anat_registered); axis image; colormap('gray')
    % Save registered anatomical image
    saveas(h, fullfile(IOI.dir.dir_subj_raw,[IOI.subj_name '.fig']), 'fig')
end

% EOF


