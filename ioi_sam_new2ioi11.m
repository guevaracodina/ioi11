function ioi_sam_new2ioi11(IOI)
% Pre-processing stage of OIS images (IOI).
% Reads Hb & HbR images, saves anatomical image, brain Mask
% These images are already band-pass filtered between 1/120 & 1/10Hz
% and spatially filtered 3x3 pixels FWHM
%% dir field
doShrinkage = true;     % True when images are 1024x1024
shrink_factor = 2;
doFlip = true;          % True for July 2016 experiments
fprintf('Processing started for subject %s\n',IOI.subj_name);
IOI.warning = {};
IOI.subj_OK = 1;
IOI.dir.dir_subj_raw = fullfile(IOI.dir.dir_group_raw, IOI.subj_name);
IOI.dir.dir_subj_res = fullfile(IOI.dir.dir_group_res, IOI.subj_name);
if ~exist(IOI.dir.dir_subj_res,'dir')
    mkdir(IOI.dir.dir_subj_res);
end
if ~exist(fullfile(IOI.dir.dir_subj_res,'S01'),'dir')
    mkdir(fullfile(IOI.dir.dir_subj_res,'S01'));
end
IOImat = fullfile(IOI.dir.dir_subj_res,'IOI.mat');

%% Read HbO Data
tic
load(fullfile(IOI.dir.dir_subj_raw,'Dim_binFile.mat'));
fileID=fopen(fullfile(IOI.dir.dir_subj_raw,'HbO.bin'),'r');
HbO = fread(fileID,'int32');
fclose(fileID);
if rem (numel(HbO), X_d2*Y_d3) ~= 0
    Temps_d1 = fix(numel(HbO) ./ (X_d2*Y_d3));
    HbO = HbO(1:Temps_d1*X_d2*Y_d3);
end
HbO = reshape(HbO,Temps_d1,X_d2,Y_d3); %Temps_d1,� are stored in Dim_binFile.mat
disp(['Read HbO Data in: ' datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS')]);

%% Do flip (HbO)
if doFlip
    tic
    for iFrames=1:size(HbO,1)
        HbO(iFrames,:,:) = rot90(squeeze(HbO(iFrames,:,:)), 2);
    end
    disp(['HbO Data Flipped in: ' datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS')]);
end

%% Compute transformation to realign HbO images
tic
% Load registered image
load(fullfile(IOI.dir.dir_subj_raw,[IOI.subj_name '_coregistration.mat']))
transformationType = 'projective';
tform = fitgeotrans(movingPoints,fixedPoints,transformationType);
for iFrames=1:size(HbO,1)
    % Apply transformation & display registered images
    HbO(iFrames,:,:) = imwarp(squeeze(HbO(iFrames,:,:)),tform,'OutputView',imref2d(size(atlas_fixed)));
end
disp(['Realignment of HbO Data in: ' datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS')]);

%% Shrinkage preparation and computation of original brainmask
n_frames = size(HbO,1);
if n_frames > 2000
    n_frames = 2000;
end
IOI.sess.n_frames = n_frames;
IOI.sess_res{1}.n_frames = n_frames;
if doShrinkage
    nx = round(size(HbO,2)/shrink_factor);
    ny = round(size(HbO,3)/shrink_factor);
    HbO_resize = zeros([n_frames nx ny]);
else
    nx = size(HbO,2);
    ny = size(HbO,3);
end

% Get edge surrounding the brain
tic
brainMask = mat2gray(squeeze(mean(HbO,1)));
if ~doShrinkage
    brainMask = ioi_MYimresize(brainMask, 2*[nx ny]);
end
% Invert pixels
brainMask = ~im2bw(brainMask, max(brainMask(:))-eps);
IOI.res.shrinkageOn = 1;
% Only if shrinkage is chosen //EGC
IOI.res.shrink_x = shrink_factor;
IOI.res.shrink_y = shrink_factor;
save(IOImat,'IOI');
disp(['Brain Mask computed in: ' datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS')]);

%% Actual Shrinkage
if doShrinkage
    tic
    for iFrames = 1:n_frames,
        HbO_resize(iFrames,:,:) = ioi_MYimresize(squeeze(HbO(iFrames,:,:)), [nx ny]);
    end
    HbO = HbO_resize;
    clear HbO_resize
    disp(['Shrink HbO Data in: ' datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS')]);
end

%% Display images
% maxVal = max(HbO(:));
% minVal = min(HbO(:));
% figure;
% for iTime=1:size(HbO,1)
%     imagesc(squeeze(HbO(iTime,:,:)),[minVal, maxVal]);
%     axis image
%     drawnow;
%     pause(0.1)
% end

%% Color field
IOI.color.eng = 'RGYLOD';
IOI.color.red = 'R';
IOI.color.green = 'G';
IOI.color.yellow = 'Y';
IOI.color.laser = 'L';
IOI.color.HbO = 'O';
IOI.color.HbR = 'D';

%% dev field
IOI.dev.TR = 0.2;

%% conc field
IOI.conc.baseline_hbt = 100;
IOI.conc.baseline_hbo = 60;
IOI.conc.baseline_hbr = 40;
save(IOImat,'IOI');

%% Anatomical
% Bregma est environ � la ligne 250, colonne 400 et Lambda: 250, 150.
% Seed Radius should be about 0.25mm = 17 pixels
% LPF = 5x5 pixels
suffix_for_anat_file = 'anat'; %to build anatomical image name
sess_label = 'S'; %prefix for name of directories for each session
sess_str = [sess_label gen_num_str(1,2)];
%leave voxel size in arbitrary units for now, for anatomical image
if doShrinkage
    vx_anat = [1 1 1];
else
%     vx_anat = [2 2 1];
    vx_anat = [1 1 1];
end
% h = open(fullfile(IOI.dir.dir_subj_res,[IOI.subj_name '.fig']));
h = open(fullfile(IOI.dir.dir_subj_raw,[IOI.subj_name '.fig']));
set(h,'units','inch')
% Image is double, range: [0 4096]
image_anat = 2^12*mat2gray(getimage);
if ~doShrinkage
    image_anat = ioi_MYimresize(image_anat, 2*[nx ny]);
end
% Do flip (anatomical) --> Not necessary, anatomical image is flipped in
% the realignment step
% if doFlip
%     image_anat = rot90(image_anat, 2);
% end
imwrite(mat2gray(getimage), fullfile(IOI.dir.dir_subj_res,[IOI.subj_name '_anat_S01.png']),...
    'BitDepth',16)
% First, create S01 folder inside IOI.dir.dir_subj_res
% Use first green colored image as the anatomy
anat_fname = fullfile(IOI.dir.dir_subj_res,sess_str, ...
    [IOI.subj_name '_' suffix_for_anat_file '_' sess_str]);
ioi_save_images(image_anat, anat_fname, vx_anat,[],sprintf('%s Anatomical image',IOI.subj_name));
% It will always point to the last anatomical image
IOI.res.file_anat = [anat_fname '.nii'];
% Save IOI matrix
save(IOImat,'IOI');

%% Brain Mask
if ~isfield(IOI, 'fcIOS')
    % Create fcIOS field to contain the whole structure of fcIOS
    % utilities
    IOI.fcIOS = struct([]);
    IOI.fcIOS(1).mask = struct([]);
    IOI.fcIOS(1).LPF  = struct([]);
    IOI.fcIOS(1).filtNdown = struct([]);
    IOI.fcIOS(1).SPM = struct([]);
    IOI.fcIOS(1).corr = struct([]);
end
% Not necessary, brain mask is createed from flipped HbO images
% if doFlip
%     brainMask = rot90(brainMask, 2);
% end
imwrite(brainMask, fullfile(IOI.dir.dir_subj_res,[IOI.subj_name '_brainMask.png']),...
    'BitDepth',1)
[dirName, fileName, fileExt] = fileparts(IOI.res.file_anat);
fileName = strcat(fileName, fileExt);
% Create filename according the existing nomenclature at subject level
brainMaskName = [IOI.subj_name '_anat_brainmask.nii'];

% Create and write a NIFTI file in the subject folder
hdr = spm_vol(fullfile(dirName,fileName));
ioi_create_vol(fullfile(dirName, brainMaskName), ...
    hdr.dim, hdr.dt, hdr.pinfo, hdr.mat, hdr.n, brainMask);
if isempty(IOI.fcIOS.mask)
    % To avoid emptyDotAssignment we create a field
    IOI.fcIOS.mask = struct('fname', []);
end

% Identifier dans IOI le nom du fichier masque
IOI.fcIOS.mask.fname = fullfile(dirName, brainMaskName);
% Mask created succesfully!
IOI.fcIOS.mask.maskOK = true;
% Save IOI matrix
save(IOImat,'IOI');

%% Save concentrations (HbO)
IOI.sess_res{1, 1}.missingFrames = [];
nzero_padding = 5;
image_str_start = gen_num_str(1,nzero_padding);
image_str_last = gen_num_str(IOI.sess.n_frames,nzero_padding);
image_str = [image_str_start 'to' image_str_last];
str_HbO = 'O';
% IOI.sess_res{1} = sess;
% IOI.sess_res{s1}.fname{IOI.color.eng==str_HbO} = fname_new_HbO_list;
% fname_new_HbO = regexprep(fname,['_OD',tmp_str_avail_color] ,tmp_str_HbO);
fname_new_HbO = fullfile(IOI.dir.dir_subj_res,'S01',...
    [IOI.subj_name '_OD_' str_HbO '_S01_' image_str '.nii']);
fname_new_HbO_list = {};
fname_new_HbO_list = [fname_new_HbO_list; fname_new_HbO];
HbO = permute(HbO,[2 3 4 1]);
oNaN = sum(isnan(HbO(:)));
% rNaN = sum(isnan(image_hbr(:)));
oInf = sum(isinf(HbO(:)));
% rInf = sum(isinf(image_hbr(:)));
omax = max(HbO(:));
% rmax = max(image_hbr(:));
if oNaN
    IOI = disp_msg(IOI,[int2str(oNaN) ' NaN in HbO in session ' int2str(s1)]);
    HbO(isnan(HbO(:))) = omax;
end
% if rNaN
%     IOI = disp_msg(IOI,[int2str(rNaN) ' NaN in HbR in session ' int2str(s1)]);
%     image_hbr(isnan(image_hbr(:))) = rmax;
% end
if oInf
    IOI = disp_msg(IOI,[int2str(oInf) ' Inf in HbO in session ' int2str(s1)]);
    HbO(isinf(HbO(:))) = omax;
end
% if rInf
%     IOI = disp_msg(IOI,[int2str(rInf) ' Inf in HbR in session ' int2str(s1)]);
%     image_hbr(isinf(image_hbr(:))) = rmax;
% end
vx_Hb = [2 2 1];
% Normalize between 0 and 1, and convert to single
HbO = single(mat2gray(HbO));
tic
ioi_save_nifti(HbO, fname_new_HbO, vx_Hb);
IOI.sess_res{1}.fname{IOI.color.eng==str_HbO} = fname_new_HbO_list;
IOI.res.concOK = 1;
% Save IOI matrix
save(IOImat,'IOI');
clear HbO
disp(['Saved HbO Data in: ' datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS')]);

%% Read HbR Data
tic
load(fullfile(IOI.dir.dir_subj_raw,'Dim_binFile.mat'));
fileID=fopen(fullfile(IOI.dir.dir_subj_raw,'HbR.bin'),'r');
HbR = fread(fileID,'int32');
fclose(fileID);
if rem (numel(HbR), X_d2*Y_d3) ~= 0
    Temps_d1 = fix(numel(HbR) ./ (X_d2*Y_d3));
    HbR = HbR(1:Temps_d1*X_d2*Y_d3);
end
HbR = reshape(HbR,Temps_d1,X_d2,Y_d3); %Temps_d1,� are stored in Dim_binFile.mat
disp(['Read HbR Data in: ' datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS')]);

%% Do flip (HbR)
if doFlip
    tic
    for iFrames=1:size(HbR,1)
        HbR(iFrames,:,:) = rot90(squeeze(HbR(iFrames,:,:)), 2);
    end
    disp(['HbR Data Flipped in: ' datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS')]);
end

%% Compute transformation to realign HbO images
tic
% Load registered image
load(fullfile(IOI.dir.dir_subj_raw,[IOI.subj_name '_coregistration.mat']))
transformationType = 'projective';
tform = fitgeotrans(movingPoints,fixedPoints,transformationType);
for iFrames=1:size(HbR,1)
    % Apply transformation & display registered images
    HbR(iFrames,:,:) = imwarp(squeeze(HbR(iFrames,:,:)),tform,'OutputView',imref2d(size(atlas_fixed)));
end
disp(['Realignment of HbR Data in: ' datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS')]);

%% Shrinkage preparation
% Only necessary one time
% IOImat = 'D:\Edgar\OIS_Results\FB25E02\IOI.mat';
% load(IOImat)
% n_frames = IOI.sess.n_frames ;
if doShrinkage
    nx = round(size(HbR,2)/IOI.res.shrink_x);
    ny = round(size(HbR,3)/IOI.res.shrink_y);
else
    nx = size(HbR,2);
    ny = size(HbR,3);
end
HbR_resize = zeros([n_frames nx ny]);
save(IOImat,'IOI');

%% Actual Shrinkage
if doShrinkage
    tic
    for iFrames = 1:n_frames,
        HbR_resize(iFrames,:,:) = ioi_MYimresize(squeeze(HbR(iFrames,:,:)), [nx ny]);
    end
    HbR = HbR_resize;
    clear HbR_resize
    disp(['Shrink HbR Data in: ' datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS')]);
end

%% Save concentrations (HbR)
nzero_padding = 5;
image_str_start = gen_num_str(1,nzero_padding);
image_str_last = gen_num_str(IOI.sess.n_frames,nzero_padding);
image_str = [image_str_start 'to' image_str_last];
str_HbR = 'D';
fname_new_HbR = fullfile(IOI.dir.dir_subj_res,'S01',...
    [IOI.subj_name '_OD_' str_HbR '_S01_' image_str '.nii']);
fname_new_HbR_list = {};
fname_new_HbR_list = [fname_new_HbR_list; fname_new_HbR];
HbR = permute(HbR,[2 3 4 1]);
oNaN = sum(isnan(HbR(:)));
oInf = sum(isinf(HbR(:)));
omax = max(HbR(:));
if oNaN
    IOI = disp_msg(IOI,[int2str(oNaN) ' NaN in HbR in session ' int2str(s1)]);
    HbR(isnan(HbR(:))) = omax;
end
if oInf
    IOI = disp_msg(IOI,[int2str(oInf) ' Inf in HbR in session ' int2str(s1)]);
    HbR(isinf(HbR(:))) = omax;
end
vx_Hb = [2 2 1];
% Normalize between 0 and 1, and convert to single
HbR = single(mat2gray(HbR));
tic
ioi_save_nifti(HbR, fname_new_HbR, vx_Hb);
% Check if cell no. 6 has deoxy filename
IOI.sess_res{1}.fname{IOI.color.eng==str_HbR} = fname_new_HbR_list;
IOI.res.concOK = 1;
% Save IOI matrix
save(IOImat,'IOI');
clear HbR
disp(['Saved HbR Data in: ' datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS')]);

%% Missing fields
% Save available colors
IOI.sess_res{1}.availCol = 'RGYLOD';
IOI.res.flowOK = true;
% Save IOI matrix
save(IOImat,'IOI');
fprintf('Processing ended for subject %s\n',IOI.subj_name);
end

% EOF