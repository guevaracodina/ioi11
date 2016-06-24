function ioi_sam_new2ioi11(IOI)
%% Path (after running spm8)
addpath(genpath('D:\spm8\toolbox\ioi'))

%% dir field
IOI.warning = {};
IOI.subj_OK = 1;
% IOI.subj_name = '16_02_25,NC01';
% IOI.dir.dir_group_all = 'D:\Edgar\';
% IOI.dir.dir_group_raw = 'D:\Edgar\OIS_Data\';
% IOI.dir.dir_group_res = 'D:\Edgar\OIS_Results\';
IOI.dir.dir_subj_raw = fullfile(IOI.dir.dir_group_raw, IOI.subj_name);
IOI.dir.dir_subj_res = fullfile(IOI.dir.dir_group_res, IOI.subj_name);
IOImat = fullfile(IOI.dir.dir_subj_res,'IOI.mat');

%% Read HbO Data
tic
load(fullfile(IOI.dir.dir_subj_raw,'Dim_binFile.mat'));
fileID=fopen(fullfile(IOI.dir.dir_subj_raw,'HbO.bin'),'r');
HbO = fread(fileID,'int32');
HbO = reshape(HbO,Temps_d1,X_d2,Y_d3); %Temps_d1,… are stored in Dim_binFile.mat
toc

%% Shrinkage preparation and original brainmask
shrink_factor = 2;
n_frames = size(HbO,1);
if n_frames > 2000
    n_frames = 2000;
end
IOI.sess.n_frames = n_frames;
IOI.sess_res{1}.n_frames = n_frames;
nx = round(size(HbO,2)/shrink_factor);
ny = round(size(HbO,3)/shrink_factor);
HbO_resize = zeros([n_frames nx ny]);
% Get edge surrounding the brain
brainMask = mat2gray(squeeze(mean(HbO,1)));
% Invert pixels
brainMask = ~im2bw(brainMask, max(brainMask(:))-eps);
IOI.res.shrinkageOn = 1;
% Only if shrinkage is chosen //EGC
IOI.res.shrink_x = shrink_factor;
IOI.res.shrink_y = shrink_factor;
save(IOImat,'IOI');

%% Actual Shrinkage
tic
for iFrames = 1:n_frames,
    HbO_resize(iFrames,:,:) = ioi_MYimresize(squeeze(HbO(iFrames,:,:)), [nx ny]);
end
toc
HbO = HbO_resize;
clear HbO_resize

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
% Bregma est environ à la ligne 250, colonne 400 et Lambda: 250, 150.
% Seed Radius should be about 0.25mm = 17 pixels
% LPF = 5x5 pixels
suffix_for_anat_file = 'anat'; %to build anatomical image name
sess_label = 'S'; %prefix for name of directories for each session
sess_str = [sess_label gen_num_str(1,2)];
%leave voxel size in arbitrary units for now, for anatomical image
vx_anat = [1 1 1];
h = open(fullfile(IOI.dir.dir_subj_res,[IOI.subj_name '.fig']));
set(h,'units','inch')
% Image is double, range: [0 4096]
image_anat = 2^12*mat2gray(getimage);
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
imwrite(brainMask, fullfile(IOI.dir.dir_subj_res,[IOI.subj_name '_brainMask.png']),...
    'BitDepth',1)
[dirName fileName fileExt] = fileparts(IOI.res.file_anat);
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
toc
IOI.sess_res{1}.fname{IOI.color.eng==str_HbO} = fname_new_HbO_list;
IOI.res.concOK = 1;
% Save IOI matrix
save(IOImat,'IOI');

%% Read HbR Data
tic
load(fullfile(IOI.dir.dir_subj_raw,'Dim_binFile.mat'));
fileID=fopen(fullfile(IOI.dir.dir_subj_raw,'HbR.bin'),'r');
HbR = fread(fileID,'int32');
HbR = reshape(HbR,Temps_d1,X_d2,Y_d3); %Temps_d1,… are stored in Dim_binFile.mat
toc

%% Shrinkage preparation
% Only necessary one time
% IOImat = 'D:\Edgar\OIS_Results\FB25E02\IOI.mat';
% load(IOImat)
% n_frames = IOI.sess.n_frames ;
nx = round(size(HbR,2)/IOI.res.shrink_x);
ny = round(size(HbR,3)/IOI.res.shrink_y);
HbR_resize = zeros([n_frames nx ny]);
save(IOImat,'IOI');

%% Actual Shrinkage
tic
for iFrames = 1:n_frames,
    HbR_resize(iFrames,:,:) = ioi_MYimresize(squeeze(HbR(iFrames,:,:)), [nx ny]);
end
toc
HbR = HbR_resize;
clear HbR_resize

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
toc
% Check if cell no. 6 has deoxy filename
IOI.sess_res{1}.fname{IOI.color.eng==str_HbR} = fname_new_HbR_list;
IOI.res.concOK = 1;
% Save IOI matrix
save(IOImat,'IOI');

%% Missing fields
% Save available colors
IOI.sess_res{1}.availCol = 'RGYLOD';
IOI.res.flowOK = true;
% Save IOI matrix
save(IOImat,'IOI');

end

% EOF