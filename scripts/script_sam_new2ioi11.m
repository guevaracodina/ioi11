%% Path (after running spm8)
addpath(genpath('D:\spm8\toolbox\ioi'))

%% dir field
IOI.warning = {};
IOI.subj_name = 'FB25E02';
IOI.subj_OK = 1;
IOI.dir.dir_group_all = 'D:\Edgar\';
IOI.dir.dir_group_raw = 'D:\Edgar\OIS_Data\';
IOI.dir.dir_group_res = 'D:\Edgar\OIS_Results\';
IOI.dir.dir_subj_raw = fullfile(IOI.dir.dir_group_raw, IOI.subj_name);
IOI.dir.dir_subj_res = fullfile(IOI.dir.dir_group_res, IOI.subj_name);

%% Read HbO Data
load(fullfile(IOI.dir.dir_subj_raw,'Dim_binFile.mat'));
fileID=fopen(fullfile(IOI.dir.dir_subj_raw,'HbO.bin'),'r');
HbO = fread(fileID,'int32');
HbO = reshape(HbO,Temps_d1,X_d2,Y_d3); %Temps_d1,… are stored in Dim_binFile.mat

%% Shrinkage preparation and original brainmask
shrink_factor = 2;
n_frames = size(HbO,1);
IOI.sess.n_frames = n_frames;
nx = round(size(HbO,2)/shrink_factor);
ny = round(size(HbO,3)/shrink_factor);
HbO_resize = zeros([size(HbO,1) nx ny]);
brainMask = mat2gray(squeeze(mean(HbO,1)));
brainMask = im2bw(brainMask, max(brainMask(:))-eps);
IOI.res.shrinkageOn = 1;
% Only if shrinkage is chosen //EGC
IOI.res.shrink_x = shrink_factor;
IOI.res.shrink_y = shrink_factor;
save(IOI.dir.dir_subj_res,'IOI');

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
save(IOI.dir.dir_subj_res,'IOI');
%% Anatomical
% Bregma est environ à la ligne 250, colonne 400 et Lambda: 250, 150.
suffix_for_anat_file = 'anat'; %to build anatomical image name
sess_label = 'S'; %prefix for name of directories for each session
sess_str = [sess_label gen_num_str(1,2)];
%leave voxel size in arbitrary units for now, for anatomical image
vx_anat = [1 1 1];
% image_anat = mat2gray(squeeze(mean(HbO,1)));
h = open(fullfile(IOI.dir.dir_subj_res,[IOI.subj_name '.fig']));
set(h,'units','inch')
% Image is double, range: [0 4096]
% anat_fname =
% D:\Edgar\OIS_Results\12_10_18,NC09\S01\12_10_18,NC09_anat_S01
% vx_anat =
% 
%      1     1     1
% 
% 12_10_18,NC09 Anatomical image

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
save(IOI.dir.dir_subj_res,'IOI');

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
% brainMask = (ioi_MYimresize(squeeze(mean(HbO,1)), shrink_factor*[nx ny]));
% brainMask=im2bw(brainMask, max(brainMask(:))-eps);
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
save(IOI.dir.dir_subj_res,'IOI');

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
% size(image_hbo) = 303   321     1   128
% fname_new_HbO =
% D:\Edgar\OIS_Results\12_10_18,NC09\S01\12_10_18,NC09_O_S01_00001to00128.nii
vx_Hb = [2 2 1];
tic
ioi_save_nifti(HbO, fname_new_HbO, vx_Hb);
toc
IOI.sess_res{1}.fname{IOI.color.eng==str_HbO} = fname_new_HbO_list;
IOI.res.concOK = 1;
% Save IOI matrix
save(IOI.dir.dir_subj_res,'IOI');