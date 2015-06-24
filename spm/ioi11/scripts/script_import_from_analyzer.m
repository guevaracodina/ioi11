%% script_import_from_analyzer
%% Loading data after filtering and downsampling
% Run spm8 shortcut first
clc
addpath(genpath('C:\spm8\toolbox\nirs10')); % After running spm8 shortcut
pathNameParent = 'C:\Edgar\Data\IOIResults\j03\ROItest';
IOImat = fullfile(pathNameParent, 'IOI.mat');
load(IOImat)
if IOI.fcIOS.filtNdown.filtNdownOK
	% Load filtered & downsampled image time-course
    vol = spm_vol(IOI.fcIOS.filtNdown.fnameWholeImage{1});
    y = spm_read_vols(vol);
end
% Get image dimensions
[nX, nY, nZ, nT] = size(y);
% Reshape as 2-D matrix
y = reshape(y,[nX*nY nT]);
% Intrinsic Imaging frequency (after downsampling to 1 Hz)
Fs = IOI.fcIOS.filtNdown.fs;
Ts = 1/Fs;
load(fullfile(pathNameParent, 'j03_resp_noCB_Raw_Data.mat'));
%% Downsample files
pathNamePhysio = 'C:\Edgar\Data\IOIData20141211\Physio Monitoring\121423';

%% Compute resampling ratio
% Free up memory
clear b physio
pack
% resampling at p/q times the original sampling rate
q = round(resp_freq);
p = round(Fs);
% Cut data to experiment length
nT = size(y,2)*resp_freq;
respFilt = respFilt(1:nT);
t = t(1:nT);

%% Memory-mapping
tic
clear im_obj
fmem_name=fullfile(pathNameParent, [IOI.subj_name '_resampled.dat']);
fid = fopen(fmem_name,'w');
nPix = size(y,1)+1;             % Number of pixel time-traces
fwrite(fid, zeros([nPix nT]), 'double');
fclose(fid);
im_obj = memmapfile(fmem_name,...
        'Format',{'double' [nPix nT] 'Frames'}, ...
        'writable',true);
toc    

%% Interpolate & Resample
% [ECGresamp,b] = resample(ECGfilt,p,q);
tic
ioi_text_waitbar(0, 'Please wait...');
% Concatenate data
im_obj.Data.Frames(1,:) = respFilt;
% [yResamp,b] = resample(y,q,p);

for iPix = 2:nPix
    im_obj.Data.Frames(iPix,:) = resample(y(iPix-1,:),q,p);
    % Update progress bar
      ioi_text_waitbar(iPix/nPix, sprintf('Processing event %d from %d', iPix, nPix));
end
ioi_text_waitbar('Clear');
toc

%% Write data to be imported in Analyzer2
fwrite_NIR(fullfile(pathNameParent,[IOI.subj_name '_resp_noCB.nir']), im_obj.Data.Frames);

% EOF