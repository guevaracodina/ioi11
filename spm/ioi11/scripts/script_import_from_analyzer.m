%% script_import_from_analyzer
%% Loading data after filtering and downsampling
% Run spm8 shortcut first
clc
addpath(genpath('C:\spm8\toolbox\nirs10')); % After running spm8 shortcut
pathNameParent = 'C:\Edgar\Data\IOIResults\j03\ROItest';
IOImat = fullfile(pathNameParent, 'IOI.mat');
load(IOImat)
clear y
if IOI.fcIOS.filtNdown.filtNdownOK
	% Load filtered & downsampled image time-course
    % vol = spm_vol(IOI.fcIOS.filtNdown.fnameWholeImage{1});
    vols = spm_vol(IOI.sess_res{1, 1}.fname{1});
    for c1 = 1:length(vols) %for each color
        y(c1,:,:,:) = squeeze(spm_read_vols(vols{c1}));
    end
%     %normalization_choice
%     if normalization_choice
%         y = ioi_normalization_choice(y,med0,med1);
%     end
%     y=reshape(y,[length(vols),nx*ny*nt]);
%     y = spm_read_vols(vols);
end
% Get image dimensions
[nZ, nX, nY, nT] = size(y);
% Reshape as 2-D matrix
y = reshape(y,[nX*nY*nZ nT]);
% Intrinsic Imaging frequency (after downsampling to 1 Hz)
% Fs = IOI.fcIOS.filtNdown.fs;
Fs = IOI.dev.TR;
Ts = 1/Fs;
tic
load(fullfile(pathNameParent, 'j03_resp_noCB_Raw_Data.mat'));
toc

%% Downsample files
pathNamePhysio = 'C:\Edgar\Data\IOIData20141211\Physio Monitoring\121423';

%% Compute resampling ratio
% Free up memory and repack
clear b physio
pack
% resampling at p/q times the original sampling rate
% q = round(Fs);
q = round(Fs);
p = round(Fs*SampleRate);
% Cut data to experiment length
nT = size(y,2)*IOI.fcIOS.filtNdown.fs;
nPix = size(y,1)+1;             % Number of pixel time-traces
% respFilt = respFilt(1:nT);
t = t(1:nT);


%% Memory-mapping
% tic
% clear im_obj
% fmem_name=fullfile(pathNameParent, [IOI.subj_name '_resampled.dat']);
% fid = fopen(fmem_name,'w');
% nPix = size(y,1)+1;             % Number of pixel time-traces
% fwrite(fid, zeros([nPix nT]), 'double');
% fclose(fid);
% im_obj = memmapfile(fmem_name,...
%         'Format',{'double' [nPix nT] 'Frames'}, ...
%         'writable',true);
% toc    


%% Interpolate & Resample
% [ECGresamp,b] = resample(ECGfilt,p,q);
tic
ioi_text_waitbar(0, 'Please wait...');
yResamp = zeros(size(y));
for iPix = 2:nPix
    tmpChannel = eval(sprintf('C%d', iPix));
    tmpChannel = resample(tmpChannel,q,p);
    yResamp(iPix-1, :) = tmpChannel';
    % Update progress bar
    ioi_text_waitbar(iPix/nPix, sprintf('Processing pixel timecourse %d from %d', iPix, nPix));
end
yResamp = reshape (yResamp, [nX, nY, nZ, nT]);
% clear channels that begin with C and are followed by a capital C
clearvars -regexp ^C\d
ioi_text_waitbar('Clear');
toc

%% Backup original data and rewrite data
% currentFile = IOI.fcIOS.filtNdown.fnameWholeImage{1};
currentFile = IOI.sess_res{1, 1}.fname{1}{1};
[pathstr,name,ext] = fileparts(currentFile);
newFile = fullfile(pathstr, [name '.bak']);
% Rename .nii to .bak
try
    % Fast alternative (undocumented java feature)
    java.io.File(currentFile).renameTo(java.io.File(newFile));
catch exception
    % Slow alternative...
    movefile(currentFile, newFile);
    disp(exception.identifier)
    disp(exception.stack(1))
end
% Write resampled and corrected data in NIfTI format
ioi_save_nifti(yResamp, currentFile, [1 1 Ts]);

%% Write data to be imported in Analyzer2
% fwrite_NIR(fullfile(pathNameParent,[IOI.subj_name '_resp_noCB.nir']), im_obj.Data.Frames);

% EOF