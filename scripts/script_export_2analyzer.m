%% script_export_2analyzer
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

%% Load EKG & Respiration Data
%  pathNamePhysio = 'C:\Edgar\Data\IOIData20141016\Physio_Monitoring\150642_J03';
% pathNamePhysio = 'C:\Edgar\Data\IOIData20141023\Physio Monitoring\111141_K04';
pathNamePhysio = 'C:\Edgar\Data\IOIData20141211\Physio Monitoring\121423';
% cd (pathNamePhysio)
physio = ioi_Bin2Mat(pathNamePhysio);
% fid=fopen(fullfile(pathNamePhysio,'ECG1.bin'),'r');
% data = uint8(fread(fid));
% ECG1 = double(swapbytes(typecast(data,'int32')))/8388608*2400;
% fclose(fid);
% fid=fopen(fullfile(pathNamePhysio,'ECG2.bin'),'r');
% data = uint8(fread(fid));
% ECG2 = double(swapbytes(typecast(data,'int32')))/8388608*2400;
% fclose(fid);
% fid=fopen(fullfile(pathNamePhysio,'ECG3.bin'),'r');
% data = uint8(fread(fid));
% ECG3 = double(swapbytes(typecast(data,'int32')))/8388608*2400;
% fclose(fid);
% fid=fopen(fullfile(pathNamePhysio,'ECG4.bin'),'r');
% data = uint8(fread(fid));
% ECG4 = double(swapbytes(typecast(data,'int32')))/8388608*2400;
% fclose(fid);
% 
% %Read resp data
% fid=fopen(fullfile(pathNamePhysio,'resp.bin'),'r');
% data = uint8(fread(fid));
% Resp = double(swapbytes(typecast(data,'int16')));
% fclose(fid);

ECG_freq = 1000;        % sampling rate (EKG)
ECG_Ts = 1/ECG_freq;    % sampling time (EKG) interval
resp_freq = 250;        % sampling rate (respiration)
resp_Ts = 1/resp_freq;  % sampling time (respiration) interval
t = 0:resp_Ts:(numel(physio.Resp)-1).*resp_Ts;

%% ECG filtering
% % 10-200 Hz, Butterworth three-order filter
% fType = 'butter';
% BPFfreq = [2 100];
% filterOrder = 4;
% % Band-pass filter configuration
% [z, p, k] = temporalBPFconfig(fType, ECG_freq, BPFfreq, filterOrder);
% ECGfilt = temporalBPFrun(ECG2, z, p, k);

%%  Optimal Minimum Order Designs
% Fp1 — frequency at the edge of the start of the pass band. Specified in
% normalized frequency units. Also called Fpass1. }
Fp1 = 0.03*2/resp_freq;
% Fp2 — frequency at the edge of the end of the pass band. Specified in
% normalized frequency units. Also called Fpass2.
Fp2 = 28*2/resp_freq;
% Fst1 — frequency at the edge of the start of the first stop band.
% Specified in normalized frequency units. Also called Fstop1.
Fst1 = 0.01*2/resp_freq;
% Fst2 — frequency at the edge of the start of the second stop band.
% Specified in normalized frequency units. Also called Fstop2.
Fst2 = 30*2/resp_freq;
% Ap — amount of ripple allowed in the pass band. Also called Apass.
Ap = 1;
% Ast1 — attenuation in the first stop band in decibels (the default
% units). Also called Astop1.
Ast1 = 30;
Ast2 = Ast1;
Hf = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2);
% Equiripple designs result in the lowpass
% filter with the smallest possible order to meet a set of specifications.
Hd = design(Hf,'equiripple','systemobject',true);
measure(Hd)
hfvt = fvtool(Hd,'Color','White');
legend(hfvt,'Equiripple design')

% -------------------- Cosmetic changes -----------------------------------
% My Red Color
myColor = [204 0 0]/255;
% set(gca,'FontSize', 12)
legend({'Equiripple design'}, 'FontSize', 14)
figureHandle = gcf;
%# make all text in the figure to size 14 and bold
set(findall(figureHandle,'type','text'),'fontSize',14)
gcachild = get(gca,'Children');
set(gca, 'Xcolor', 'k');
set(gca, 'Ycolor', 'k');
set(gcachild, 'Linewidth', 2);
% Change plot marker
hchildren = get(hfvt,'children');
haxes = hchildren(strcmpi(get(hchildren,'type'),'axes'));
hline = get(haxes,'children');
set(hline(1), 'Color', [0 0 204]/255, 'LineWidth',2); 
set(hline(2), 'Color',  myColor, 'LineWidth',2);
% set(hline(3), 'Color', 'k', 'LineWidth',2);

%% Apply filter to the signal
b =  get(Hd,'Numerator');   % N = 75
a = 1;  
respFilt = filtfilt(b, a, physio.Resp);

%% Plot filter and unfiltered signals
figure;
plot(t,physio.Resp)
hold on
plot(t,respFilt)
xlim([115 118]);

%% Compute resampling ratio
% Free up memory
clear b physio
pack
% resampling at p/q times the original sampling rate
q = round(resp_freq);
p = round(Fs);
% Cut data to experiment length
respFilt = respFilt(1:size(y,2)*resp_freq+1);
t = t(1:size(y,2)*resp_freq+1);

%% Memory-mapping
tic
fmem_name=fullfile(pathNameParent, [IOI.subj_name '_resampled.dat']);
fid = fopen(fmem_name,'w');
fwrite(fid, zeros([size(y,1) size(y,2)*resp_freq+1]), 'int16');
fclose(fid);
im_obj = memmapfile(fmem_name,...
        'Format',{'double' [size(y,1) size(y,2)*resp_freq+1] 'data'}, ...
        'writable',true);
toc    
%% Interpolate & Resample
% [ECGresamp,b] = resample(ECGfilt,p,q);
[yResamp,b] = resample(y,q,p);

% Concatenate data
y = [respFilt'; y];
% y = y'; % [nChannels x nTimePoints] e.g. [5562 x 360]
%% Write data to be imported in Analyzer2
fwrite_NIR(fullfile(pathNameParent,[IOI.subj_name '_ECG_noCB.nir']), y);

% EOF