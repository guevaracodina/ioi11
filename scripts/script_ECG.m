%% ECG example
addpath(genpath('C:\spm8\toolbox\ioi11')); clc
% 1. Crear un vector tiempo de n puntos muestreado a 1000 Hz
Fs = 1000; %sampling rate
Ts = 1/Fs; %sampling time interval

% Change to signal of interest (ECG2-ECG4)
% ECG.mat is just an example from mice
load('C:\Edgar\Data\OIS_Data\12_10_18,NC09\ECG.mat');

% Read ECG data and convert to mV
% Change folder as needed
% pathName = 'C:\Edgar\Data\IOIData20141211\Physio Monitoring\104650';
pathName = 'C:\Edgar\Data\IOIData20141211\Physio Monitoring\121423'; % looks promising!

cd (pathName)
fid=fopen('ECG1.bin','r');
data = uint8(fread(fid));
ECG1 = double(swapbytes(typecast(data,'int32')))/8388608*2400;
fclose(fid);
fid=fopen('ECG2.bin','r');
data = uint8(fread(fid));
ECG2 = double(swapbytes(typecast(data,'int32')))/8388608*2400;
fclose(fid);
fid=fopen('ECG3.bin','r');
data = uint8(fread(fid));
ECG3 = double(swapbytes(typecast(data,'int32')))/8388608*2400;
fclose(fid);
fid=fopen('ECG4.bin','r');
data = uint8(fread(fid));
ECG4 = double(swapbytes(typecast(data,'int32')))/8388608*2400;
fclose(fid);

%Read resp data
fid=fopen('resp.bin','r');
data = uint8(fread(fid));
Resp = double(swapbytes(typecast(data,'int16')));
fclose(fid);

%Read SpO2 data
fid=fopen('SPO2_IR.bin','r');
data = uint8(fread(fid));
SPO2_IR = double(swapbytes(typecast(data,'int32')));
fclose(fid);
fid=fopen('SPO2_R.bin','r');
data = uint8(fread(fid));
SPO2_R = double(swapbytes(typecast(data,'int32')));
fclose(fid);

%Read temperature data
fid=fopen('temperature1.bin','r');
data = uint8(fread(fid));
Temp1 = double(swapbytes(typecast(data,'int16')));
B25_50 = 3380; voltage = Temp1./25./2.^10.*3.3/2; res = voltage./((3.3-voltage)./10e3); num = (log(10e3)-log(res));
Temp1 = ( 1./(273.15+25) - num./B25_50 ).^-1 - 273.15;  
fclose(fid);
fid=fopen('temperature2.bin','r');
data = uint8(fread(fid));
Temp2 = double(swapbytes(typecast(data,'int16')));
B25_50 = 3950; voltage = Temp2./25./2.^10.*3.3/2; res = voltage./((3.3-voltage)./10e3); num = (log(10e3)-log(res));
Temp2 = ( 1./(273.15+25) - num./B25_50 ).^-1 - 273.15;  
fclose(fid);
fid=fopen('temperature3.bin','r');
data = uint8(fread(fid));
Temp3 = double(swapbytes(typecast(data,'int16')));
B25_50 = 3380; voltage = Temp3./25./2.^10.*3.3/2; res = voltage./((3.3-voltage)./10e3); num = (log(10e3)-log(res));
Temp3 = ( 1./(273.15+25) - num./B25_50 ).^-1 - 273.15;  
fclose(fid);

%% Choose ECG signal ECG or ECG2/3/4
% Read sample ECG as column (only first 360 seconds)
% ECGsignal = ECG(1:360000);      % Example ECG
ECGsignal = ECG2;             % Rat pups ECG (ECG2/3/4)
% Convert to column
ECGsignal = ECGsignal';
% clear ECG;

%% FFT of ECG signal
n = length(ECGsignal); %number of samples
t = (0:n-1)*Ts; %time vector
N=length(ECGsignal);
%this part of the code generates that frequency axis
if mod(N,2)==0
    k=-N/2:N/2-1; % N even
else
    k=-(N-1)/2:(N-1)/2; % N odd
end
T=N/Fs;
freq=k/T;  %the frequency axis

figure; set(gcf, 'color', 'w')
subplot(221)
plot(t,ECGsignal)
title('Raw ECG')
xlabel('t (s)')
xlim([202 203.5])

%% takes the fft of the signal, and adjusts the amplitude accordingly
ECGfft=fft(ECGsignal)/N; % normalize the data
ECGfft=fftshift(ECGfft); %shifts the fft data so that it is centered
subplot(223)
plot(freq,abs(ECGfft))
title('FFTSHIFT(FFT(ECGsignal)))')
xlabel('f (Hz)')
xlim([0 Fs/2]);

%% Filter
fType = 'butter';
BPFfreq = [2 40]; % BPFfreq = [0.2 100];
filterOrder = 4;
% Band-pass filter configuration
[z, p, k] = temporalBPFconfig(fType, Fs, BPFfreq, filterOrder);
ECGfilt = temporalBPFrun(ECGsignal, z, p, k);
subplot(222)
plot(t, ECGfilt)
title('Filtered ECG')
xlabel('t (s)')
xlim([202 203.5])

%% FFT of filtered ECG
ECGsignal = ECGfilt;
n = length(ECGsignal); %number of samples
t = (0:n-1)*Ts; %time vector
N=length(ECGsignal);
%this part of the code generates that frequency axis
if mod(N,2)==0
    k=-N/2:N/2-1; % N even
else
    k=-(N-1)/2:(N-1)/2; % N odd
end
T=N/Fs;
freq=k/T;  %the frequency axis

%takes the fft of the signal, and adjusts the amplitude accordingly
ECGfft=fft(ECGsignal)/N; % normalize the data
ECGfft=fftshift(ECGfft); %shifts the fft data so that it is centered
subplot(224)
plot(freq,abs(ECGfft))
title('FFTSHIFT(FFT(Filtered)))')
xlabel('f (Hz)')
xlim(BPFfreq)

%% Retrieving bpm
tic
% addpath(genpath('C:\spm8\toolbox\oct12'))
% Wrapper function to use OCT tool subfunction_find_ECG_peaks
acqui_info.ecg_signal = ECGsignal;

% ECG sampling period in micro-seconds
acqui_info.line_period_us = Ts*1e6;

% Number of A-lines per ramp
acqui_info.ramp_length = 100;

% Do not display plot
display_plots = false;

% Compute instantaneous heart rate
[~, bpm, qrs] = ioi_find_ECG_peaks(acqui_info,display_plots);

% Display heart rate
figure; set(gcf, 'color', 'w')
plot(qrs.peaks_position, qrs.peaks_rate)
% Save qrs to a file
ylabel('Heart rate (Hz)');
xlabel('t (s)')
toc
% EOF

