%% ECG example
addpath(genpath('C:\spm8\toolbox\ioi11')); clc
% 1. Crear un vector tiempo de n puntos muestreado a 1000 Hz
Fs = 1000; %sampling rate
Ts = 1/Fs; %sampling time interval
% Change to signal of interest (ECG1-ECG4)
load('C:\Edgar\Data\OIS_Data\12_10_18,NC09\ECG.mat');
% Read ECG as column (only first 360 seconds)
ECGsignal = ECG(1:360000);
% Convert to column
ECGsignal = ECGsignal';
clear ECG;
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
BPFfreq = [0.2 100];
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
ylabel('Heart rate (Hz)');
xlabel('t (s)')
toc
% EOF

