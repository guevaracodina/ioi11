%Read ECG data and convert to mV
pathName = 'C:\Edgar\Data\IOIData20141023\Physio Monitoring\103931_K01';
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

%% Plot FFT
% 1. 1.	Crear un vector tiempo de n puntos muestreado a 1000 Hz
Fs = 1000; %sampling rate
Ts = 1/Fs; %sampling time interval
% Change to signal of interest (ECG1-ECG4)
ECGsignal = ECG3;
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
Y3=fft(ECGsignal)/N; % normalize the data
Y3=fftshift(Y3); %shifts the fft data so that it is centered
subplot(211)
plot(freq,abs(Y3))
title('FFTSHIFT(FFT(ECGsignal)))')
xlabel('f (Hz)')

%% Filter
fType = 'butter';
BPFfreq = [0.2 400];
filterOrder = 4;
% Band-pass filter configuration
[z, p, k] = temporalBPFconfig(fType, Fs, BPFfreq, filterOrder);
Yfilt = temporalBPFrun(ECGsignal, z, p, k);

%% FFT of filtered ECG
ECGsignal = Yfilt;
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
Y3=fft(ECGsignal)/N; % normalize the data
Y3=fftshift(Y3); %shifts the fft data so that it is centered
subplot(212)
plot(freq,abs(Y3))
title('FFTSHIFT(FFT(Filtered)))')
xlabel('f (Hz)')


