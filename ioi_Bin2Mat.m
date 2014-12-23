%Read ECG data and convert to mV
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

