% For HbO
dataName = 'C:\Edgar\Data\OIS_Data\25Fev\E02';
resultsName = 'C:\Edgar\Data\OIS_Results';
load(fullfile(dataName,'Dim_binFile.mat'));
fileID=fopen(fullfile(dataName,'HbO.bin'),'r');
HbO=fread(fileID,'int32');
HbO=reshape(HbO,Temps_d1,X_d2,Y_d3); %Temps_d1,… are stored in Dim_binFile.mat
