%% Read HbO Data
dataName = 'D:\Edgar\OIS_Data\FB25E02';
resultsName = 'D:\Edgar\OIS_Results';
subjectName = 'FB25E02';
load(fullfile(dataName,'Dim_binFile.mat'));
fileID=fopen(fullfile(dataName,'HbO.bin'),'r');
HbO = fread(fileID,'int32');
HbO = reshape(HbO,Temps_d1,X_d2,Y_d3); %Temps_d1,… are stored in Dim_binFile.mat

%%
maxVal = max(HbO(:));
minVal = min(HbO(:));
figure;
for iTime=1:size(HbO,1)
    imagesc(squeeze(HbO(iTime,:,:)),[minVal, maxVal]);
    axis image
    drawnow;
    pause(0.1)
end

%% Anatomical
% Bregma est environ à la ligne 250, colonne 400 et Lambda: 250, 150.
anat = mat2gray(squeeze(mean(HbO,1)));
imwrite(anat, fullfile(resultsName,[subjectName '_anat_S01.png']),...
    'BitDepth',16)

%% Brain Mask
brainMask = im2bw(anat,1-eps);
imwrite(brainMask, fullfile(resultsName,[subjectName '_brainMask.png']),...
    'BitDepth',1)