%% Read BLK files in acquisition folder
clear; clc
% Info parameters
verbose = true;
% Pre-processing parameters
% decoupage temporel [MIN MAX], si [0] = pas de decoupage
bz = 0;
% decoupage spatial [xMIN xMAX yMIN yMAX], si [0] = pas de decoupage
bxy = 0;
% binning temporel, si [1] = pas de binning
r = 1;
% soustraction de la ligne de base 
base = -1;
% filtrage spatial 
filtre = 0;
% limite des fréquences conservées lors de la Transformée de Fourier
FM = 0;
% nombre d'étapes de traitement
N = 1;
% sens du signal (-1 = neg)
sens = -1;
% 1 si mise en pourcentage, 0 autrement
pourcent = 1;
% Acquisition path
acqPath = 'F:\Edgar\Data\Visual_RS_Data\13_10_08,RS01\RS01E04';
%% Reading file
dirStruct = dir(fullfile(acqPath,'*.blk'));
for iFiles = 1:numel(dirStruct),
    fileName = fullfile(acqPath, dirStruct(iFiles).name);
    [width(iFiles),height(iFiles),nFramesPerCond(iFiles),nConditions(iFiles),...
        nTrialsPerBlock(iFiles),nBlocks(iFiles),nRepetitions(iFiles),...
        segmentLength(iFiles),dataType(iFiles),fileNameOut{iFiles}] = ioi_read_blk_info(fileName, verbose);
    [currAcqPath, currFileNameOut, currExtension] = fileparts(fileNameOut{iFiles});
    cible{iFiles} = '';
    rawImage = ioi_read_blk_images(fileNameOut{iFiles},bz,bxy,r,base,filtre,FM,N,cible{iFiles},sens,pourcent);
end
figure; imagesc(squeeze(rawImage(:,:,1))); colorbar
%% 

