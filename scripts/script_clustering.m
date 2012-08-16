%% Singular Value Decomposition script
% loading data
tic
IOImat = 'E:\Edgar\Data\IOS_Results\12_07_17,CC01\GLMfcIOS\IOI.mat';
load(IOImat)
% Session number
s1 = 1;
% HbO (5) HbR (6) Flow (7)
c1 = 5;
% Get colors to iIOI.fcIOS.masknclude information
colorNames = fieldnames(IOI.color);

% Check if regression was correctly performed
if IOI.fcIOS.SPM.wholeImageRegressOK{s1, c1}
    fprintf('Loading NIfTI data for subject %s, session: %02d (%s)...\n',...
        IOI.subj_name,s1,colorNames{1+c1})
    % Load filtered, downsampled & regressed image time-course
    vol = spm_vol(IOI.fcIOS.SPM.fname{s1, c1});
    y = spm_read_vols(vol);
else
    fprintf('No regressed time-series for subject %s, session: %02d (%s)\n',...
        IOI.subj_name,s1,colorNames{1+c1})
end

% Check if brain mask was performed
if IOI.fcIOS.mask.maskOK
    % Load brain mask
    brainMaskVol = spm_vol(IOI.fcIOS.mask.fname);
    brainMask = logical(spm_read_vols(brainMaskVol));
    if size(brainMask,1)~= size(y,1)|| size(brainMask,2)~= size(y,2)
        brainMask = ioi_MYimresize(brainMask, [size(y,1) size(y,2)]);
    end
else
    fprintf('No regressed time-series for subject %s\n',IOI.subj_name)
end

% Load SVD
tic
fprintf('Loading SVD data for subject %s, session: %02d (%s)...\n',...
        IOI.subj_name,s1,colorNames{1+c1})
pathName = fullfile(IOI.dir.dir_subj_res, 'SVD');
load(fullfile(pathName,sprintf('SVD_S%02d_C%02d.mat',s1,c1)),'U','S','V');
disp(['Data loaded in: ' datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS')]);

%% Preallocate memory and make 2D matrix from 4d tensor
yMatrix = zeros([nnz(brainMask) size(y,4)]);
for iCols = 1:size(yMatrix,2)
    ySlice = squeeze(y(:,:,1,iCols));
    yMatrix(:, iCols) = ySlice(brainMask);
end
% cleanup
clear y vol

%% k-means clustering
% nClusters = size(U,2);
nClusters = 6;
nRep = 10;
fprintf('Computing k-means clustering for subject %s, session: %02d (%s)...\n',...
        IOI.subj_name,s1,colorNames{1+c1})
tic
IDX = kmeans(yMatrix, nClusters, 'distance', 'correlation', 'replicates', nRep);
fprintf('k-means computed in: %s\n', datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS'));

%% Display clustering
yTmp = zeros(size(ySlice));
yTmp(brainMask) = IDX;
figure; 
imagesc(yTmp); axis image; colormap(jet(nClusters)); colorbar

%% Correlation between pixels and 10 firstSVD
fprintf('Computing pixel/SVD correlation for subject %s, session: %02d (%s)...\n',...
        IOI.subj_name,s1,colorNames{1+c1})
tic
corrPixSVD = corr(yMatrix, U);
fprintf('Pixel-SVD correlation computed in: %s\n', datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS'));
fprintf('Writing Pixel-SVD correlation matrix to disk for subject %s, session: %02d (%s)...\n',...
    IOI.subj_name,s1,colorNames{1+c1})
tic
save(fullfile(pathName,sprintf('corrPixSVD_S%02d_C%02d.mat',s1,c1)),'corrPixSVD');
disp(['Data written to disk in: ' datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS')]);
