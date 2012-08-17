%% Singular Value Decomposition script
% loading data
tic
IOImat = 'E:\Edgar\Data\IOS_Results\12_07_17,NC01\IOI.mat';
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
load(fullfile(pathName,sprintf('SVD_S%02d_C%d.mat',s1,c1)),'U','S','V');
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
% Number of clusters to form
nClusters       = 10;
% Number of iterations
nRep            = 10;
% Type of distance metric
distanceType    = 'correlation';
fprintf('Computing k-means clustering for subject %s, session: %02d (%s)...\n',...
        IOI.subj_name,s1,colorNames{1+c1})
tic
IDX = kmeans(yMatrix, nClusters, 'distance', distanceType, 'replicates', nRep);
fprintf('k-means computed in: %s\n', datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS'));

%% Display k-means clustering
yTmp = zeros(size(ySlice));
yTmp(brainMask) = IDX;
figure; 
imagesc(yTmp); axis image; colormap(jet(nClusters)); colorbar
title(sprintf('%s S%02d C%d (%s) k-means\n',...
        IOI.subj_name,s1,c1,colorNames{1+c1}),...
        'interpreter', 'none', 'FontSize', 14)

%% linkage
% Compute pairwise distance first
D = pdist(yMatrix, distanceType);
linkageMethod = 'weighted';
% Perform linkage
Z = linkage(D, linkageMethod);
figure;
% Display binary tree
[H, T] = dendrogram(Z, nClusters);
title('Dendrogram computed from linkage','interpreter', 'none', 'FontSize', 14)
set(H, 'LineWidth', 2);

%% Using clusterdata 
% (performs all of the necessary steps for you. You do not need to execute the
% pdist, linkage, or cluster functions separately.)
addpath(genpath('D:\spm8\toolbox\ioi11'))
[T2, Z2] = ioi_clusterdata(yMatrix,'distance', distanceType, 'linkage',...
    linkageMethod, 'maxclust', nClusters);
figure;
[H2, T3] = dendrogram(Z2, nClusters);
title('Dendrogram recomputed from ioi_clusterdata','interpreter', 'none', 'FontSize', 14)
set(H2, 'LineWidth', 2);

%% Display hierarchical clustering
yTmp = zeros(size(ySlice));
yTmp(brainMask) = T;
figure; 
imagesc(yTmp); axis image; colormap(jet(nClusters)); colorbar
title(sprintf('%s S%02d C%d (%s) original dendrogram\n',...
        IOI.subj_name,s1,c1,colorNames{1+c1}),...
        'interpreter', 'none', 'FontSize', 14)

yTmp = zeros(size(ySlice));
yTmp(brainMask) = T2;
figure; 
imagesc(yTmp); axis image; colormap(jet(nClusters)); colorbar
title(sprintf('%s S%02d C%d (%s) ioi_clusterdata\n',...
        IOI.subj_name,s1,c1,colorNames{1+c1}),...
        'interpreter', 'none', 'FontSize', 14)
    
yTmp = zeros(size(ySlice));
yTmp(brainMask) = T3;
figure; 
imagesc(yTmp); axis image; colormap(jet(nClusters)); colorbar
title(sprintf('%s S%02d C%d (%s) dendrogram recomputed\n',...
        IOI.subj_name,s1,c1,colorNames{1+c1}),...
        'interpreter', 'none', 'FontSize', 14)
    


%% Correlation between pixels and 10 firstSVD
fprintf('Computing pixel/SVD correlation for subject %s, session: %02d (%s)...\n',...
        IOI.subj_name,s1,colorNames{1+c1})
tic
% the first initial condition was derived from the first ten singular modes of
% the connectivity matrix (after SVD) with pixels assigned to whichever singular
% vector with which they have the highest positive/negative coefficient.
% not sure of this part??
corrPixSVD = corr(yMatrix, U*S*V');
fprintf('Pixel-SVD correlation computed in: %s\n', datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS'));
fprintf('Writing Pixel-SVD correlation matrix to disk for subject %s, session: %02d (%s)...\n',...
    IOI.subj_name,s1,colorNames{1+c1})
tic
save(fullfile(pathName,sprintf('corrPixSVD_S%02d_C%02d.mat',s1,c1)),'corrPixSVD');
fprintf('Data written to disk in: %s\n', datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS'));
