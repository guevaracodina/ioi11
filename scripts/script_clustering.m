%% Clustering/linkage script
addpath(genpath('D:\spm8\toolbox\ioi11'))
% loading data
tic
IOImat = 'E:\Edgar\Data\IOS_Results\12_07_24,RS01\IOI.mat';
load(IOImat)
% Session number
s1 = 1;
% HbO (5) HbR (6) Flow (7)
c1 = 6;
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
% tic
% fprintf('Loading SVD data for subject %s, session: %02d (%s)...\n',...
%         IOI.subj_name,s1,colorNames{1+c1})
% pathName = fullfile(IOI.dir.dir_subj_res, 'SVD');
% load(fullfile(pathName,sprintf('SVD_S%02d_C%d.mat',s1,c1)),'U','S','V');
% disp(['Data loaded in: ' datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS')]);

[pathName, ~, ~] = fileparts(IOImat);
pathName = fullfile(pathName, 'SVD');
if ~exist(pathName, 'file'), mkdir(pathName); end

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
nClusters       = 12;
% Number of iterations
nRep            = 10;
% Type of distance metric
distanceType    = 'correlation';
% Linkage method
linkageMethod   = 'weighted';
% Number of independent runs
nIter = 1;

% IDX = zeros([size(yMatrix,1) nIter]);
% groupsK = zeros([nIter, nClusters]);
% 
% for i = 1:nIter,
%     fprintf('Computing k-means clustering for subject %s, session: %02d (%s)...\n',...
%         IOI.subj_name,s1,colorNames{1+c1})
%     tic
%     IDX(:,i) = kmeans(yMatrix, nClusters, 'distance', distanceType, 'replicates', nRep);
%     IDX(:,i) = ioi_sort_clusters(IDX(:,i));
%     for j = 1:nClusters,
%         groupsK(i,j)=numel(find(IDX(:,i)==j));
%     end
%     groupsK(i,:) = sort(groupsK(i,:));
%     fprintf('k-means iter%d computed in: %s\n', i, datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS'));
% end

%% Display k-means clustering
% h = figure; 
% set(gcf,'color','w')
% set(gcf,'name',sprintf('%s S%02d C%d (%s) k-means\n',...
%     IOI.subj_name,s1,c1,colorNames{1+c1}))
% colormap([0 0 0; hsv(nClusters)]); 
% for i = 1:nIter,
%     yTmp = zeros(size(ySlice));
%     yTmp(brainMask) = IDX(:,i);
%     subplot(ceil(sqrt(nIter)), round(nIter/sqrt(nIter)), i)
%     imagesc(yTmp); axis image; axis off
%     title(sprintf('Run: %d',i))
%         if i == nIter
%             colorbar
%         end
% end
% newName = sprintf('%s_S%02d_C%d_(%s)_kMeansRuns.PNG',IOI.subj_name,s1,c1,colorNames{1+c1});
%     % Save as PNG
% print(h, '-dpng', '-r300', fullfile(pathName,newName));
    

%% Using clusterdata (hierarchical clustering)
% (performs all of the necessary steps for you. You do not need to execute the
% pdist, linkage, or cluster functions separately.)
T2 = zeros([size(yMatrix,1) nIter]);
Z2 = zeros([size(yMatrix,1)-1 3 nIter]);
groupsClusterData = zeros([nIter, nClusters]);


for i = 1:nIter,
    fprintf('Computing hierarchical clustering for subject %s, session: %02d (%s)...\n',...
        IOI.subj_name,s1,colorNames{1+c1})
    [T2(:,i), Z2(:,:,i)] = ioi_clusterdata(yMatrix,'distance', distanceType,...
        'linkage', linkageMethod, 'maxclust', nClusters);
    T2(:,i) = ioi_sort_clusters(T2(:,i));
    for j = 1:nClusters,
        groupsClusterData(i,j) = numel(find(T2(:,i)==j));
    end
    groupsClusterData(i,:) = sort(groupsClusterData(i,:));
    fprintf('clusterData iter%d computed in: %s\n', i, datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS'));
end

%% Display dendogram
% h = figure;
% set(h, 'color', 'w')
% [H2, T3] = dendrogram(Z2, nClusters);
% title('Dendrogram','interpreter', 'none', 'FontSize', 14)
% set(H2, 'LineWidth', 2);
% set(gca, 'FontSize', 12);
% xlabel('Cluster', 'FontSize', 12)
% ylabel('Distance', 'FontSize', 12)
% newName = sprintf('%s_S%02d_C%d_(%s)_HierarClustDendogram.PNG',IOI.subj_name,s1,c1,colorNames{1+c1});
% % Save as PNG
% print(h, '-dpng', '-r300', fullfile(pathName,newName));

%% Display hierarchical clustering
newName = sprintf('%s_S%02d_C%d_(%s)_HierarClustRuns.PNG',IOI.subj_name,s1,c1,colorNames{1+c1});
h = figure; 
set(h,'color','w')
set(h,'name',sprintf('%s S%02d C%d (%s) hierarchical clustering\n',...
    IOI.subj_name,s1,c1,colorNames{1+c1}))
colormap([0 0 0; hsv(nClusters)]); 
for i = 1:nIter,
    yTmp = zeros(size(ySlice));
    yTmp(brainMask) = T2(:,i);
    subplot(ceil(sqrt(nIter)), round(nIter/sqrt(nIter)), i)
    imagesc(yTmp); axis image; axis off
    if nIter == 1
        title(newName, 'interpreter', 'none');
    else
        title(sprintf('Run: %d',i))
    end
    if i == nIter
        colorbar
    end
end

% Save as PNG
print(h, '-dpng', '-r300', fullfile(pathName,newName));

% yTmp = zeros(size(ySlice));
% yTmp(brainMask) = ioi_sort_clusters(T3);
% figure; set(gcf,'color','w')
% imagesc(yTmp); axis image; colormap([0 0 0; hsv(nClusters)]); colorbar
% title(sprintf('%s S%02d C%d (%s) dendrogram recomputed from ioi_clusterdata\n',...
%         IOI.subj_name,s1,c1,colorNames{1+c1}),...
%         'interpreter', 'none', 'FontSize', 14)
%     


%% Correlation between pixels and 10 firstSVD
% fprintf('Computing pixel/SVD correlation for subject %s, session: %02d (%s)...\n',...
%         IOI.subj_name,s1,colorNames{1+c1})
% tic
% % the first initial condition was derived from the first ten singular modes of
% % the connectivity matrix (after SVD) with pixels assigned to whichever singular
% % vector with which they have the highest positive/negative coefficient.
% % not sure of this part??
% corrPixSVD = corr(yMatrix, U*S*V');
% fprintf('Pixel-SVD correlation computed in: %s\n', datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS'));
% fprintf('Writing Pixel-SVD correlation matrix to disk for subject %s, session: %02d (%s)...\n',...
%     IOI.subj_name,s1,colorNames{1+c1})
% tic
% save(fullfile(pathName,sprintf('corrPixSVD_S%02d_C%02d.mat',s1,c1)),'corrPixSVD');
% fprintf('Data written to disk in: %s\n', datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS'));
