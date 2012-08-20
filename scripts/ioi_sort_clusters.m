function IDXout = ioi_sort_clusters(IDX)
% Sorts the cluster indices in ascending order by the number of pixels belonging
% to each cluster. Helps visualizing independent runs of clustering algorithms.
% SYNTAX
% IDXout = ioi_sort_clusters(IDX)
% INPUT
% IDX       Contains the cluster indices of each pixel
% OUTPUT
% IDXout    Cluster indices sorted by ascending number of pixels
%_______________________________________________________________________________
% Copyright (C) 2012 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________

% Preallocate memory
IDXout = zeros(size(IDX));
% Find the number of clusters
nClusters = numel(unique(IDX));
% Preallocate memory
pixUnsorted = zeros([nClusters 1]);

% Find number of pixels per cluster
for i = 1:nClusters,
    pixUnsorted(i) = numel(find(IDX==i));
end

% Sort clusters in ascending order
[~, sortedIndex]  = sort(pixUnsorted);

% Replace cluster indices by their sorted equivalent
for iClusters = 1:nClusters,
    IDXout(IDX == iClusters) = sortedIndex(iClusters);
end

% EOF
