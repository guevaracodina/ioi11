%% Singular Value Decomposition script
% loading data
tic
IOImat = 'E:\Edgar\Data\IOS_Results\12_07_17,CC01\GLMfcIOS\IOI.mat';
load(IOImat)
% Session number
s1 = 2;
% HbO (5) HbR (6) Flow (7)
c1 = 7;
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
disp(['Data loaded in: ' datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS')]);

%% Test: Create matrix with only the brain pixels
% % ySlice is the first slice of the time-series
% ySlice = squeeze(y(:,:,1,1));
% figure; imagesc(ySlice);
% % yindex is a column vector containing only non-masked pixels of ySlice
% yindex = ySlice(brainMask);
% % ynew will contain the pixels from yindex in their right place
% ynew = zeros(size(ySlice));
% ynew(brainMask) = yindex;
% figure; imagesc(ynew);

%% Preallocate memory and make 2D matrix from 4d tensor
yMatrix = zeros([nnz(brainMask) size(y,4)]);
for iCols = 1:size(yMatrix,2)
    ySlice = squeeze(y(:,:,1,iCols));
    yMatrix(:, iCols) = ySlice(brainMask);
end
% cleanup
clear y vol

%% Create NxN connectivity (auto-correlation) matrix 
% (N = number of pixels defined as brain)
% Size of rMatrix in Gb  => 8*nnz(brainMask)^2 / 2^30; (~5.75 Gb)
tic
fprintf('Computing connectivity matrix for subject %s, session: %02d (%s)...\n',...
    IOI.subj_name,s1,colorNames{1+c1})
[rMatrix] = corrcoef(yMatrix');
disp(['Connectivity (auto-correlation) matrix computed in: ' datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS')]);
[pathName, ~, ~] = fileparts(IOImat);
pathName = fullfile(pathName, 'SVD');
if ~exist(pathName, 'file'), mkdir(pathName); end
% cleanup
clear yMatrix
fprintf('Writing connectivity matrix to disk for subject %s, session: %02d (%s)...\n',...
    IOI.subj_name,s1,colorNames{1+c1})
tic
save(fullfile(pathName,sprintf('corrMatrix_S%02d_C%d.mat',s1,c1)),'rMatrix');
disp(['Data written to disk in: ' datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS')]);

%% PCA on correlation matrix
% not really necessary, equivalent to SVD
% tic
% fprintf('Computing PCA for subject %s, session: %02d (%s)...\n',...
%     IOI.subj_name,s1,colorNames{1+c1})
% [coeff, latent, explained] = pcacov(rMatrix);
% fprintf('PCA computed in: %s \n', datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS'));
% tic
% save(fullfile(pathName,sprintf('PCA_C%2d.mat',c1)),'coeff','latent','explained');
% fprintf('Data written to disk in: %s \n', datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS'));

%% PCA display

% % Preallocate
% yRGB = zeros([size(ySlice), 3]);
% yTmp = zeros(size(ySlice));
% 
% % First 3PC
% pcaR = coeff(:,1);
% pcaG = coeff(:,2);
% pcaB = coeff(:,3);
% 
% % Indexing
% yTmp(brainMask) = pcaR;
% yRGB(:,:,1) = yTmp;
% yTmp(brainMask) = pcaG;
% yRGB(:,:,2) = yTmp;
% yTmp(brainMask) = pcaB;
% yRGB(:,:,3) = yTmp;
% 
% % Mask out non-brain voxels
% brainMaskRep = repmat(brainMask,[1 1 3]);
% yRGB = yRGB .* double(brainMaskRep);
% 
% % PCA
% for iPCA = 1:3,
%     h = figure;
%     imagesc(squeeze(yRGB(:,:,iPCA))); axis image; colorbar
%     title(sprintf('%s S%02d C%d (%s) PC%d\n',...
%         IOI.subj_name,s1,c1,colorNames{1+c1},iPCA),...
%         'interpreter', 'none', 'FontSize', 14)
%     set(h,'color','w')
%     newName = sprintf('%s_S%02d_C%d_(%s)_PC%02d.PNG',IOI.subj_name,s1,c1,colorNames{1+c1},iPCA);
%     % Save as PNG
%     print(h, '-dpng', '-r300', fullfile(pathName,newName));
%     close(h)
% end

%% Load correlation matrix
% tic
% load(fullfile(pathName,sprintf('corrMatrix_S%02d_C%d.mat',s1,c1)),'rMatrix');
% disp(['Correlation matrix loaded in: ' datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS')]);
% cleanup
clear y vol

%% Create memory/mapped objects
% only the nSVD largest vectors
nSVD = 10;

% clear im_obj*
% delete(fullfile(IOI.dir.dir_subj_res, 'svd*.dat'))

% %create large memmapfile for S(several GB)
% fmemS_name = fullfile(IOI.dir.dir_subj_res, 'svdS.dat');
% fid = fopen(fmemS_name,'w');
% fwrite(fid, zeros([nSVD nSVD]), 'double');
% fclose(fid);
% im_objS = memmapfile(fmemS_name,'format',...
%     {'double' [nSVD nSVD] 'svdS'},...
%     'writable',true);
% 
% %create large memmapfile for U(several GB)
% fmemU_name = fullfile(IOI.dir.dir_subj_res, 'svdU.dat');
% fid = fopen(fmemU_name,'w');
% fwrite(fid, zeros([nnz(brainMask) nSVD]), 'double');
% fclose(fid);
% im_objU = memmapfile(fmemU_name,'format',...
%     {'double' [nnz(brainMask) nSVD] 'svdU'},...
%     'writable',true);
% 
% %create large memmapfile for V(several GB)
% fmemV_name = fullfile(IOI.dir.dir_subj_res, 'svdV.dat');
% fid = fopen(fmemV_name,'w');
% fwrite(fid, zeros([nnz(brainMask) nSVD]), 'double');
% fclose(fid);
% im_objV = memmapfile(fmemV_name,'format',...
%     {'double' [nnz(brainMask) nSVD] 'svdV'},...
%     'writable',true);

%create large memmapfile for rMatrix (several GB)
fmemR_name = fullfile(IOI.dir.dir_subj_res, 'svdR.dat');
fid = fopen(fmemR_name,'w');
fwrite(fid, rMatrix, 'double');
fclose(fid);
im_objR = memmapfile(fmemR_name,'format',...
    {'double' size(rMatrix) 'rMatrix'},...
    'writable', false);
% cleanup
clear rMatrix
% pack

%% Compute largest SVD vectors
% rMatrix = USV'
tic
fprintf('Computing the %d largest SV...\n', nSVD);
%[im_objU.Data.svdU im_objS.Data.svdS im_objV.Data.svdV] = svds(im_objR.Data.rMatrix, nSVD);
[U S V] = svds(im_objR.Data.rMatrix, nSVD);
fprintf('Largest %d SVD computed in: %s\n', nSVD, datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS'));
% pathName = fullfile(IOI.dir.dir_subj_res, 'SVD');
% if ~exist(pathName, 'file'), mkdir(pathName); end
fprintf('Writing SVD data to disk for subject %s, session: %02d (%s)...\n',...
    IOI.subj_name,s1,colorNames{1+c1})
tic
save(fullfile(pathName,sprintf('SVD_S%02d_C%d.mat',s1,c1)),'U','S','V');
disp(['Data written to disk in: ' datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS')]);

% cleanup
clear im_obj*
delete(fullfile(IOI.dir.dir_subj_res, 'svd*.dat'))

%% False RGB image created from first 3 SV
iSVD = 1;
%svdR = U(:,iSVD)*S(iSVD,:)*V(iSVD,:)';
svdR = U(:,iSVD);
iSVD = 2;
%svdG = U(:,iSVD)*S(iSVD,:)*V(iSVD,:)';
svdG = U(:,iSVD);
iSVD = 3;
%svdB = U(:,iSVD)*S(iSVD,:)*V(iSVD,:)';
svdB = U(:,iSVD);
iSVD = 4;
%svdB = U(:,iSVD)*S(iSVD,:)*V(iSVD,:)';
svdY = U(:,iSVD);

%% Arrange the first 4 singular modes in their correspondent place
% Preallocate
yRGB = zeros([size(ySlice), 4]);
yTmp = zeros(size(ySlice));
% Normalize each channel between 0 and 1
% svdR = svdR - min(svdR);
% svdR = svdR ./ max(svdR);
% svdG = svdG - min(svdG);
% svdG = svdG ./ max(svdG);
% svdB = svdB - min(svdB);
% svdB = svdB ./ max(svdB);
% Indexing
yTmp(brainMask) = svdR;
yRGB(:,:,1) = yTmp;
yTmp(brainMask) = svdG;
yRGB(:,:,2) = yTmp;
yTmp(brainMask) = svdB;
yRGB(:,:,3) = yTmp;
yTmp(brainMask) = svdY;
yRGB(:,:,4) = yTmp;

%% Mask out non-brain voxels
brainMaskRep = repmat(brainMask,[1 1 4]);
yRGB = yRGB .* double(brainMaskRep);

% SVD
for iSVD = 1:4,
    h = figure;
    imagesc(squeeze(yRGB(:,:,iSVD))); axis image; colorbar
    title(sprintf('%s S%02d C%d (%s) SV%d\n',...
        IOI.subj_name,s1,c1,colorNames{1+c1},iSVD),...
        'interpreter', 'none', 'FontSize', 14)
    set(h,'color','w')
    newName = sprintf('%s_S%02d_C%d_(%s)_SVD%02d.PNG',IOI.subj_name,s1,c1,colorNames{1+c1},iSVD);
    % Save as PNG
    print(h, '-dpng', '-r300', fullfile(pathName,newName));
    close(h)
end

% % Maximize figure
% set(0, 'Units', 'normalized');
% screenSize = get(0,'Screensize');
% screenSize = [0  0.0370 screenSize(3) screenSize(4)- 0.0370];

% % Display false RGB image
% h = figure; imshow(yRGB);
% title(sprintf('%s S%02d C%d (%s)\n',...
%     IOI.subj_name,s1,c1,colorNames{1+c1}),...
%     'interpreter', 'none', 'FontSize', 14)
% set(h,'color','w')
% % Normalized units
% set(h, 'Units', 'normalized');
% % Maximize figure
% set(h, 'OuterPosition', screenSize);

% newName = sprintf('%s_S%02d_C%d_(%s)',IOI.subj_name,s1,c1,colorNames{1+c1});
% Save as PNG (Maximize before printing)
% print(h, '-dpng', fullfile(pathName,newName), '-r300', 'WriteMode', 'overwrite');
% close(h)
fprintf('SVD done!\n\n')