%% Load connectivity data
resultsFolder = 'C:\Edgar\Dropbox\PostDoc\Newborn\OIS_Results\svm';
% NaClCO2 = [59.2; 49.4; 41.9; NaN; 62.6; NaN; 56.7; 50];
% LPSCO2 = [NaN; 39.5; 54; 40.3; 80.3];
onlyBilateral = false;
if onlyBilateral
    bilatROIsIdx = [(1:2:10)' (2:2:10)'];
end
% -------------------------------------------------------------------------
% Load HbR data 
% load('D:\Edgar\OIS_Results\networkResOut\results_S01_HbR.mat')
load('C:\Edgar\Dropbox\PostDoc\Newborn\OIS_Results\networkResOut\results_S01_HbR.mat');
% Extract Z for HbR
ZNaClHbR = results.Z(3:end,3:end,controlGroupIdx);
ZLPSHbR = results.Z(3:end,3:end,treatmentGroupIdx);
nNaCl = size(ZNaClHbR, 3);
nLPS = size(ZLPSHbR, 3);
ZLPSVecHbR=[];
ZNaClVecHbR=[];
% Indices of all ROI-to-ROI fc values
% [idxI, idxJ] =  ind2sub(size(squeeze(ZLPSHbR(:,:,1))),...
%     find(~tril(ones(size(squeeze(ZLPSHbR(:,:,1)))))));
idx = find(~tril(ones(size(squeeze(ZLPSHbR(:,:,1))))));
for iLPS = 1:nLPS,
    if onlyBilateral
        for iROI = 1:size(bilatROIsIdx,1)
            ZLPSVecHbR = [ZLPSVecHbR; ZLPSHbR(bilatROIsIdx(iROI, 1),bilatROIsIdx(iROI, 2), iLPS)];
        end
    else
        tmpVec = squeeze(ZLPSHbR(:,:,iLPS));
        ZLPSVecHbR = [ZLPSVecHbR; tmpVec(idx)];
%         ZLPSVecHbR = [ZLPSVecHbR; reshape(ZLPSHbR(idxI, idxJ, iLPS),[numel(idxI) 1])];
    end
end
for iNaCl = 1:nNaCl,
    if onlyBilateral
        for iROI = 1:size(bilatROIsIdx,1)
            ZNaClVecHbR = [ZNaClVecHbR; ZNaClHbR(bilatROIsIdx(iROI, 1),bilatROIsIdx(iROI, 2), iNaCl)];
        end
    else
        tmpVec = squeeze(ZNaClHbR(:,:,iNaCl));
        ZNaClVecHbR = [ZNaClVecHbR; tmpVec(idx)];
    end
end

% -------------------------------------------------------------------------
% Load HbO data
% load('D:\Edgar\OIS_Results\networkResOut\results_S01_HbO.mat')
load('C:\Edgar\Dropbox\PostDoc\Newborn\OIS_Results\networkResOut\results_S01_HbO.mat')
% Extract Z for HbO
ZNaCl = results.Z(3:end,3:end,controlGroupIdx);
ZLPS = results.Z(3:end,3:end,treatmentGroupIdx);

ZLPSVec=[];
ZNaClVec=[];
for iLPS = 1:nLPS,
    if onlyBilateral
        for iROI = 1:size(bilatROIsIdx,1)
            ZLPSVec = [ZLPSVec; ZLPS(bilatROIsIdx(iROI, 1),bilatROIsIdx(iROI, 2), iLPS)];
        end
    else
        tmpVec = squeeze(ZLPS(:,:,iLPS));
        ZLPSVec = [ZLPSVec; tmpVec(idx)];
    end
end
for iNaCl = 1:nNaCl,
    if onlyBilateral
        for iROI = 1:size(bilatROIsIdx,1)
            ZNaClVec = [ZNaClVec; ZNaCl(bilatROIsIdx(iROI, 1),bilatROIsIdx(iROI, 2), iNaCl)];
        end
    else
        tmpVec = squeeze(ZNaCl(:,:,iNaCl));
        ZNaClVec = [ZNaClVec; tmpVec(idx)];
    end
end

%% SVM training with newborn data

% xdata = [NaClCO2, reshape(ZNaClVec, [nNaCl numel(ZNaClVec)/size(ZNaCl,3)]),...
%     reshape(ZNaClVecHbR, [nNaCl numel(ZNaClVecHbR)/size(ZNaClHbR,3)]);...
%     LPSCO2, reshape(ZLPSVec, [nLPS numel(ZLPSVec)/size(ZLPS,3)]),...
%     reshape(ZLPSVecHbR, [nLPS numel(ZLPSVecHbR)/size(ZLPSHbR,3)])];

xdata = [reshape(ZNaClVec, [nNaCl numel(ZNaClVec)/size(ZNaCl,3)]),...
    reshape(ZNaClVecHbR, [nNaCl numel(ZNaClVecHbR)/size(ZNaClHbR,3)]);...
    reshape(ZLPSVec, [nLPS numel(ZLPSVec)/size(ZLPS,3)]),...
    reshape(ZLPSVecHbR, [nLPS numel(ZLPSVecHbR)/size(ZLPSHbR,3)])];


% Remove columns containing NaNs
[~, NANc] = find(isnan(xdata));
xdata(:,NANc)=[];

group = {};
% First 8/11 rows of xdata have NaCl samples
group(1:8,:) = {'NaCl'};
% Last 5/7 rows contain LPS samples
group(9:13,:) = {'LPS'};

%% Histology data
% LP03, LP04, LP05 & LP02
LPS = [2332350; 1121757; 11627024; 1363882];
% NC08, NC09 & NC07
NaCl = [258777; 272253; 273222];
