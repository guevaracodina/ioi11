%% Load connectivity data
resultsFolder = 'C:\Edgar\Dropbox\PostDoc\Newborn\OIS_Results\ANN';
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
load('C:\Edgar\Dropbox\PostDoc\Newborn\OIS_Results\networkResOutNoVis\resultsROI_Condition01_HbR.mat')
% Extract Z for HbR
% ZNaClHbR = results.Z(3:end,3:end,controlGroupIdx);
% ZLPSHbR = results.Z(3:end,3:end,treatmentGroupIdx);
ZNaClHbR = r(3:end,3:end,controlGroupIdx);
ZLPSHbR = r(3:end,3:end,treatmentGroupIdx);

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
load('C:\Edgar\Dropbox\PostDoc\Newborn\OIS_Results\networkResOutNoVis\resultsROI_Condition01_HbO.mat')
% Extract Z for HbO
% ZNaCl = results.Z(3:end,3:end,controlGroupIdx);
% ZLPS = results.Z(3:end,3:end,treatmentGroupIdx);
ZNaCl = r(3:end,3:end,controlGroupIdx);
ZLPS = r(3:end,3:end,treatmentGroupIdx);

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

%% SVM data preparation with newborn data

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

%% Find correlation (HbR)
clear x y;
i = 1;
for iRow=1:size(ZLPSHbR,1)
    for iCol=1:size(ZLPSHbR,2)
        for iLPS=1:size(ZLPSHbR,3)
            xHbR(i,:) = squeeze(ZLPSHbR(iRow, iCol, iLPS));
            yHbR(i,:) = squeeze(ZNaClHbR(iRow, iCol, iLPS));
            x(i,:) = squeeze(ZLPS(iRow, iCol, iLPS));
            y(i,:) = squeeze(ZNaCl(iRow, iCol, iLPS));
            i = i + 1;
        end
    end
end
figure; subplot(121)
plot(xHbR,yHbR,'k.')
subplot(122)
plot(x,y,'r.')
[rhoHbR, pHbR] = corr(xHbR,yHbR,'rows','complete');
[rho, p] = corr(x,y,'rows','complete');

%% Compute average of the correlation matrices
addpath(genpath('C:\Edgar\Dropbox\Matlab\'))
alphaVal = 0.05;
LPSavg = mean(ZLPS,3);
NaClavg = mean(ZNaCl,3);
LPSHbRavg = mean(ZLPSHbR,3);
NaClHbRavg = mean(ZNaClHbR,3);

% Stat test
NaNmask = tril(ones(size(squeeze(ZNaCl(:,:,1)))));
NaNmask = repmat(NaNmask, [1 1 size(ZNaCl, 3)]);
idx = find(NaNmask);
% NaCl (HbO)
ZNaCl(idx) = NaN;
[tMap, pMapFDR, pMapFDRalpha] = ioi_stat_map(ZNaCl, alphaVal);
NaClavg(isnan(pMapFDRalpha)) = NaN;

% NaCl (HbR)
ZNaClHbR(idx) = NaN;
[tMap, pMapFDR, pMapFDRalpha] = ioi_stat_map(ZNaClHbR, alphaVal);
NaClHbRavg(isnan(pMapFDRalpha)) = NaN;

NaNmask = tril(ones(size(squeeze(ZLPS(:,:,1)))));
NaNmask = repmat(NaNmask, [1 1 size(ZLPS, 3)]);
idx = find(NaNmask);

% LPS (HbO)
ZLPS(idx) = NaN;
[tMap, pMapFDR, pMapFDRalpha] = ioi_stat_map(ZLPS, alphaVal);
LPSavg(isnan(pMapFDRalpha)) = NaN;

% LPS (HbR)
ZLPSHbR(idx) = NaN;
[tMap, pMapFDR, pMapFDRalpha] = ioi_stat_map(ZLPSHbR, alphaVal);
LPSHbRavg(isnan(pMapFDRalpha)) = NaN;

% circular graphs
% h1 = schemaball(LPSavg, names(3:10), [0 0 1; 1 0 0], [1 1 1]);
% h2 = circularGraph(LPSavg,'Colormap',ioi_get_colormap('bipolar', size(LPSavg, 1)),'Label',names(3:10));

%% Circos preparation for LPS (HbO)
clear circosLPS
% Rescale [-1, 1] to [0, 2000]
offset = 1000;
scale = 1000;
LPSavgScaled = fix(scale * LPSavg + offset);
% Add '-' for missing values
idx = find((~triu(ones(size(LPSavgScaled) + [3 3]))));
circosLPS = cell(size(LPSavgScaled) + [3 3]);
circosLPS(4:end, 4:end) = num2cell(LPSavgScaled);
circosLPS(idx) = {'-'};
circosLPS(1:3, 1:3) = {'-'};

% Specify column order
colOrder = {3 8 4 7 5 6 2 1};
circosLPS(1, 4:end) = colOrder;
circosLPS(4:end, 1) = colOrder';

% Set segments color
segColor = {'0,148,255'	'0,148,255'	'255,106,0'	'255,106,0'	...
    '178,0,255'	'178,0,255'	'76,255,0'	'76,255,0'};
circosLPS(2, 4:end) = segColor;
circosLPS(4:end, 2) = segColor';

% Add labels to segments
circosLPS(3, 4:end) = names(3:10);
circosLPS(4:end, 3) = names(3:10)';      

%% Circos preparation for LPS (HbR)
clear circosHbRLPS
% Rescale [-1, 1] to [0, 2000]
offset = 1000;
scale = 1000;
LPSHbRavgScaled = fix(scale * LPSHbRavg + offset);
% Add '-' for missing values
idx = find((~triu(ones(size(LPSHbRavgScaled) + [3 3]))));
circosHbRLPS = cell(size(LPSHbRavgScaled) + [3 3]);
circosHbRLPS(4:end, 4:end) = num2cell(LPSHbRavgScaled);
circosHbRLPS(idx) = {'-'};
circosHbRLPS(1:3, 1:3) = {'-'};

%  Specify column order
circosHbRLPS(1, 4:end) = colOrder;
circosHbRLPS(4:end, 1) = colOrder';

% Set segments color
circosHbRLPS(2, 4:end) = segColor;
circosHbRLPS(4:end, 2) = segColor';

% Add labels to chromosomes
circosHbRLPS(3, 4:end) = names(3:10);
circosHbRLPS(4:end, 3) = names(3:10)';

%% Circos preparation for NaCl (HbR)
clear circosHbRNaCl
% Rescale [-1, 1] to [0, 2000]
offset = 1000;
scale = 1000;
NaClHbRavgScaled = fix(scale * NaClHbRavg + offset);
% Add '-' for missing values
idx = find((~triu(ones(size(NaClHbRavgScaled) + [3 3]))));
circosHbRNaCl = cell(size(NaClHbRavgScaled) + [3 3]);
circosHbRNaCl(4:end, 4:end) = num2cell(NaClHbRavgScaled);
circosHbRNaCl(idx) = {'-'};
circosHbRNaCl(1:3, 1:3) = {'-'};

% Specify column order
circosHbRNaCl(1, 4:end) = colOrder;
circosHbRNaCl(4:end, 1) = colOrder';

% Set segments color
circosHbRNaCl(2, 4:end) = segColor;
circosHbRNaCl(4:end, 2) = segColor';

% Add labels to segments
circosHbRNaCl(3, 4:end) = names(3:10);
circosHbRNaCl(4:end, 3) = names(3:10)';

%% Circos preparation for NaCl (HbO)
clear circosNaCl
% Rescale [-1, 1] to [0, 2000]
offset = 1000;
scale = 1000;
NaClavgScaled = fix(scale * NaClavg + offset);
% Add '-' for missing values
idx = find((~triu(ones(size(NaClavgScaled) + [3 3]))));
circosNaCl = cell(size(NaClavgScaled) + [3 3]);
circosNaCl(4:end, 4:end) = num2cell(NaClavgScaled);
circosNaCl(idx) = {'-'};
circosNaCl(1:3, 1:3) = {'-'};

% Specify column order
circosNaCl(1, 4:end) = colOrder;
circosNaCl(4:end, 1) = colOrder';

% Set segments color
circosNaCl(2, 4:end) = segColor;
circosNaCl(4:end, 2) = segColor';

% Add labels to segments
circosNaCl(3, 4:end) = names(3:10);
circosNaCl(4:end, 3) = names(3:10)';

%% Create colormap (12 colors)
fix(linspace(0, 2000, 12 + 1))
fix(255*ioi_get_colormap('bipolar',12))
% EOF