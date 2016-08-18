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

%% Histology data
% LP02, LP03, LP04 & LP05
LPShist.L = [2002507; 2149179; 957905; 6718538];
LPShist.R = [168662; 1142728; 909365; 6718538];
LPShist.T = LPShist.L + LPShist.R;
% NC07, NC08 & NC09
NaClhist.L = [122725; 135333; 194805];
NaClhist.R = [110830; 159154; 198182];
NaClhist.T = NaClhist.L + NaClhist.R;
% Indices of subjects to keep
idx2keep = [6:8, 10:13]';

%% Crop X-data
groupLabels = {'NC01'; 'NC02'; 'NC03'; 'NC05'; 'NC06'; 'NC07'; 'NC08'; 'NC09';...
    'LP01'; 'LP02'; 'LP03'; 'LP04'; 'LP05'};
groupLabels = groupLabels(idx2keep);
xdata = xdata(idx2keep,:);
ydata = [NaClhist.T; LPShist.T];
group = group(idx2keep);

%% k-fold Cross validation using SVM classifier
clc
k = 5;                                     % Number of folds
cvFolds = crossvalind('Kfold', group, k);   %# get indices of k-fold CV
% Initialize performance trackers
pred = [];
groundTruth = [];
targetLabels = {};
for i = 1:k                                 %# for each fold
    testIdx = (cvFolds == i);               %# get indices of test instances
    trainIdx = ~testIdx;                    %# get indices training instances
    
    xtrain = xdata(trainIdx,:);
    ytrain = ydata(trainIdx,:);
    
    xtest = xdata(testIdx,:);
    ytest = ydata(testIdx,:);
    
    % Bayesian regularization (takes a Q x 56 matrix, where Q = No. samples)
    pred = [pred; myNeuralNetworkFunction(xtest')'];

    % Keep target values to create correlation plot
    groundTruth = [groundTruth; ytest];
    
    % Keep target labels 
    targetLabels = [targetLabels; groupLabels(testIdx)];
end
        
%% Plot correlation between measured and predicted values
p = polyfit(groundTruth, pred, 1);
yfit = p(1)*groundTruth + p(2);
r = corr(groundTruth, pred);
RMSEP = sqrt(sum((groundTruth-pred).^2)/size(xdata,1));

% Find indices of LPS
idxLPS = find(~cellfun(@isempty, regexp(targetLabels,'LP')));
idxNaCl = find(~cellfun(@isempty, regexp(targetLabels,'NC')));

figure; hold on;
% Plot Predicted = Measured
plot(groundTruth, groundTruth, 'k:', 'LineWidth', 3)
% Plot fit
plot(groundTruth, yfit, 'b-', 'LineWidth', 3)
% Plot NaCl
plot(groundTruth(idxNaCl), pred(idxNaCl), 'ko', 'LineWidth', 3, 'MarkerSize', 12)
% Plot LPS
plot(groundTruth(idxLPS), pred(idxLPS), 'rx', 'LineWidth', 3, 'MarkerSize', 12)

axis equal
legend({'Predicted = Measured' 'Linear fit' 'NaCl' 'LPS' }, 'Location', 'SouthEast')
set(gca, 'FontSize', 14)
xlim([min(groundTruth) max(groundTruth)])
ylim([min(groundTruth) max(groundTruth)])
xlabel('Ventricular Lesion Size (pixels)', 'FontSize', 14); 
ylabel('Predicted Ventricular Size (pixels)', 'FontSize', 14)
title(sprintf('r = %0.4f', r), 'FontSize', 14); 

% EOF