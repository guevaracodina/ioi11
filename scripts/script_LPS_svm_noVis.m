%% Load CO2 data
clear; clc;
resultsFolder = 'C:\Edgar\Dropbox\PostDoc\Newborn\OIS_Results\svm';
NaClCO2 = [59.2; 49.4; 41.9; NaN; 62.6; NaN; 56.7; 50];
LPSCO2 = [NaN; 39.5; 54; 40.3; 80.3];
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

figure; set(gcf,'color','w')
hold on
% Plot HbR
plot(NaClCO2,reshape(ZNaClVecHbR, [nNaCl numel(ZNaClVecHbR)/size(ZNaClHbR,3)]),...
    'bo','MarkerSize',12,'LineWidth',2)
plot(LPSCO2,reshape(ZLPSVecHbR, [nLPS numel(ZLPSVecHbR)/size(ZLPSHbR,3)]),...
    'bx','MarkerSize',12,'LineWidth',2)

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
hold on
plot(NaClCO2,reshape(ZNaClVec, [nNaCl numel(ZNaClVec)/size(ZNaCl,3)]),...
    'ro','MarkerSize',12,'LineWidth',2)
plot(LPSCO2,reshape(ZLPSVec, [nLPS numel(ZLPSVec)/size(ZLPS,3)]),...
    'rx','MarkerSize',12,'LineWidth',2)
legend
ylabel('z(r)','FontSize',14)
xlabel('pCO_2 values','FontSize',14);
set(gca,'FontSize', 12)

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

%% k-fold Cross validation using SVM classifier
% addpath(genpath('D:\Edgar\biolearning'))
clc
k = 10;                                     % Number of folds

% CVO = cvpartition(group,'k',10);
% boxconstraint — One strategy is to try a geometric sequence of the box
% constraint parameter. For example, take 11 values, from 1e-5 to 1e5 by a
% factor of 10.
% 
% rbf_sigma — One strategy is to try a geometric sequence of the RBF sigma
% parameter. For example, take 11 values, from 1e-5 to 1e5 by a factor of
% 10.
SigmaVals = logspace(-5, 5, 11);
BoxVals = logspace(-5, 5, 11);
for iSigma = 1:11,
    for iBox = 1:11,
        cvFolds = crossvalind('Kfold', group, k);   %# get indices of 10-fold CV
        cp = classperf(group);                      %# init performance tracker
        
        for i = 1:k                                 %# for each fold
            testIdx = (cvFolds == i);               %# get indices of test instances
            trainIdx = ~testIdx;                    %# get indices training instances
            
            xtrain = xdata(trainIdx,:);
            ytrain = group(trainIdx);
            xtest = xdata(testIdx,:);
            ytest = group(testIdx);
            
            rbf_sigma = SigmaVals(iSigma);
            % boxconstraint = z(2); % 2e-1
            boxconstraint = BoxVals(iBox);
            %# train an SVM model over training instances
            svmModel = svmtrain(xtrain, ytrain, ...
                'Autoscale',true, 'Showplot',false, 'Method','SMO', ...
                'BoxConstraint',boxconstraint, 'Kernel_Function','rbf', ...
                'RBF_Sigma',rbf_sigma, 'tolkkt', 1e-6);
            
            %# test using test instances, save output labels to create confusion matrix
            pred{i} = svmclassify(svmModel, xtest, 'Showplot',false);
            % Keep target labels to create confusion matrix
            groundTruth{i} = ytest;
            
            %# evaluate and update performance object
            cp = classperf(cp, pred{i}, testIdx);
        end
        
        %# get accuracy
        fprintf('Accuracy  = %0.2f %%\n', 100*cp.CorrectRate);
        accMat(iSigma, iBox) = cp.CorrectRate;
        % disp(cp.ClassLabels');
        %# get confusion matrix
        %# columns:actual, rows:predicted, last-row: unclassified instances
        % disp(cp.CountingMatrix);
        fprintf('\t\t%s\t%s\n %s\t\t%d\t%d\n %s\t\t%d\t%d\n Unclass\t%d\t%d\n',cp.ClassLabels{1}, cp.ClassLabels{2}, cp.ClassLabels{1},...
            cp.CountingMatrix(1,:), cp.ClassLabels{2}, cp.CountingMatrix(2,:), cp.CountingMatrix(3,:));
    end
end

%% Create confusion plot
outputsCell = {};
targetsCell = {};
clear targets outputs
for iSample = 1:k,
    if ~isempty(pred{iSample})
        outputsCell = [outputsCell; pred{iSample}];
    end
    if ~isempty(groundTruth{iSample})
        targetsCell = [targetsCell; groundTruth{iSample}];
    end
end
for iSample = 1:nLPS + nNaCl
    % LPS are class 1 and NaCl are class 0
    targets(iSample) = strcmp(targetsCell{iSample},'LPS');
    outputs(iSample) = strcmp(outputsCell{iSample},'LPS');
end
% cleanup
clear i iBox iLPS iNaCl iSigma idx results ss tmpVec
% Save results
% save(fullfile(resultsFolder,'LPS_svm.mat'))

%% Plot confusion matrix
hConf = figure; set(hConf, 'color', 'w')
plotconfusion(targets, outputs)
set(gca, 'YTickLabel', {cp.ClassLabels{1} cp.ClassLabels{2} ' ' },...
    'XTickLabel', {cp.ClassLabels{1} cp.ClassLabels{2} ' ' },...
    'FontWeight','bold',...
	'FontSize',16);
% Specify window units
set(hConf, 'units', 'inches')
% Change figure and paper size
set(hConf, 'Position', [0.1 0.1 4 4])
set(hConf, 'PaperPosition', [0.1 0.1 4 4])
% Save as PNG at the user-defined resolution
% print(hConf, '-dpng', ...
%     fullfile(resultsFolder, 'LPS_svm_confusionplot.png'),...
%     sprintf('-r%d',1200));
hROC = figure; set(hROC, 'color', 'w')
plotroc(targets, outputs)

% EOF
