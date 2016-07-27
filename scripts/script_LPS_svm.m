%% Load CO2 data
clear; clc;
NaClCO2 = [59.2; 49.4; 41.9; NaN; 62.6; NaN; 56.7; 50];
LPSCO2 = [NaN; 39.5; 54; 40.3; 80.3];
bilatROIsIdx = [(1:2:10)' (2:2:10)'];

%% Load HbR data 
% load('D:\Edgar\OIS_Results\networkResOut\results_S01_HbR.mat')
load('C:\Edgar\Dropbox\PostDoc\Newborn\OIS_Results\networkResOut\results_S01_HbR.mat');
% Extract Z for HbR
ZNaClHbR = results.Z(:,:,controlGroupIdx);
ZLPSHbR = results.Z(:,:,treatmentGroupIdx);
nNaCl = size(ZNaClHbR, 3);
nLPS = size(ZLPSHbR, 3);
ZLPSVecHbR=[];
ZNaClVecHbR=[];
[idxI, idxJ] =  ind2sub(size(squeeze(ZLPSHbR(:,:,1))),...
    find(~tril(ones(size(squeeze(ZLPSHbR(:,:,1)))))));
for iLPS = 1:nLPS,
%     for iROI = 1:size(bilatROIsIdx,1)
%         ZLPSVecHbR = [ZLPSVecHbR; ZLPSHbR(bilatROIsIdx(iROI, 1),bilatROIsIdx(iROI, 2), iLPS)];
%     end
ZLPSVecHbR = [ZLPSVecHbR; ZLPSHbR(idxI, idxJ, iLPS)];
end
for iNaCl = 1:nNaCl,
%     for iROI = 1:size(bilatROIsIdx,1)
%         ZNaClVecHbR = [ZNaClVecHbR; ZNaClHbR(bilatROIsIdx(iROI, 1),bilatROIsIdx(iROI, 2), iNaCl)];
%     end
ZNaClVecHbR = [ZNaClVecHbR; ZNaClHbR(idxI, idxJ, iNaCl)];
end

figure; set(gcf,'color','w')
hold on
% Plot HbR
plot(NaClCO2,reshape(ZNaClVecHbR, [numel(NaClCO2) numel(ZNaClVecHbR)/size(ZNaClHbR,3)]),...
    'bo','MarkerSize',12,'LineWidth',2)
plot(LPSCO2,reshape(ZLPSVecHbR, [numel(LPSCO2) numel(ZLPSVecHbR)/size(ZLPSHbR,3)]),...
    'bx','MarkerSize',12,'LineWidth',2)

%% Load HbO data
% load('D:\Edgar\OIS_Results\networkResOut\results_S01_HbO.mat')
load('C:\Edgar\Dropbox\PostDoc\Newborn\OIS_Results\networkResOut\results_S01_HbO.mat')
% Extract Z for HbO
ZNaCl = results.Z(:,:,controlGroupIdx);
ZLPS = results.Z(:,:,treatmentGroupIdx);

ZLPSVec=[];
ZNaClVec=[];
for iLPS = 1:nLPS,
%     for iROI = 1:size(bilatROIsIdx,1)
%         ZLPSVec = [ZLPSVec; ZLPS(bilatROIsIdx(iROI, 1),bilatROIsIdx(iROI, 2), iLPS)];
%     end
ZLPSVec = [ZLPSVec; ZLPS(idxI, idxJ, iLPS)];
end
for iNaCl = 1:nNaCl,
%     for iROI = 1:size(bilatROIsIdx,1)
%         ZNaClVec = [ZNaClVec; ZNaCl(bilatROIsIdx(iROI, 1),bilatROIsIdx(iROI, 2), iNaCl)];
%     end
ZNaClVec = [ZNaClVec; ZNaCl(idxI, idxJ, iLPS)];
end
hold on
plot(NaClCO2,reshape(ZNaClVec, [numel(NaClCO2) numel(ZNaClVec)/size(ZNaCl,3)]),...
    'ro','MarkerSize',12,'LineWidth',2)
plot(LPSCO2,reshape(ZLPSVec, [numel(LPSCO2) numel(ZLPSVec)/size(ZLPS,3)]),...
    'rx','MarkerSize',12,'LineWidth',2)
legend
ylabel('z(r)','FontSize',14)
xlabel('pCO_2 values','FontSize',14);
set(gca,'FontSize', 12)
% legend({'NaCl_{HbR}  '; 'LPS_{HbR}   '; 'NaCl_{HbO_2}'; 'LPS_{HbO_2} '},'Location','NorthWest')

%% SVM example - Training
clear; clc
load fisheriris
% xdata: matrix with 100 rows (samples) and 2 columns (observations)
xdata = meas(51:end,3:4);
% group: column vector (cell) with 100 rows with labels
group = species(51:end);
figure;
svmStruct = svmtrain(xdata,group,'ShowPlot',true, 'kernel_function', 'rbf',...
    'autoscale', true);

% SVM example - Classification
% Classify a new flower with petal length 5 and petal width 2, and circle the new point:
species = svmclassify(svmStruct,[5 2],'showplot',true)
hold on;plot(5,2,'ro','MarkerSize',12);hold off

%% SVM training with newborn data
% 1st predictor is CO2, x2-x6 are HbO bilateral correlations and x7-x11 are
% HbR bilteral correlations
% CO2, HbO & HbR Class loss: 40%
% xdata = [NaClCO2, reshape(ZNaClVec, [numel(NaClCO2) numel(ZNaClVec)/size(ZNaCl,3)]),...
%     reshape(ZNaClVecHbR, [numel(NaClCO2) numel(ZNaClVecHbR)/size(ZNaClHbR,3)]);...
%     LPSCO2, reshape(ZLPSVec, [numel(LPSCO2) numel(ZLPSVec)/size(ZLPS,3)]),...
%     reshape(ZLPSVecHbR, [numel(LPSCO2) numel(ZLPSVecHbR)/size(ZLPSHbR,3)])];
% HbO & HbR Class loss: 46%
% xdata = [ reshape(ZNaClVec, [numel(NaClCO2) numel(ZNaClVec)/size(ZNaCl,3)]),...
%     reshape(ZNaClVecHbR, [numel(NaClCO2) numel(ZNaClVecHbR)/size(ZNaClHbR,3)]);...
%      reshape(ZLPSVec, [numel(LPSCO2) numel(ZLPSVec)/size(ZLPS,3)]),...
%     reshape(ZLPSVecHbR, [numel(LPSCO2) numel(ZLPSVecHbR)/size(ZLPSHbR,3)])];
% HbR class loss: 53%
% xdata = [reshape(ZNaClVecHbR, [numel(NaClCO2) numel(ZNaClVecHbR)/size(ZNaClHbR,3)]);...
%     reshape(ZLPSVecHbR, [numel(LPSCO2) numel(ZLPSVecHbR)/size(ZLPSHbR,3)])];
% CO2 & all seed-to-seed correlations class loss: 23%

xdata = [NaClCO2, reshape(ZNaClVec, [numel(NaClCO2) numel(ZNaClVec)/size(ZNaCl,3)]),...
    reshape(ZNaClVecHbR, [numel(NaClCO2) numel(ZNaClVecHbR)/size(ZNaClHbR,3)]);...
    LPSCO2, reshape(ZLPSVec, [numel(LPSCO2) numel(ZLPSVec)/size(ZLPS,3)]),...
    reshape(ZLPSVecHbR, [numel(LPSCO2) numel(ZLPSVecHbR)/size(ZLPSHbR,3)])];
group = {};
xdata = xdata
% Remove columns containing NaNs
[~, NANc] = find(isnan(xdata));
xdata(:,NANc)=[];

% First 8 rows of xdata have NaCl samples
group(1:8,:) = {'NaCl'};
% Last 5 rows contain LPS samples
group(9:13,:) = {'LPS'};
svmStruct = svmtrain(xdata,group,'ShowPlot', false, 'kernel_function', 'rbf',...
    'autoscale', true);

%% SVM classification with newborn data
groupID = svmclassify(svmStruct,xdata(10,:),'showplot',false)

%% SVM cross validation
SVMModel = fitcsvm(zscore(xdata), group,'Standardize',true,'KernelFunction','RBF',...
    'KernelScale','auto');
% Cross validate the SVM classifier. By default, the software uses 10-fold
% cross validation.
CVSVMModel = crossval(SVMModel);
classLoss = kfoldLoss(CVSVMModel)

%% linear discriminant analysis confusion matrix
order = unique(group); % Order of the group labels
cp = cvpartition(group,'k',10); % Stratified cross-validation

f = @(xtr,ytr,xte,yte)confusionmat(yte,...
classify(xte,xtr,ytr),'order',order);

cfMat = crossval(f,zscore(xdata),group,'partition',cp);
disp(order')
cfMat = reshape(sum(cfMat), numel(order), numel(order))

% EOF
