%% Load CO2 data
clear; clc;
NaClCO2 = [59.2; 49.4; 41.9; NaN; 62.6; NaN; 56.7; 50];
LPSCO2 = [NaN; 39.5; 54; 40.3; 80.3];
onlyBilateral = false;
if onlyBilateral
    bilatROIsIdx = [(1:2:10)' (2:2:10)'];
end
% -------------------------------------------------------------------------
% Load HbR data 
load('D:\Edgar\OIS_Results\networkResOut\results_S01_HbR.mat')
% load('C:\Edgar\Dropbox\PostDoc\Newborn\OIS_Results\networkResOut\results_S01_HbR.mat');
% Extract Z for HbR
ZNaClHbR = results.Z(:,:,controlGroupIdx);
ZLPSHbR = results.Z(:,:,treatmentGroupIdx);
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
load('D:\Edgar\OIS_Results\networkResOut\results_S01_HbO.mat')
% load('C:\Edgar\Dropbox\PostDoc\Newborn\OIS_Results\networkResOut\results_S01_HbO.mat')
% Extract Z for HbO
ZNaCl = results.Z(:,:,controlGroupIdx);
ZLPS = results.Z(:,:,treatmentGroupIdx);

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
% legend({'NaCl_{HbR}  '; 'LPS_{HbR}   '; 'NaCl_{HbO_2}'; 'LPS_{HbO_2} '},'Location','NorthWest')

%% SVM example - Training
% clear; clc
% load fisheriris
% % xdata: matrix with 100 rows (samples) and 2 columns (observations)
% xdata = meas(51:end,3:4);
% % group: column vector (cell) with 100 rows with labels
% group = species(51:end);
% figure;
% svmStruct = svmtrain(xdata,group,'ShowPlot',true, 'kernel_function', 'rbf',...
%     'autoscale', true);
% 
% % SVM example - Classification
% % Classify a new flower with petal length 5 and petal width 2, and circle the new point:
% species = svmclassify(svmStruct,[5 2],'showplot',true)
% hold on;plot(5,2,'ro','MarkerSize',12);hold off

%% SVM training with newborn data
% 1st predictor is CO2, x2-x6 are HbO bilateral correlations and x7-x11 are
% HbR bilteral correlations
% CO2, HbO & HbR Class loss: 40%
% xdata = [NaClCO2, reshape(ZNaClVec, [nNaCl numel(ZNaClVec)/size(ZNaCl,3)]),...
%     reshape(ZNaClVecHbR, [nNaCl numel(ZNaClVecHbR)/size(ZNaClHbR,3)]);...
%     LPSCO2, reshape(ZLPSVec, [nLPS numel(ZLPSVec)/size(ZLPS,3)]),...
%     reshape(ZLPSVecHbR, [nLPS numel(ZLPSVecHbR)/size(ZLPSHbR,3)])];
% HbO & HbR Class loss: 46%
% xdata = [ reshape(ZNaClVec, [nNaCl numel(ZNaClVec)/size(ZNaCl,3)]),...
%     reshape(ZNaClVecHbR, [nNaCl numel(ZNaClVecHbR)/size(ZNaClHbR,3)]);...
%      reshape(ZLPSVec, [nLPS numel(ZLPSVec)/size(ZLPS,3)]),...
%     reshape(ZLPSVecHbR, [nLPS numel(ZLPSVecHbR)/size(ZLPSHbR,3)])];
% HbR class loss: 53%
% xdata = [reshape(ZNaClVecHbR, [nNaCl numel(ZNaClVecHbR)/size(ZNaClHbR,3)]);...
%     reshape(ZLPSVecHbR, [nLPS numel(ZLPSVecHbR)/size(ZLPSHbR,3)])];
% CO2 & all seed-to-seed correlations class loss: 23%

xdata = [NaClCO2, reshape(ZNaClVec, [nNaCl numel(ZNaClVec)/size(ZNaCl,3)]),...
    reshape(ZNaClVecHbR, [nNaCl numel(ZNaClVecHbR)/size(ZNaClHbR,3)]);...
    LPSCO2, reshape(ZLPSVec, [nLPS numel(ZLPSVec)/size(ZLPS,3)]),...
    reshape(ZLPSVecHbR, [nLPS numel(ZLPSVecHbR)/size(ZLPSHbR,3)])];
% clear
% load('D:\Edgar\OIS_Results\alff\alff_vals.mat')
% xdata = [HbONaCl HbRNaCl; HbOLPS HbRLPS];
% load('D:\Edgar\OIS_Results\lateralization\lateral_idx.mat')
% xdata = [ fcliNaCl.HbO' fcliNaCl.HbR'; fcliLPS.HbO' fcliLPS.HbR'];
% load('D:\Edgar\OIS_Results\averaged_maps\HbO\5\stats_R5_C5.mat')
% xdata = [reshape(NaCl,[11, 512*512]); reshape(LPS,[7, 512*512])];
% clear 'NaCl' 'LPS'

% load('D:\Edgar\OIS_Results\averaged_maps\stats_C5.mat')
% xdata = [NaCl_spatial_extension'; LPS_spatial_extension'];
% load('D:\Edgar\OIS_Results\averaged_maps\stats_C6.mat')
% xdata = [xdata [NaCl_spatial_extension'; LPS_spatial_extension']];


% load('D:\Edgar\OIS_Results\so2\avg_vals.mat')
% xdata = [HbONaCl HbRNaCl HbTNaCl sO2NaCl; HbOLPS HbRLPS HbTLPS sO2LPS];
% load('D:\Edgar\OIS_Results\networkResOut\resultsROI_Condition01_HbO.mat');

% load('D:\Edgar\OIS_Results\averaged_maps\HbO\10\spatial_extension_R10_C5.mat')
% xdata = [NaCl_spatial_extension; LPS_spatial_extension];
% xdata = [];
% load('D:\Edgar\OIS_Results\averaged_maps\HbR\10\spatial_extension_R10_C6.mat')
% xdata = [xdata [NaCl_spatial_extension; LPS_spatial_extension]];

% Remove columns containing NaNs
[~, NANc] = find(isnan(xdata));
xdata(:,NANc)=[];

group = {};
% First 8/11 rows of xdata have NaCl samples
group(1:8,:) = {'NaCl'};
% Last 5/7 rows contain LPS samples
group(9:13,:) = {'LPS'};
% svmStruct = svmtrain(xdata,group,'ShowPlot', false, 'kernel_function', 'rbf',...
%     'autoscale', true);

%% SVM classification with newborn data
% groupID = svmclassify(svmStruct,xdata(4,:),'showplot',false)

%% k-fold Cross validation using SVM classifier
addpath(genpath('D:\Edgar\biolearning'))
clc
k = 10;                                     % Number of folds

% CVO = cvpartition(group,'k',10);
% boxconstraint — One strategy is to try a geometric sequence of the box
% constraint parameter. For example, take 11 values, from 1e-5 to 1e5 by a
% factor of 10.
% 
% 
% rbf_sigma — One strategy is to try a geometric sequence of the RBF sigma
% parameter. For example, take 11 values, from 1e-5 to 1e5 by a factor of
% 10.
SigmaVals = logspace(-5, 5, 11);
BoxVals = logspace(-5, 5, 11);
for iSigma = 8,
    for iBox = 4,
        cvFolds = crossvalind('Kfold', group, k);   %# get indices of 10-fold CV
        cp = classperf(group);                      %# init performance tracker
        
        for i = 1:k                                 %# for each fold
            testIdx = (cvFolds == i);               %# get indices of test instances
            trainIdx = ~testIdx;                    %# get indices training instances
            
            xtrain = xdata(trainIdx,:);
            ytrain = group(trainIdx);
            xtest = xdata(testIdx,:);
            
            %     minfn = @(z)crossval('mcr',xdata,group,'Predfun', ...
            %     @(xtrain,ytrain,xtest)crossfun(xtrain,ytrain,...
            %     xtest,exp(z(1)),exp(z(2))),'partition',CVO);
            
            % opts = optimset('TolX',5e-4,'TolFun',5e-4);
            % [searchmin, fval] = fminsearch(minfn,randn(2,1),opts);
            
            % z = exp(searchmin);
            
            % rbf_sigma = z(1); %0.8
            rbf_sigma = SigmaVals(iSigma);
            % boxconstraint = z(2); % 2e-1
            boxconstraint = BoxVals(iBox);
            %# train an SVM model over training instances
            svmModel = svmtrain(xtrain, ytrain, ...
                'Autoscale',true, 'Showplot',false, 'Method','SMO', ...
                'BoxConstraint',boxconstraint, 'Kernel_Function','rbf', ...
                'RBF_Sigma',rbf_sigma, 'tolkkt', 1e-6);
            
            %# test using test instances
            pred{i} = svmclassify(svmModel, xtest, 'Showplot',false);
            
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

%% SVM cross validation using fitcsvm
% SVMModel = fitcsvm(zscore(xdata), group,'Standardize',true,'KernelFunction','RBF',...
%     'KernelScale','auto');
% % Cross validate the SVM classifier. By default, the software uses 10-fold
% % cross validation.
% CVSVMModel = crossval(SVMModel);
% classLoss = kfoldLoss(CVSVMModel)

%% linear discriminant analysis confusion matrix
% order = unique(group); % Order of the group labels
% cp = cvpartition(group,'k',10); % Stratified cross-validation
% 
% f = @(xtr,ytr,xte,yte)confusionmat(yte,...
% classify(xte,xtr,ytr),'order',order);
% 
% cfMat = crossval(f,zscore(xdata),group,'partition',cp);
% disp(order')
% cfMat = reshape(sum(cfMat), numel(order), numel(order))

% EOF
