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
LPS = [1363882; 2332350; 1121757; 11627024];
% NC07, NC08 & NC09
NaCl = [273222; 258777; 272253];

%% Crop X-data
groupLabels = {'NC01'; 'NC02'; 'NC03'; 'NC05'; 'NC06'; 'NC07'; 'NC08'; 'NC09';...
    'LP01'; 'LP02'; 'LP03'; 'LP04'; 'LP05'};
idx2keep = [6:8, 10:13]';
% idx2keep = [1:13]';
xdata = xdata(idx2keep,:);
ydata = [NaCl; LPS];
group = group(idx2keep);

%% Use svr_trainer
addpath('C:\Edgar\Dropbox\Matlab')
svrobj = svr_trainer(xdata(1:6,:),ydata(1:6,:),400,0.000000025,'gaussian',0.5);
y = svrobj.predict(xdata(7,:));

%% k-fold Cross validation using SVM classifier
% addpath(genpath('D:\Edgar\biolearning'))
clc
k = 5;                                     % Number of folds

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
ioi_text_waitbar(0, 'Please wait...');
for iSigma = 1:numel(SigmaVals),
    for iBox = 1:11,
        cvFolds = crossvalind('Kfold', group, k);   %# get indices of k-fold CV
%         cp = classperf(ydata);                      %# init performance tracker
        pred = [];
        groundTruth = [];
        for i = 1:k                                 %# for each fold
            testIdx = (cvFolds == i);               %# get indices of test instances
            trainIdx = ~testIdx;                    %# get indices training instances
            
            xtrain = xdata(trainIdx,:);
%             ytrain = group(trainIdx);
            ytrain = ydata(trainIdx,:);
            xtest = xdata(testIdx,:);
%             ytest = group(testIdx);
            ytest = ydata(testIdx,:);
            
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
%             svmModel = svmtrain(xtrain, ytrain, ...
%                 'Autoscale',true, 'Showplot',false, 'Method','SMO', ...
%                 'BoxConstraint',boxconstraint, 'Kernel_Function','rbf', ...
%                 'RBF_Sigma',rbf_sigma, 'tolkkt', 1e-6);
            svmModel = svr_trainer(xtrain,ytrain,boxconstraint,1e-6,'gaussian',rbf_sigma);
            %# test using test instances, save output labels to create confusion matrix
%             pred{i} = svmclassify(svmModel, xtest, 'Showplot',false);
            pred  = [pred; svrobj.predict(xtest)];
            % Keep target labels to create confusion matrix
            groundTruth = [groundTruth; ytest];
            
            %# evaluate and update performance object
%             cp = classperf(cp, pred{i}, testIdx);
            close all
        end
        corrIdx(iSigma, iBox) = corr(groundTruth, pred);
        %# get accuracy
%         fprintf('Accuracy  = %0.2f %%\n', 100*cp.CorrectRate);
%         accMat(iSigma, iBox) = cp.CorrectRate;
        % disp(cp.ClassLabels');
        %# get confusion matrix
        %# columns:actual, rows:predicted, last-row: unclassified instances
        % disp(cp.CountingMatrix);
%         fprintf('\t\t%s\t%s\n %s\t\t%d\t%d\n %s\t\t%d\t%d\n Unclass\t%d\t%d\n',cp.ClassLabels{1}, cp.ClassLabels{2}, cp.ClassLabels{1},...
%             cp.CountingMatrix(1,:), cp.ClassLabels{2}, cp.CountingMatrix(2,:), cp.CountingMatrix(3,:));
    end
    ioi_text_waitbar(iSigma/numel(SigmaVals), sprintf('Processing event %d from %d', ...
        iSigma, numel(SigmaVals)));
end
ioi_text_waitbar('Clear');

%% 
pred = -(pred - max(pred));
pred = pred .* max(groundTruth)/max(pred);
p = polyfit(groundTruth, pred, 1);
yfit = p(1)*groundTruth + p(2);
figure; hold on;
plot(groundTruth([1 4 6]), pred([1 4 6]), 'ko')
plot(groundTruth([2 3 5 7]), pred([2 3 5 7]), 'rx')
plot(groundTruth, yfit, 'b-')
legend({'NaCl' 'LPS' 'Linear fit'}, 'Location', 'SouthEast')
% xlim([0 12e6])
% ylim([0 12e6])
xlabel('Measured Lesion size (pixels)'); ylabel('Predicted lesion size (pixels)')
r = corr(groundTruth, pred);
RMSEP = sqrt(sum((groundTruth-pred).^2)/size(xdata,1));
% EOF