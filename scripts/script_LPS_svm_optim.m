%% script_LPS_svm_optim
clear; close all; clc;

%% 1. Generate the 10 base points for each class:
grnpop = mvnrnd([1,0],eye(2),10);
redpop = mvnrnd([0,1],eye(2),10);

%% 2. View the base points:
plot(grnpop(:,1),grnpop(:,2),'go')
hold on
plot(redpop(:,1),redpop(:,2),'ro')
hold off

%% 3. Generate the 100 data points of each class:
redpts = zeros(100,2);grnpts = redpts;
for i = 1:100
    grnpts(i,:) = mvnrnd(grnpop(randi(10),:),eye(2)*0.2);
    redpts(i,:) = mvnrnd(redpop(randi(10),:),eye(2)*0.2);
end

%% 4. View the data points:
figure
plot(grnpts(:,1),grnpts(:,2),'go')
hold on
plot(redpts(:,1),redpts(:,2),'ro')
hold off

%% 5. Put the data into one matrix, and make a vector grp that labels the class of each point:
cdata = [grnpts;redpts];
grp = ones(200,1);
% green label 1, red label -1
grp(101:200) = -1;

%% 6. Check the basic classification of all the data using the default parameters:
svmStruct = svmtrain(cdata,grp,'Kernel_Function','rbf',...
'showplot',true);

%% 7. Write a function called crossfun to calculate the predicted classification 
% yfit from a test vector xtest, when the SVM is trained on a sample xtrain
% that has classification ytrain. Since you want to find the best
% parameters rbf_sigma and boxconstraint, include those in the function.

%% 8. Set up a partition for cross validation. This step causes the cross validation to be fixed. 
% Without this step, the cross validation is random, so a minimization
% procedure can find a spurious local minimum.
c = cvpartition(200,'kfold',10);

%% 9. Set up a function that takes an input z=[rbf_sigma,boxconstraint], 
% and returns the cross-validation value of exp(z). The reason to take
% exp(z) is twofold: rbf_sigma and boxconstraint must be positive. You
% should look at points spaced approximately exponentially apart. This
% function handle computes the cross validation at parameters
% exp([rbf_sigma,boxconstraint]):
minfn = @(z)crossval('mcr',cdata,grp,'Predfun', ...
    @(xtrain,ytrain,xtest)crossfun(xtrain,ytrain,...
    xtest,exp(z(1)),exp(z(2))),'partition',c);

%% 10. Search for the best parameters [rbf_sigma,boxconstraint] with fminsearch, 
% setting looser tolerances than the defaults. Tip   If you have a Global
% Optimization Toolbox license, use patternsearch for faster, more reliable
% minimization. Give bounds on the components of z to keep the optimization
% in a sensible region, such as [–5,5], and give a relatively loose TolMesh
% tolerance.
 
opts = optimset('TolX',5e-4,'TolFun',5e-4);
[searchmin fval] = fminsearch(minfn,randn(2,1),opts)

% searchmin =
%     0.9758
%    -0.1569
% 
% fval =
%     0.3350

% The best parameters [rbf_sigma;boxconstraint] in this run are:
z = exp(searchmin)
% z =
%     2.6534
%     0.8548

%% 11. Since the result of fminsearch can be a local minimum, not a global minimum, 
% try again with a different starting point to check that your result is meaningful:
[searchmin fval] = fminsearch(minfn,randn(2,1),opts)

% searchmin =
%     0.2778
%     0.6395
% 
% fval =
%     0.3100
% 
% The best parameters [rbf_sigma;boxconstraint] in this run are:
z = exp(searchmin)
% z =
%     1.3202
%     1.8956

%% 12. Try another search:
[searchmin fval] = fminsearch(minfn,randn(2,1),opts)

% searchmin =
%    -0.0749
%     0.6085
% 
% fval =
%     0.2850

% The third search obtains the lowest function value. The final parameters are:
z = exp(searchmin)
% z =
%     0.9278
%     1.8376

% The default parameters [1,1] are close to optimal for this data and partition.

%% 13. Use the z parameters to train a new SVM classifier:
svmStruct = svmtrain(cdata,grp,'Kernel_Function','rbf',...
'rbf_sigma',z(1),'boxconstraint',z(2),'showplot',true);

%% 14. Generate and classify some new data points:
grnobj = gmdistribution(grnpop,.2*eye(2));
redobj = gmdistribution(redpop,.2*eye(2));

newData = random(grnobj,10);
newData = [newData;random(redobj,10)];
grpData = ones(20,1);
grpData(11:20) = -1; % red = -1

v = svmclassify(svmStruct,newData,'showplot',true);

%% 15. See which new data points are correctly classified. 
% Circle the correctly classified points in red, and the incorrectly
% classified points in black.
mydiff = (v == grpData); % classified correctly
hold on
for ii = mydiff % plot red circles around correct pts
    plot(newData(ii,1),newData(ii,2),'ro','MarkerSize',12)
end

for ii = not(mydiff) % plot black circles around incorrect pts
    plot(newData(ii,1),newData(ii,2),'ko','MarkerSize',12)
end
hold off

%% Test newborn data
close all
addpath(genpath('D:\Edgar\biolearning'))
resultsFolder = 'D:\Edgar\OIS_Results\svm';
load(fullfile(resultsFolder,'LPS_svm.mat'))

%% Update minimization function
c = cvpartition(group,'k',10);

minfn = @(z)crossval('mcr',xdata,group,'Predfun', ...
    @(xtrain,ytrain,xtest)crossfun(xtrain,ytrain,...
    xtest,exp(z(1)),exp(z(2))),'partition',c);

%% Search for the best parameters (run a few times to get minimal z values)
opts = optimset('TolX',5e-6,'TolFun',5e-6);
[searchmin, fval] = fminsearch(minfn,randn(2,1),opts)
z = exp(searchmin)

%% k-fold crossval

cvFolds = crossvalind('Kfold', group, k);   %# get indices of 10-fold CV
cp = classperf(group);                      %# init performance tracker
for i = 1:k                                 %# for each fold
    testIdx = (cvFolds == i);               %# get indices of test instances
    trainIdx = ~testIdx;                    %# get indices training instances
    
    xtrain = xdata(trainIdx,:);
    ytrain = group(trainIdx);
    xtest = xdata(testIdx,:);
    ytest = group(testIdx);
    
    % Update minimization function
    c = cvpartition(group,'k',10);
    minfn = @(z)crossval('mcr',xdata,group,'Predfun', ...
        @(xtrain,ytrain,xtest)crossfun(xtrain,ytrain,...
        xtest,exp(z(1)),exp(z(2))),'partition',c);
    
    % Search for the best parameters (run a few times to get minimal z values)
    opts = optimset('TolX',5e-6,'TolFun',5e-6);
    [searchmin, fval] = fminsearch(minfn,randn(2,1),opts);
    z = exp(searchmin);
    SigmaVals = z(1);
    BoxVals = z(2);
    % rbf_sigma = z(1); %0.8
    rbf_sigma = SigmaVals;
    % boxconstraint = z(2); % 2e-1
    boxconstraint = BoxVals;
    %# train an SVM model over training instances
    svmModel = svmtrain(xtrain, ytrain, ...
        'Autoscale',true, 'Showplot',false, 'Method','SMO', ...
        'BoxConstraint',boxconstraint, 'Kernel_Function','rbf', ...
        'RBF_Sigma',rbf_sigma); % , 'tolkkt', 1e-6
    
    %# test using test instances, save output labels to create confusion matrix
    pred{i} = svmclassify(svmModel, xtest, 'Showplot',false);
    % Keep target labels to create confusion matrix
    groundTruth{i} = ytest;
    
    %# evaluate and update performance object
    cp = classperf(cp, pred{i}, testIdx);
end

%# get accuracy
fprintf('Accuracy  = %0.2f %%\n', 100*cp.CorrectRate);
% accMat(iSigma, iBox) = cp.CorrectRate;
% disp(cp.ClassLabels');
%# get confusion matrix
%# columns:actual, rows:predicted, last-row: unclassified instances
% disp(cp.CountingMatrix);
fprintf('\t\t%s\t%s\n %s\t\t%d\t%d\n %s\t\t%d\t%d\n Unclass\t%d\t%d\n',cp.ClassLabels{1}, cp.ClassLabels{2}, cp.ClassLabels{1},...
    cp.CountingMatrix(1,:), cp.ClassLabels{2}, cp.CountingMatrix(2,:), cp.CountingMatrix(3,:));

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
% Plot confusion matrix
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