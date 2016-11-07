%% Load connectivity data
% resultsFolder = 'C:\Users\Ramón\Desktop\Edgar\ANN'; 
% resultsFolder = 'D:\Edgar\OIS_Results\ANN';
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
% load('D:\Edgar\OIS_Results\networkResOut\resultsROI_Condition01_HbR.mat')
% load('C:\Users\Ramón\Desktop\Edgar\networkResOutNoVis\results_S01_HbR.mat');
% load('C:\Users\Ramón\Desktop\Edgar\networkResOutNoVis\resultsROI_Condition01_HbR.mat')
load('C:\Edgar\Dropbox\PostDoc\Newborn\OIS_Results\networkResOut\results_S01_HbR.mat')
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
% load('D:\Edgar\OIS_Results\networkResOut\resultsROI_Condition01_HbO.mat')
% load('C:\Users\Ramón\Desktop\Edgar\networkResOutNoVis\results_S01_HbO.mat')
% load('C:\Users\Ramón\Desktop\Edgar\networkResOutNoVis\resultsROI_Condition01_HbO.mat')
load('C:\Edgar\Dropbox\PostDoc\Newborn\OIS_Results\networkResOutNoVis\results_S01_HbO.mat')
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

%% ANN data preparation with newborn data

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
% LP01; LP02, LP03, LP04 & LP05
% LPShist.L = [2002507; 2149179; 957905; 6718538];
% LPShist.R = [168662; 1142728; 909365; 6718538];
% LPShist.T = LPShist.L + LPShist.R;
LPShist.T = [11.98; 4.00; 5.98; 3.86; 20.58];
% NC06, NC07, NC08 & NC09
% NaClhist.L = [122725; 135333; 194805];
% NaClhist.R = [110830; 159154; 198182];
% NaClhist.T = NaClhist.L + NaClhist.R;
NaClhist.T = [0.58; 0.70; 0.53; 1.01];
% Indices of subjects to keep
idx2keep = [5:8, 9:13]';

%% Crop X-data
groupLabels = {'NC01'; 'NC02'; 'NC03'; 'NC05'; 'NC06'; 'NC07'; 'NC08'; 'NC09';...
    'LP01'; 'LP02'; 'LP03'; 'LP04'; 'LP05'};
groupLabels = groupLabels(idx2keep);
xdata = xdata(idx2keep,:);
ydata = [NaClhist.T; LPShist.T];
group = group(idx2keep);
myVarsPred = [xdata ydata];
myVarsClass = [xdata strcmp(group, 'NaCl')];

% clearvars -except xdata ydata group myVarsPred myVarsClass

%% k-fold Cross validation using ANN classifier
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
    pred = [pred; myNeuralNetworkFunction(xtest)];

    % Keep target values to create correlation plot
    groundTruth = [groundTruth; ytest];
    
    % Keep target labels 
    targetLabels = [targetLabels; groupLabels(testIdx)];
end
        
%% Plot correlation between measured and predicted values
p = polyfit(groundTruth, pred, 1);
yfit = p(1)*groundTruth + p(2);
[rho, pVal] = corr(groundTruth, pred);
RMSEP = sqrt(sum((groundTruth-pred).^2)/size(xdata,1));

% Find indices of LPS
idxLPS = find(~cellfun(@isempty, regexp(targetLabels,'LP')));
idxNaCl = find(~cellfun(@isempty, regexp(targetLabels,'NC')));

hFig = figure; hold on;
% Plot Predicted = Measured
plot(groundTruth, groundTruth, 'k:', 'LineWidth', 3)
% Plot fit
plot(groundTruth, yfit, 'b-', 'LineWidth', 3)
% Plot NaCl
plot(groundTruth(idxNaCl), pred(idxNaCl), 'ko', 'LineWidth', 3, 'MarkerSize', 12)
% Plot LPS
plot(groundTruth(idxLPS), pred(idxLPS), 'rx', 'LineWidth', 3, 'MarkerSize', 12)

axis equal
% legend({'Predicted = Measured' 'Linear fit' 'NaCl' 'LPS' }, 'Location', 'SouthEast')
set(gca, 'FontSize', 14)
xlim([min(groundTruth) max(groundTruth)])
ylim([min(groundTruth) max(groundTruth)])
xlabel('Measured (%)', 'FontSize', 14); 
ylabel('Predicted(pixels)', 'FontSize', 14)
title(sprintf('Fractional Lesion Volume'), 'FontSize', 14); 
text(1,18, sprintf('RMSEP=%0.4f%%\nr = %0.4f\np = %0.4e',RMSEP, rho, pVal), 'FontSize', 14)

% Specify window units
set(hFig, 'units', 'inches')
% Change figure and paper size
set(hFig, 'Position', [0.1 0.1 3.5 3.5])
set(hFig, 'PaperPosition', [0.1 0.1 3.5 3.5])

if false
    % Save as PNG at the user-defined resolution
    print(hFig, '-dpng', ...
        fullfile(resultsFolder,...
        sprintf('ANN_pred_meas_ventricular_frac.png')),...
        sprintf('-r%d',300));
    % Return the property to its default
    set(hFig, 'units', 'pixels')
    close(hFig)
end

%% ANN diagram
% %# neural net, and view it
% load('C:\Edgar\Dropbox\PostDoc\Newborn\OIS_Results\ANN\ANNresults.mat')
% jframe = view(ANNresults.net);
% 
% %# create it in a MATLAB figure
% % Specify window units
% hFig = figure('Menubar','none');
% set(hFig, 'units','inches', 'position',[0.1 0.1 6.5 2])
% % hFig = figure('Menubar','none');
% jpanel = get(jframe,'ContentPane');
% [~,h] = javacomponent(jpanel);
% set(h, 'units','inches', 'position',[0.1 0.1 6.5 2])
% 
% %# close java window
% jframe.setVisible(false);
% jframe.dispose();
% 
% %# print to file
% % set(hFig, 'PaperPositionMode', 'auto')
% % saveas(hFig, 'out.png')
% 
% %# close figure
% % close(hFig)
% 
% 
% % Change figure and paper size
% % set(hFig, 'Position', [0.1 0.1 3.5 1])
% set(hFig, 'PaperPosition', [0.1 0.1 6.5 2])
% 
% if false
%     % Save as PNG at the user-defined resolution
%     print(hFig, '-dpng', ...
%         fullfile(resultsFolder,...
%         sprintf('ANN_diagram.png')),...
%         sprintf('-r%d',300));
%     % Return the property to its default
%     set(hFig, 'units', 'pixels')
% %     close(hFig)
% end

%% Find correlation
% clear x y;
% i = 1;
% for iRow=1:size(ZLPSHbR,1)
%     for iCol=1:size(ZLPSHbR,2)
%         for iLPS=1:size(ZLPSHbR,3)
% %             plot(squeeze(ZLPSHbR(iRow, iCol, iLPS)), squeeze(ZNaClHbR(iRow, iCol, :)), 'k.');
% %             hold on
%             x(i,:) = squeeze(ZLPS(iRow, iCol, iLPS));
%             y(i,:) = squeeze(ZNaCl(iRow, iCol, iLPS));
%             i = i + 1;
%         end
%     end
% end
% plot(x,y,'k.')
% [rho, p] = corr(x,y,'rows','complete')


%% 
%  fc vs lesion
plot(ydata, xdata, 'ko')
% Load spatial extent (HbO)
load('D:\Edgar\OIS_Results\averaged_maps\stats_C5.mat')
LPS_extent.HbO = LPS_spatial_extension;
NaCl_extent.HbO = NaCl_spatial_extension;
% Load spatial extent (HbR)
load('D:\Edgar\OIS_Results\averaged_maps\stats_C6.mat')
LPS_extent.HbR = LPS_spatial_extension;
NaCl_extent.HbR = NaCl_spatial_extension;

% Find smallest difference in spatial extent
vol = spm_vol('D:\Edgar\OIS_Results\averaged_maps\16_02_25,NC01_anat_brainmask.nii');
brainMask = spm_read_vols(vol);
diffPixels = min(abs(NaCl_extent.HbO - LPS_extent.HbO))*nnz(brainMask);
% 70 pixels/mm
pix_per_mm = 70;
diff_mm2 = diffPixels / (pix_per_mm*pix_per_mm);

% EOF