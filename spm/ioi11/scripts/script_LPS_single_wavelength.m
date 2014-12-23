%% script_LPS_single_wavelength
clear; close all; clc
% LPS group
seed2seedLPS = {
                'C:\Edgar\Data\IOIResults\j02\ROItest\GLMfcIOS\corrMap\s2sCorrMat.mat'
                'C:\Edgar\Data\IOIResults\j03\ROItest\GLMfcIOS\corrMap\s2sCorrMat.mat'
                'C:\Edgar\Data\IOIResults\j04\ROItest\GLMfcIOS\corrMap\s2sCorrMat.mat'
                'C:\Edgar\Data\IOIResults\j05s01\ROItest\GLMfcIOS\corrMap\s2sCorrMat.mat'
                %'C:\Edgar\Data\IOIResults\j05s02\ROItest\GLMfcIOS\corrMap\s2sCorrMat.mat'
                'C:\Edgar\Data\IOIResults\j06\ROItest\GLMfcIOS\corrMap\s2sCorrMat.mat'
                'C:\Edgar\Data\IOIResults\j07\ROItest\GLMfcIOS\corrMap\s2sCorrMat.mat'
                };
% Control (NaCl) group  
seed2seedCtrl = {  
                'C:\Edgar\Data\IOIResults\k01\ROItest\GLMfcIOS\corrMap\s2sCorrMat.mat'
                'C:\Edgar\Data\IOIResults\k02\ROItest\GLMfcIOS\corrMap\s2sCorrMat.mat'
                'C:\Edgar\Data\IOIResults\k03\ROItest\GLMfcIOS\corrMap\s2sCorrMat.mat'
                'C:\Edgar\Data\IOIResults\k04\ROItest\GLMfcIOS\corrMap\s2sCorrMat.mat'
                'C:\Edgar\Data\IOIResults\k05\ROItest\GLMfcIOS\corrMap\s2sCorrMat.mat'
                'C:\Edgar\Data\IOIResults\k06s01\ROItest\GLMfcIOS\corrMap\s2sCorrMat.mat'
                % 'C:\Edgar\Data\IOIResults\k06s02\ROItest\GLMfcIOS\corrMap\s2sCorrMat.mat'
                };
%% Total population            
seed2seedTotal = {
                'C:\Edgar\Data\IOIResults\j02\ROItest\GLMfcIOS\corrMap\s2sCorrMat.mat'
                'C:\Edgar\Data\IOIResults\j03\ROItest\GLMfcIOS\corrMap\s2sCorrMat.mat'
                'C:\Edgar\Data\IOIResults\j04\ROItest\GLMfcIOS\corrMap\s2sCorrMat.mat'
                'C:\Edgar\Data\IOIResults\j05s01\ROItest\GLMfcIOS\corrMap\s2sCorrMat.mat'
                %'C:\Edgar\Data\IOIResults\j05s02\ROItest\GLMfcIOS\corrMap\s2sCorrMat.mat'
                'C:\Edgar\Data\IOIResults\j06\ROItest\GLMfcIOS\corrMap\s2sCorrMat.mat'
                'C:\Edgar\Data\IOIResults\j07\ROItest\GLMfcIOS\corrMap\s2sCorrMat.mat'
                'C:\Edgar\Data\IOIResults\k01\ROItest\GLMfcIOS\corrMap\s2sCorrMat.mat'
                'C:\Edgar\Data\IOIResults\k02\ROItest\GLMfcIOS\corrMap\s2sCorrMat.mat'
                'C:\Edgar\Data\IOIResults\k03\ROItest\GLMfcIOS\corrMap\s2sCorrMat.mat'
                'C:\Edgar\Data\IOIResults\k04\ROItest\GLMfcIOS\corrMap\s2sCorrMat.mat'
                'C:\Edgar\Data\IOIResults\k05\ROItest\GLMfcIOS\corrMap\s2sCorrMat.mat'
                'C:\Edgar\Data\IOIResults\k06s01\ROItest\GLMfcIOS\corrMap\s2sCorrMat.mat'
                % 'C:\Edgar\Data\IOIResults\k06s02\ROItest\GLMfcIOS\corrMap\s2sCorrMat.mat'
                'C:\Edgar\Data\IOIResults\L01\ROItest\GLMfcIOS\corrMap\s2sCorrMat.mat'
                'C:\Edgar\Data\IOIResults\L02\ROItest\GLMfcIOS\corrMap\s2sCorrMat.mat'
                'C:\Edgar\Data\IOIResults\L03\ROItest\GLMfcIOS\corrMap\s2sCorrMat.mat'
                'C:\Edgar\Data\IOIResults\L04\ROItest\GLMfcIOS\corrMap\s2sCorrMat.mat'
                'C:\Edgar\Data\IOIResults\L05\ROItest\GLMfcIOS\corrMap\s2sCorrMat.mat'
                'C:\Edgar\Data\IOIResults\L06\ROItest\GLMfcIOS\corrMap\s2sCorrMat.mat'
                'C:\Edgar\Data\IOIResults\L07\ROItest\GLMfcIOS\corrMap\s2sCorrMat.mat'
                'C:\Edgar\Data\IOIResults\L08\ROItest\GLMfcIOS\corrMap\s2sCorrMat.mat'
                };
nTot = numel(seed2seedTotal);
seed2seedMatTot = zeros([nTot 2]);
for i=1:nTot,
    load(seed2seedTotal{i});
    seed2seedMatTot(i,1) = seed2seedCorrMat{1}{1}(1,2);
    seed2seedMatTot(i,2) = seed2seedCorrMat{1}{1}(3,4);
end
            
%% Load values
nLPS = numel(seed2seedLPS);
nCtrl = numel(seed2seedCtrl);
seed2seedMatLPS = zeros([nLPS 2]);
seed2seedMatCtrl = zeros([nCtrl 2]);
for i=1:nLPS,
    load(seed2seedLPS{i});
    seed2seedMatLPS(i,1) = seed2seedCorrMat{1}{1}(1,2);
    seed2seedMatLPS(i,2) = seed2seedCorrMat{1}{1}(3,4);
    load(seed2seedCtrl{i});
    seed2seedMatCtrl(i,1) = seed2seedCorrMat{1}{1}(1,2);
    seed2seedMatCtrl(i,2) = seed2seedCorrMat{1}{1}(3,4);
end

%% Convert to Fisher z, do stats
alpha = 0.05;
seed2seedMatLPSZ = fisherz(seed2seedMatLPS);
seed2seedMatCtrlZ = fisherz(seed2seedMatCtrl);
[motor.P, motor.H, motor.STATS] = ranksum...
                        (seed2seedMatLPS(:,1), seed2seedMatCtrl(:,1), 'alpha', alpha);
motor.id = 'Wilcoxon rank sum test(raw data)';
[somato.P, somato.H, somato.STATS] = ranksum...
                        (seed2seedMatLPS(:,2), seed2seedMatCtrl(:,2), 'alpha', alpha);
somato.id = 'Wilcoxon rank sum test(raw data)';

% mean
motor.LPS.avg = mean(seed2seedMatLPS(:,1));
motor.Ctrl.avg = mean(seed2seedMatCtrl(:,1));
somato.LPS.avg = mean(seed2seedMatLPS(:,2));
somato.Ctrl.avg = mean(seed2seedMatCtrl(:,2));
% std dev
motor.LPS.stdev = std(seed2seedMatLPS(:,1));
motor.Ctrl.stdev = std(seed2seedMatCtrl(:,1));
somato.LPS.stdev = std(seed2seedMatLPS(:,2));
somato.Ctrl.stdev = std(seed2seedMatCtrl(:,2));
% std error
motor.LPS.stderror = motor.LPS.stdev / sqrt(nLPS);
somato.LPS.stderror = somato.LPS.stdev / sqrt(nLPS);
motor.Ctrl.stderror = motor.Ctrl.stdev / sqrt(nCtrl);
somato.Ctrl.stderror = somato.Ctrl.stdev / sqrt(nCtrl);

%% plot bar with error
e(1,:) = [motor.LPS.stderror motor.Ctrl.stderror];
e(2,:) = [somato.LPS.stderror somato.Ctrl.stderror];

y(1,:) = [motor.LPS.avg motor.Ctrl.avg];
y(2,:) = [somato.LPS.avg somato.Ctrl.avg];

figSize = [3.5 3.5];
close all;
h = figure; set(h,'color','w')
% Specify window units
set(h, 'units', 'inches')
% Change figure and paper size
set(h, 'Position', [0.1 0.1 figSize(1) figSize(2)])
set(h, 'PaperPosition', [0.1 0.1 figSize(1) figSize(2)])
% Custom bar graphs with error bars (1st arg: error)
barwitherr(e, y);

axisFontSize    = 12;
starFontSize    = 22;
figRes = 300;
colormap([0.5 0.5 0.5; 1 1 1]);
set(gca, 'XTickLabel', {'Motor' 'Somatosensory'}, 'FontSize', axisFontSize)
ylabel('Bilateral functional correlation z(r)', 'FontSize', axisFontSize);
legend({'LPS' 'NaCl'});
text(0.98, 0.7, '*', 'FontSize', starFontSize, 'FontWeight', 'b');
set(gca,'FontSize', axisFontSize)

%%
% Save as PNG

%print(h, '-dpng', fullfile('C:\Edgar\Dropbox\PostDoc\Newborn','z_r_LPS_Ctrl.png'), sprintf('-r%d',figRes));
% EOF