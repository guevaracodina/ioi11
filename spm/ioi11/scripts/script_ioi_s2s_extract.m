%% script_ioi_add_physio
addpath(genpath('C:\spm8\toolbox\ioi11'))
clear; close all; clc
topDir{1} = 'C:\Edgar\Data\IOIResults\j02\ROItest\GLMheart\';
topDir{2} = 'C:\Edgar\Data\IOIResults\j03\ROItest\GLMheart\';
topDir{3} = 'C:\Edgar\Data\IOIResults\j04\ROItest\GLMheart\';
% topDir{4} = 'C:\Edgar\Data\IOIResults\j05s01\'; IMPOSSIBLE TO EXTRACT
topDir{4} = 'C:\Edgar\Data\IOIResults\j06\ROItest\GLMheart\';
% topDir{6} = 'C:\Edgar\Data\IOIResults\j07\';  IMPOSSIBLE TO EXTRACT
topDir{5} = 'C:\Edgar\Data\IOIResults\k01\ROItest\GLMheart\';
topDir{6} = 'C:\Edgar\Data\IOIResults\k02\ROItest\GLMheart\';
topDir{7} = 'C:\Edgar\Data\IOIResults\k03\ROItest\GLMheart\';
topDir{8} = 'C:\Edgar\Data\IOIResults\k04\ROItest\GLMheart\';
topDir{9} = 'C:\Edgar\Data\IOIResults\k05\ROItest\GLMheart\';
topDir{10} = 'C:\Edgar\Data\IOIResults\k06s01\ROItest\GLMheart\';
topDir{11} = 'C:\Edgar\Data\IOIResults\L01\ROItest\GLMheart\';
topDir{12} = 'C:\Edgar\Data\IOIResults\L02\ROItest\GLMheart\';
topDir{13} = 'C:\Edgar\Data\IOIResults\L03\ROItest\GLMheart\';
topDir{14} = 'C:\Edgar\Data\IOIResults\L04\ROItest\GLMheart\';
topDir{15} = 'C:\Edgar\Data\IOIResults\L05\ROItest\GLMheart\';
topDir{16} = 'C:\Edgar\Data\IOIResults\L06\ROItest\GLMheart\';
topDir{17} = 'C:\Edgar\Data\IOIResults\L07\ROItest\GLMheart\';
topDir{18} = 'C:\Edgar\Data\IOIResults\L08\ROItest\GLMheart\';

currentDir = 'corrMap';

%%
clc
motorBilatCorr = [];
somatoBilatCorr = [];
for iDirs = 1:numel(topDir)
    % Choose heart rate data
%     [dirPhysio, sts] = cfg_getfile([1 1],'dir','Select folder',{fullfile(topDir{iDirs},currentDir)}, topDir, '.*');
    dirListIOI = dir(fullfile(topDir{iDirs},[currentDir filesep 'IOI.mat']));
%     [heartRateFile, sts] = cfg_getfile(1,'mat','Select heart rate file',[], dirPhysio{1}, 'heartRate*');
    
    % Backup IOI
    IOImat = fullfile(topDir{iDirs},[currentDir filesep dirListIOI.name]);
%     copyfile(IOImat,...
%         fullfile(topDir{iDirs},[currentDir filesep 'IOI.bak']));
%     
    % Point to the heart rate file and save IOI mat
    load(IOImat);
%     IOI.fcIOS.SPM(1).physioHeartFile = heartRateFile;
%     save(IOImat, 'IOI')
    load(IOI.fcIOS.corr.corrMatrixFname)
    fprintf('%s: Motor: %.10f Somato: %.10f \n',IOI.subj_name, seed2seedCorrMat{1}{1}(1,2), seed2seedCorrMat{1}{1}(3,4));
    motorBilatCorr(iDirs,1) = fisherz(seed2seedCorrMat{1}{1}(1,2));
    somatoBilatCorr(iDirs,1) = fisherz(seed2seedCorrMat{1}{1}(3,4));
end
%% LPS & NaCl
% All J are LPS + L05, L06, L07, L08
idxLPS = [(1:4)'; (15:18)'];
idxNaCl = (5:14)';
LPSmotor = motorBilatCorr(idxLPS);
LPSsomato = somatoBilatCorr(idxLPS);
NaCLmotor = motorBilatCorr(idxNaCl);
NaCLsomato = somatoBilatCorr(idxNaCl);

%% Distribution plots
addpath(genpath('C:\Edgar\Dropbox\Matlab'))
addpath(genpath('C:\spm8\toolbox\ioi11'))
figure; set(gcf,'color','w');
distributionPlot([[LPSmotor; mean(LPSmotor); mean(LPSmotor)] NaCLmotor...
    [LPSsomato; mean(LPSsomato); mean(LPSsomato)] NaCLsomato],'color',...
    {[0.8500    0.3250    0.0980] [0    0.4470    0.7410] ...
    [0.8500    0.3250    0.0980] [0    0.4470    0.7410]},...
    'addSpread',true,'showMM',6,'variableWidth',true);
% set(gca, 'XTickLabel',{'LPS' 'NaCl' 'LPS' 'NaCl'}, 'FontSize',14)       
set(gca, 'XTick',[1.5 3.5])       
set(gca, 'XTickLabel',{'Motor' 'Somatosensory'}, 'FontSize',14)       
ylabel('Bilateral Correlations {\it z(r)}','FontSize',14)
legend({'LPS' 'NaCl'},'FontSize',14)
ylim([-2 2])
[pMotor,hMotor,statsMotor] = ranksum(LPSmotor,NaCLmotor);
[pSomato,hSomato,statsSomato] = ranksum(LPSsomato,NaCLsomato);
fprintf('There was no statistical difference between LPS and NaCl, neither in the motor region (p=%.4f), nor in the somatosensory cortex (p=%.4f)\n', pMotor, pSomato);
% Effect size Cohen's d
dMotor = 0.39107874467381487;
rMotor = 0.1919049817838292;
dSomato = 0.22023315516401254;
rSomato = 0.10945497023211405;

%% Considering the degree of lesion (motor cortex)
% Column 1: Lesion degree
% Column 2: Motor Bilateral correlation
Xmotor = zeros(numel(NaCLmotor)+numel(LPSmotor), 2);
% NaCl, label = 0
Xmotor(1:numel(NaCLmotor),1) = 0;
Xmotor(1:numel(NaCLmotor),2) = NaCLmotor;
% Low lesion
Xmotor(numel(NaCLmotor)+1,1) = 1;
Xmotor(numel(NaCLmotor)+2,1) = 1;
Xmotor(numel(NaCLmotor)+1,2) = motorBilatCorr(2); % j03
Xmotor(numel(NaCLmotor)+2,2) = motorBilatCorr(4); % j06
% Moderate lesion
Xmotor(numel(NaCLmotor)+3,1) = 2;
Xmotor(numel(NaCLmotor)+4,1) = 2;
Xmotor(numel(NaCLmotor)+5,1) = 2;
Xmotor(numel(NaCLmotor)+3,2) = motorBilatCorr(1); % j02
Xmotor(numel(NaCLmotor)+4,2) = motorBilatCorr(15); % l05
Xmotor(numel(NaCLmotor)+5,2) = motorBilatCorr(16); % l06
% High lesion
Xmotor(numel(NaCLmotor)+6,1) = 3;
Xmotor(numel(NaCLmotor)+7,1) = 3;
Xmotor(numel(NaCLmotor)+8,1) = 3;
Xmotor(numel(NaCLmotor)+6,2) = motorBilatCorr(17); % l07
Xmotor(numel(NaCLmotor)+7,2) = motorBilatCorr(18); % l08
Xmotor(numel(NaCLmotor)+8,2) = motorBilatCorr(3); % j04

%% plot data
plot(Xmotor(:,1), Xmotor(:,2), 'k*')
xlabel('Degree of lesion')
ylabel('Bilateral correlation z(r)')
xlim([-0.5 3.5]);
ylim([-1 1.5]);

%% k-means
k = 4;
% idx = kmeans(Xmotor,k);
opts = statset('Display','final');
[idx,C] = kmeans(Xmotor,k,'Distance','cityblock',...
    'Replicates',15,'Options',opts);
figure;
plot(Xmotor(idx==1,1),Xmotor(idx==1,2),'r*','MarkerSize',12,'LineWidth',2)
hold on
plot(Xmotor(idx==2,1),Xmotor(idx==2,2),'b*','MarkerSize',12,'LineWidth',2)
hold on
plot(Xmotor(idx==3,1),Xmotor(idx==3,2),'m*','MarkerSize',12,'LineWidth',2)
hold on
plot(Xmotor(idx==4,1),Xmotor(idx==4,2),'c*','MarkerSize',12,'LineWidth',2)

plot(C(:,1),C(:,2),'kx',...
     'MarkerSize',15,'LineWidth',3)
xlim([-0.5 3.5]);
ylim([-1 1.5]);
legend({'Cluster 1','Cluster 2','Cluster 3','Cluster 4','Centroids'},...
       'FontSize', 12,'Location','NW')
title ('Cluster Assignments and Centroids', 'FontSize', 12)
xlabel('Degree of lesion', 'FontSize', 12)
ylabel('Bilateral correlation z(r)', 'FontSize', 12)
set(gca, 'FontSize', 12, 'XTick', [0 1 2 3])
hold off

%% Considering the degree of lesion (motor cortex)
% Column 1: Lesion degree
% Column 2: Motor Bilateral correlation
Xsomato = zeros(numel(NaCLsomato)+numel(LPSsomato), 2);
% NaCl, label = 0
Xsomato(1:numel(NaCLsomato),1) = 0;
Xsomato(1:numel(NaCLsomato),2) = NaCLsomato;
% Low lesion
Xsomato(numel(NaCLsomato)+1,1) = 1;
Xsomato(numel(NaCLsomato)+2,1) = 1;
Xsomato(numel(NaCLsomato)+1,2) = somatoBilatCorr(2); % j03
Xsomato(numel(NaCLsomato)+2,2) = somatoBilatCorr(4); % j06
% Moderate lesion
Xsomato(numel(NaCLsomato)+3,1) = 2;
Xsomato(numel(NaCLsomato)+4,1) = 2;
Xsomato(numel(NaCLsomato)+5,1) = 2;
Xsomato(numel(NaCLsomato)+3,2) = somatoBilatCorr(1); % j02
Xsomato(numel(NaCLsomato)+4,2) = somatoBilatCorr(15); % l05
Xsomato(numel(NaCLsomato)+5,2) = somatoBilatCorr(16); % l06
% High lesion
Xsomato(numel(NaCLsomato)+6,1) = 3;
Xsomato(numel(NaCLsomato)+7,1) = 3;
Xsomato(numel(NaCLsomato)+8,1) = 3;
Xsomato(numel(NaCLsomato)+6,2) = somatoBilatCorr(17); % l07
Xsomato(numel(NaCLsomato)+7,2) = somatoBilatCorr(18); % l08
Xsomato(numel(NaCLsomato)+8,2) = somatoBilatCorr(3); % j04

%% plot data
plot(Xsomato(:,1), Xsomato(:,2), 'k*')
xlabel('Degree of lesion')
ylabel('Bilateral correlation z(r)')
xlim([-0.5 3.5]);
ylim([-1 1.5]);

%% k-means
k = 4;
% idx = kmeans(Xsomato,k);
opts = statset('Display','final');
[idx,C] = kmeans(Xsomato,k,'Distance','cityblock',...
    'Replicates',15,'Options',opts);
figure;
plot(Xsomato(idx==1,1),Xsomato(idx==1,2),'r*','MarkerSize',12,'LineWidth',2)
hold on
plot(Xsomato(idx==2,1),Xsomato(idx==2,2),'b*','MarkerSize',12,'LineWidth',2)
hold on
plot(Xsomato(idx==3,1),Xsomato(idx==3,2),'m*','MarkerSize',12,'LineWidth',2)
hold on
plot(Xsomato(idx==4,1),Xsomato(idx==4,2),'c*','MarkerSize',12,'LineWidth',2)

plot(C(:,1),C(:,2),'kx',...
     'MarkerSize',15,'LineWidth',3)
xlim([-0.5 3.5]);
ylim([-1 1.5]);
legend({'Cluster 1','Cluster 2','Cluster 3','Cluster 4','Centroids'},...
       'FontSize', 12,'Location','NW')
title ('Cluster Assignments and Centroids', 'FontSize', 12)
xlabel('Degree of lesion', 'FontSize', 12)
ylabel('Bilateral correlation z(r)', 'FontSize', 12)
set(gca, 'FontSize', 12, 'XTick', [0 1 2 3])
hold off