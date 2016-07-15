%% load data
clear; clc;
% load('D:\Edgar\OIS_Results\groupTest1LPS\group_corr_pair_seeds.mat')
load('D:\Edgar\OIS_Results\networkResOut\results_S01_HbO.mat')
% addpath(genpath('D:\Edgar\conn'))

%% Extract Z
ZNaCl = results.Z(:,:,controlGroupIdx);
ZLPS = results.Z(:,:,treatmentGroupIdx);
nNaCl = size(ZNaCl, 3);
nLPS = size(ZLPS, 3);
ZNaClVec = [];
ZLPSVec = [];

%% Make multiple comparisons
nComp = numel(nonzeros(triu(ZNaCl(:,:,1), 1)'));
maskROI = (triu(ZNaCl(:,:,1), 1))&true;
p = nan(size(maskROI));
H = zeros(size(maskROI));
STATS = struct([]);
alpha = 0.05;
for iRows = 1:size(ZLPS, 1)
    for iCols = 1:size(ZLPS, 2)
        %         [iRows, iCols]
        if maskROI(iRows, iCols)
            [p(iRows, iCols), H(iRows, iCols)] = ranksum...
                (squeeze(ZNaCl(iRows, iCols, :)), squeeze(ZLPS(iRows, iCols, :)), 'alpha', alpha);
        end
    end
end
figure;imagesc(p, [0 alpha]);  title('p-val'); axis image; colorbar; 
colormap(ioi_get_colormap('redbluecmap'))

%% FDR
Qvec = ioi_fdr(nonzeros(triu(p, 1)'));
Q = nan(size(maskROI));
Q(maskROI) = Qvec;
Hfdr = Q < alpha;
figure;imagesc(Q, [0 alpha]);  title('FDR corrected p-val'); axis image; colorbar; 
colormap(ioi_get_colormap('redbluecmap'))

%% Plot homotopic functional connectivity values
close all
for iLPS = 1:nLPS,
    ZLPSVec = [ZLPSVec; nonzeros(triu(ZLPS(:,:,iLPS), 1)')];
end
for iNaCl = 1:nNaCl,
    ZNaClVec = [ZNaClVec; nonzeros(triu(ZNaCl(:,:,iNaCl), 1)')];
end

h=figure; set(h,'color','w')
hold on
plot(ones(size(ZNaClVec)),ZNaClVec,'ro','MarkerSize',12,'LineWidth',2)
plot(ones(size(ZLPSVec)),ZLPSVec,'kx','MarkerSize',12,'LineWidth',2)

load('D:\Edgar\OIS_Results\networkResOut\results_S01_HbR.mat')
% addpath(genpath('D:\Edgar\conn'))

% Extract Z for HbR
ZNaCl = results.Z(:,:,controlGroupIdx);
ZLPS = results.Z(:,:,treatmentGroupIdx);

for iLPS = 1:nLPS,
    ZLPSVec = [ZLPSVec; nonzeros(triu(ZLPS(:,:,iLPS), 1)')];
end
for iNaCl = 1:nNaCl,
    ZNaClVec = [ZNaClVec; nonzeros(triu(ZNaCl(:,:,iNaCl), 1)')];
end

plot(2*ones(size(ZNaClVec)),ZNaClVec,'bo','MarkerSize',12,'LineWidth',2)
plot(2*ones(size(ZLPSVec)),ZLPSVec,'kx','MarkerSize',12,'LineWidth',2)
ylabel('z(r)','FontSize',14)
xlabel('fc values','FontSize',14);
xlim([0 2])
set(gca,'XTick',[0 1 2])
set(gca,'FontSize', 12)
set(gca,'XTickLabel',{'' 'HbO' 'HbR'})
legend({'NaCl_{HbO_2}' 'LPS' 'NaCl_{HbR}' 'LPS'},'Location','NorthWest')

%% Observe repeatability in LPS
h1=figure; set(h,'color','w')
hold on
x(1:2:9) = 1;
x(2:2:10) = 2;
y(1:2:9) = LPS_HbO(1:5);
y(2:2:10) = LPS_HbO(6:10);
plot(x,y,'r-o','MarkerSize',12,'LineWidth',2)
title('LPS01','FontSize',14)
ylabel('z(r)','FontSize',14)
xlabel('Homotopic fc','FontSize',14);
xlim([0 2])
set(gca,'XTick',[0 1 2])
set(gca,'FontSize', 12)
set(gca,'XTickLabel',{'' 'LPS_{session_1}' 'LPS_{session_2}'})

%% Observe repeatability in NaCl
h2=figure; set(h,'color','w')
hold on
x(1:2:9) = 1;
x(2:2:10) = 2;
y(1:2:9) = LPS_HbO(1:5);
y(2:2:10) = LPS_HbO(6:10);
plot(x,y,'r-o','MarkerSize',12,'LineWidth',2)
title('LPS01','FontSize',14)
ylabel('z(r)','FontSize',14)
xlabel('Homotopic fc','FontSize',14);
xlim([0 2])
set(gca,'XTick',[0 1 2])
set(gca,'FontSize', 12)
set(gca,'XTickLabel',{'' 'LPS_{session_1}' 'LPS_{session_2}'})

% EOF
