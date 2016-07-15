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
ZLPSVec=[];
ZNaClVec=[];
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

ZLPSVec=[];
ZNaClVec=[];

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

%% Plot homotopic fc values as a function of CO2
NaClCO2 = [59.2; 49.4; 41.9; NaN; 62.6; NaN; 56.7; 50];
LPSCO2 = [NaN; 39.5; 54; 40.3; 80.3];
bilatROIsIdx = [(1:2:10)' (2:2:10)'];

load('D:\Edgar\OIS_Results\networkResOut\results_S01_HbR.mat')
% Extract Z for HbR
ZNaCl = results.Z(:,:,controlGroupIdx);
ZLPS = results.Z(:,:,treatmentGroupIdx);
ZLPSVec=[];
ZNaClVec=[];
for iLPS = 1:nLPS,
    for iROI = 1:size(bilatROIsIdx,1)
        ZLPSVec = [ZLPSVec; ZLPS(bilatROIsIdx(iROI, 1),bilatROIsIdx(iROI, 2), iLPS)];
    end
end
for iNaCl = 1:nNaCl,
    for iROI = 1:size(bilatROIsIdx,1)
        ZNaClVec = [ZNaClVec; ZNaCl(bilatROIsIdx(iROI, 1),bilatROIsIdx(iROI, 2), iNaCl)];
    end
end

figure; set(gcf,'color','w')
hold on
% Plot HbR
plot(NaClCO2,reshape(ZNaClVec, [numel(NaClCO2) numel(ZNaClVec)/size(ZNaCl,3)]),...
    'bo','MarkerSize',12,'LineWidth',2)
plot(LPSCO2,reshape(ZLPSVec, [numel(LPSCO2) numel(ZLPSVec)/size(ZLPS,3)]),...
    'bx','MarkerSize',12,'LineWidth',2)

load('D:\Edgar\OIS_Results\networkResOut\results_S01_HbO.mat')
% addpath(genpath('D:\Edgar\conn'))

% Extract Z for HbO
ZNaCl = results.Z(:,:,controlGroupIdx);
ZLPS = results.Z(:,:,treatmentGroupIdx);

ZLPSVec=[];
ZNaClVec=[];
for iLPS = 1:nLPS,
    for iROI = 1:size(bilatROIsIdx,1)
        ZLPSVec = [ZLPSVec; ZLPS(bilatROIsIdx(iROI, 1),bilatROIsIdx(iROI, 2), iLPS)];
    end
end
for iNaCl = 1:nNaCl,
    for iROI = 1:size(bilatROIsIdx,1)
        ZNaClVec = [ZNaClVec; ZNaCl(bilatROIsIdx(iROI, 1),bilatROIsIdx(iROI, 2), iNaCl)];
    end
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

% EOF
