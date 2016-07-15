%% load data
clear; clc;
NaClCO2 = [59.2; 49.4; 41.9; NaN; 62.6; NaN; 56.7; 50];
LPSCO2 = [NaN; 39.5; 54; 40.3; 80.3];
bilatROIsIdx = [(1:2:10)' (2:2:10)'];

load('D:\Edgar\OIS_Results\networkResOut\results_S01_HbR.mat')
% Extract Z for HbR
ZNaClHbR = results.Z(:,:,controlGroupIdx);
ZLPSHbR = results.Z(:,:,treatmentGroupIdx);
nNaCl = size(ZNaClHbR, 3);
nLPS = size(ZLPSHbR, 3);
ZLPSVecHbR=[];
ZNaClVecHbR=[];
for iLPS = 1:nLPS,
    for iROI = 1:size(bilatROIsIdx,1)
        ZLPSVecHbR = [ZLPSVecHbR; ZLPSHbR(bilatROIsIdx(iROI, 1),bilatROIsIdx(iROI, 2), iLPS)];
    end
end
for iNaCl = 1:nNaCl,
    for iROI = 1:size(bilatROIsIdx,1)
        ZNaClVecHbR = [ZNaClVecHbR; ZNaClHbR(bilatROIsIdx(iROI, 1),bilatROIsIdx(iROI, 2), iNaCl)];
    end
end

figure; set(gcf,'color','w')
hold on
% Plot HbR
plot(NaClCO2,reshape(ZNaClVecHbR, [numel(NaClCO2) numel(ZNaClVecHbR)/size(ZNaClHbR,3)]),...
    'bo','MarkerSize',12,'LineWidth',2)
plot(LPSCO2,reshape(ZLPSVecHbR, [numel(LPSCO2) numel(ZLPSVecHbR)/size(ZLPSHbR,3)]),...
    'bx','MarkerSize',12,'LineWidth',2)

load('D:\Edgar\OIS_Results\networkResOut\results_S01_HbO.mat')
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
