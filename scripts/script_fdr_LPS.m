%% load data
% load('D:\Edgar\Data\IOS_Carotid_Res\group_results_redo\group_corr_pair_seeds.mat')
load('D:\Edgar\OIS_Results\groupTest1LPS\group_corr_pair_seeds.mat')
% addpath(genpath('D:\Edgar\conn'))

%% FDR
alfa = 0.05;
FDR_HbO.p = ioi_fdr(cell2mat(statTest.w.P(:,5)));
FDR_HbO.H = FDR_HbO.p < alfa;
FDR_HbR.p = ioi_fdr(cell2mat(statTest.w.P(:,6)));
FDR_HbR.H = FDR_HbR.p < alfa;


%% Plot homotopic functional connectivity values
close all
LPS_HbO = [];
LPS_HbR = [];
NaCl_HbO = [];
NaCl_HbR = [];
for iSeeds = 1:size(groupCorrData,1)
    LPS_HbO = [LPS_HbO; groupCorrData{iSeeds,5}(isTreatment)];
    LPS_HbR = [LPS_HbR; groupCorrData{iSeeds,6}(isTreatment)];
    NaCl_HbO = [NaCl_HbO; groupCorrData{iSeeds,5}(~isTreatment)];
    NaCl_HbR = [NaCl_HbR; groupCorrData{iSeeds,6}(~isTreatment)];
end
h=figure; set(h,'color','w')
hold on
plot(ones(size(NaCl_HbO)),NaCl_HbO,'ro','MarkerSize',12,'LineWidth',2)
plot(ones(size(LPS_HbO)),LPS_HbO,'kx','MarkerSize',12,'LineWidth',2)
plot(2*ones(size(NaCl_HbR)),NaCl_HbR,'bo','MarkerSize',12,'LineWidth',2)
plot(2*ones(size(LPS_HbR)),LPS_HbR,'kx','MarkerSize',12,'LineWidth',2)
ylabel('z(r)','FontSize',14)
xlabel('Homotopic fc','FontSize',14);
xlim([0 2])
set(gca,'XTick',[0 1 2])
set(gca,'FontSize', 12)
set(gca,'XTickLabel',{'' 'HbO' 'HbR'})
legend({'NaCl_{HbO_2}' 'LPS' 'NaCl_{HbR}' 'LPS'},'Location','NorthWest')
% EOF
