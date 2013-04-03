%% load data
% load('D:\Edgar\Data\IOS_Carotid_Res\group_results_redo\group_corr_pair_seeds.mat')
load('D:\Edgar\Data\IOS_Carotid_Res\group_results_redo_outliers\group_corr_pair_seeds.mat')
addpath(genpath('D:\Edgar\conn'))

%% FDR
alfa = 0.05;
FDR_HbO.p = conn_fdr(cell2mat(statTest.t.P(:,5)));
FDR_HbO.H = FDR_HbO.p < alfa;
FDR_HbR.p = conn_fdr(cell2mat(statTest.t.P(:,6)));
FDR_HbR.H = FDR_HbR.p < alfa;
FDR_CBF.p = conn_fdr(cell2mat(statTest.t.P(:,7)));
FDR_CBF.H = FDR_CBF.p < alfa;
FDR_CMRO2.p = conn_fdr(cell2mat(statTest.t.P(:,8)));
FDR_CMRO2.H = FDR_CMRO2.p < alfa;

%% Derivative
load('D:\Edgar\Data\IOS_Carotid_Res\group_results_redo\group_corr_pair_seeds_diff.mat')
alfa = 0.05;
FDR_HbODiff.p = conn_fdr(cell2mat(statTestDiff.t.P(:,5)));
FDR_HbODiff.H = FDR_HbODiff.p < alfa;
FDR_HbRDiff.p = conn_fdr(cell2mat(statTestDiff.t.P(:,6)));
FDR_HbRDiff.H = FDR_HbRDiff.p < alfa;
FDR_CBFDiff.p = conn_fdr(cell2mat(statTestDiff.t.P(:,7)));
FDR_CBFDiff.H = FDR_CBFDiff.p < alfa;
FDR_CMRO2Diff.p = conn_fdr(cell2mat(statTestDiff.t.P(:,8)));
FDR_CMRO2Diff.H = FDR_CMRO2Diff.p < alfa;

% EOF
