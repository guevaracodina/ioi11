
%% Load data
clear; clc
dataFolder = 'D:\Edgar\OIS_Results\alff';

LPSIOImat{1} = 'D:\Edgar\OIS_Results\16_02_25,LP01a\ROI\LPF\FiltNDown\GLM\corrMap\IOI.mat';
LPSIOImat{2} = 'D:\Edgar\OIS_Results\16_02_25,LP01b\ROI\LPF\FiltNDown\GLM\corrMap\IOI.mat';
LPSIOImat{3} = 'D:\Edgar\OIS_Results\16_07_07,LP02a\ROI\LPF\FiltNDown\GLM\corrMap\IOI.mat';
LPSIOImat{4} = 'D:\Edgar\OIS_Results\16_07_07,LP02b\ROI\LPF\FiltNDown\GLM\corrMap\IOI.mat';
LPSIOImat{5} = 'D:\Edgar\OIS_Results\16_07_08,LP03\ROI\LPF\FiltNDown\GLM\corrMap\IOI.mat';
LPSIOImat{6} = 'D:\Edgar\OIS_Results\16_07_08,LP04\ROI\LPF\FiltNDown\GLM\corrMap\IOI.mat';
LPSIOImat{7} = 'D:\Edgar\OIS_Results\16_07_08,LP05a\ROI\LPF\FiltNDown\GLM\corrMap\IOI.mat';

NaClIOImat{1} = 'D:\Edgar\OIS_Results\16_02_25,NC01\ROI\LPF\FiltNDown\GLM\corrMap\IOI.mat';
NaClIOImat{2} = 'D:\Edgar\OIS_Results\16_02_25,NC02\ROI\LPF\FiltNDown\GLM\corrMap\IOI.mat';
NaClIOImat{3} = 'D:\Edgar\OIS_Results\16_02_25,NC03a\ROI\LPF\FiltNDown\GLM\corrMap\IOI.mat';
NaClIOImat{4} = 'D:\Edgar\OIS_Results\16_02_25,NC03b\ROI\LPF\FiltNDown\GLM\corrMap\IOI.mat';
NaClIOImat{5} = 'D:\Edgar\OIS_Results\16_02_26,NC05a\ROI\LPF\FiltNDown\GLM\corrMap\IOI.mat';
NaClIOImat{6} = 'D:\Edgar\OIS_Results\16_02_26,NC05b\ROI\LPF\FiltNDown\GLM\corrMap\IOI.mat';
NaClIOImat{7} = 'D:\Edgar\OIS_Results\16_02_26,NC06a\ROI\LPF\FiltNDown\GLM\corrMap\IOI.mat';
NaClIOImat{8} = 'D:\Edgar\OIS_Results\16_02_26,NC06b\ROI\LPF\FiltNDown\GLM\corrMap\IOI.mat';
NaClIOImat{9} = 'D:\Edgar\OIS_Results\16_07_07,NC07\ROI\LPF\FiltNDown\GLM\corrMap\IOI.mat';
NaClIOImat{10} = 'D:\Edgar\OIS_Results\16_07_08,NC08\ROI\LPF\FiltNDown\GLM\corrMap\IOI.mat';
NaClIOImat{11} = 'D:\Edgar\OIS_Results\16_07_08,NC09\ROI\LPF\FiltNDown\GLM\corrMap\IOI.mat';

nLPS = numel(LPSIOImat);
nNaCl = numel(NaClIOImat);
for iLPS = 1:nLPS,
    load (LPSIOImat{iLPS})
    load(IOI.fcIOS.mask.fnameSeries)
    yGlobalHbO = ioi_alff(brainMaskSeries{1}{5});
    yGlobalHbR = ioi_alff(brainMaskSeries{1}{6});
    load(IOI.ROI.ROIfname)
    for iR = 1:numel(ROI),
        HbOLPS(iLPS, iR) = ioi_alff(ROI{iR}{5}, yGlobalHbO);
        HbRLPS(iLPS, iR) = ioi_alff(ROI{iR}{6}, yGlobalHbR);
    end
end

for iNaCl = 1:nNaCl,
    load (NaClIOImat{iNaCl})
    load(IOI.fcIOS.mask.fnameSeries)
    yGlobalHbO = ioi_alff(brainMaskSeries{1}{5});
    yGlobalHbR = ioi_alff(brainMaskSeries{1}{6});
    load(IOI.ROI.ROIfname)
    for iR = 1:numel(ROI),
        HbONaCl(iNaCl, iR) = ioi_alff(ROI{iR}{5}, yGlobalHbO);
        HbRNaCl(iNaCl, iR) = ioi_alff(ROI{iR}{6}, yGlobalHbR);
    end
end

for iR = 1:numel(ROI),
    [p.HbO(iR), h.HbO(iR)] = ranksum(HbOLPS(:, iR), HbONaCl(:, iR));
    [p.HbR(iR), h.HbR(iR)] = ranksum(HbRLPS(:, iR), HbRNaCl(:, iR));
end

q.HbO = ioi_fdr(p.HbO);
q.HbR = ioi_fdr(p.HbR);

% Save data
save(fullfile(dataFolder,'alff_vals.mat'))
