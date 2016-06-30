%%
IOImatList {1} = 'D:\Edgar\OIS_Results\16_02_25,LP01a\ROI\LPF\IOI.mat';
IOImatList {2} = 'D:\Edgar\OIS_Results\16_02_25,LP01b\ROI\LPF\IOI.mat';
IOImatList {3} = 'D:\Edgar\OIS_Results\16_02_25,NC02\ROI\LPF\IOI.mat';
IOImatList {4} = 'D:\Edgar\OIS_Results\16_02_25,NC03a\ROI\LPF\IOI.mat';
IOImatList {5} = 'D:\Edgar\OIS_Results\16_02_25,NC03b\ROI\LPF\IOI.mat';
IOImatList {6} = 'D:\Edgar\OIS_Results\16_02_26,NC05a\ROI\LPF\IOI.mat';
IOImatList {7} = 'D:\Edgar\OIS_Results\16_02_26,NC05b\ROI\LPF\IOI.mat';
IOImatList {8} = 'D:\Edgar\OIS_Results\16_02_26,NC06a\ROI\LPF\IOI.mat';
IOImatList {9} = 'D:\Edgar\OIS_Results\16_02_26,NC06b\ROI\LPF\IOI.mat';
for iMat = 1:numel(IOImatList)
    load(IOImatList{iMat});
    IOI.res.ROI = ROIbackup;
    save(IOImatList{iMat});
end