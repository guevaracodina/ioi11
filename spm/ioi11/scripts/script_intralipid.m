%% Load IOI matrices
job.IOImat = {
    'E:\Edgar\Data\IOS_Resolution\Results\13_07_30,WL01\GLMfcIOS\corrMap\IOI.mat'
    'E:\Edgar\Data\IOS_Resolution\Results\13_07_30,WL02\GLMfcIOS\corrMap\IOI.mat'
    'E:\Edgar\Data\IOS_Resolution\Results\13_07_30,WL03\GLMfcIOS\corrMap\IOI.mat'
    'E:\Edgar\Data\IOS_Resolution\Results\13_07_30,WL04\GLMfcIOS\corrMap\IOI.mat'
    'E:\Edgar\Data\IOS_Resolution\Results\13_07_30,WL05\GLMfcIOS\corrMap\IOI.mat'
    'E:\Edgar\Data\IOS_Resolution\Results\13_07_30,WL06\GLMfcIOS\corrMap\IOI.mat'
    'E:\Edgar\Data\IOS_Resolution\Results\13_07_30,WL07\GLMfcIOS\corrMap\IOI.mat'
    'E:\Edgar\Data\IOS_Resolution\Results\13_07_30,WL08\GLMfcIOS\corrMap\IOI.mat'
    'E:\Edgar\Data\IOS_Resolution\Results\13_07_30,WL09\GLMfcIOS\corrMap\IOI.mat'
    'E:\Edgar\Data\IOS_Resolution\Results\13_07_30,WL10\GLMfcIOS\corrMap\IOI.mat'
	};
job.IOImatCopyChoice = [];
% CBF only
c1 = 7;
nSubj = numel(job.IOImat);
flowCorrMat = zeros([12 12 nSubj]);

%% Compute average correlation matrix
for SubjIdx = 1:nSubj
    [IOI IOImat dir_ioimat] = ioi_get_IOI(job,SubjIdx);
    matStruct = load(IOI.fcIOS.corr.corrMatrixFname);
    flowCorrMat(:,:,SubjIdx) = matStruct.seed2seedCorrMat{1}{c1};
end
flowCorrMatAvg = mean(flowCorrMat,3);

%% Display average correlation matrix
h = figure;
set(h,'color','w');
colormap(get_colormaps('rwbdoppler'));
subplot(121)
imagesc(flowCorrMatAvg, [-1 1]);
title('r')
colorbar
axis image
subplot(122)
imagesc(fisherz(flowCorrMatAvg), [-1 1]);
title('z(r)')
colorbar
axis image
% EOF
