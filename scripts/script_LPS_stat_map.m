%% script_LPS_stat_map
clear; close all; clc
c1 = 6;
r1 = 5;
alphaVal = 0.05;
alphaValFDR = 0.0001;

%% Read files
% Get anatomical image
% anatVol             = spm_vol(fullfile('C:\Edgar\Dropbox\PostDoc\Newborn\OIS_Results\averaged_maps\',...
%     'AVG_Atlas.img'));
% anatomical          = spm_read_vols(anatVol);

% Read brain mask
% maskVol             = spm_vol(fullfile('C:\Edgar\Dropbox\PostDoc\Newborn\OIS_Results\averaged_maps\',...
%     '16_02_25,NC01_anat_brainmask.nii'));
% brainMask           = logical(spm_read_vols(maskVol));

% Load data from 1 seed, 1 contrast
switch(c1)
    case 5
        figFolder = fullfile('C:\Edgar\Dropbox\PostDoc\Newborn\OIS_Results\averaged_maps\HbO\',num2str(r1));
    case 6
        figFolder = fullfile('C:\Edgar\Dropbox\PostDoc\Newborn\OIS_Results\averaged_maps\HbR\',num2str(r1));
end
job.parent_results_dir{1} = fullfile(figFolder,'overlay');
load(fullfile(figFolder, sprintf('stats_R%d_C%d.mat', r1, c1)))
% load('C:\Edgar\Dropbox\PostDoc\Newborn\OIS_Results\averaged_maps\HbR\7\stats_R7_C6.mat')

% Get functional images from ctrl group
ctrl = LPS;

% Mask out non-brain elements
ctrl(~repmat(brainMask,[1 1 size(ctrl,3)])) =  nan;
% Mask out non-brain elements
LPS(~repmat(brainMask,[1 1 size(LPS,3)])) =  nan;
anatomical = Underlay;

%% Perform statistical test (1st level analysis)
[ tMap, pMapFDR, pMapFDRalpha ] = ioi_stat_map( ctrl, alphaVal );

%% 1. Load data
%--------------------------------------------------------------------------
% load AOD_data.mat
close all; clc
% addpath(genpath('C:\Edgar\Dropbox\Matlab\dualcode'))
% For a single axial slice (Z = 2 mm) of data, you should have:
% Bmap_N_S: 'Difference between Novel and Standard betas averaged over 28 subjects'
Bmap_N_S = rot90(squeeze(nanmean(ctrl, 3)));
Bmap_N_S(isnan(Bmap_N_S)) = 0;
% Tmap_N_S: 'T-statistics for the paired t-test comparing Novel and Standard betas'
Tmap_N_S = rot90(tMap);
% Pmap_N_S: 'Binary map indicating significance at P<0.001 (fdr corrected)'
Pmap_N_S = rot90(pMapFDR <= alphaValFDR)  & brainMask;
% Underlay: 'Structural image ch2bet from MRIcron, warped to functional data'
%--------------------------------------------------------------------------

%% 2. Plot maps hue & alpha-coded
%--------------------------------------------------------------------------
% Set the Min/Max values for hue coding
% absmax = max(abs(Bmap_N_S(:))); 
absmax = 1; % Pearson's r
H_range = [-absmax absmax]; % The colormap is symmetric around zero

% Set the Min/Max T-values for alpha coding
A_range = [0 15];
ioi_dualcodeImage (imadjust(mat2gray(Underlay)), Bmap_N_S, Tmap_N_S, Pmap_N_S, H_range, A_range)

%% Create nifti
% hdr = ioi_create_vol(fullfile('F:\Edgar\Data\PAT_Results_20130517\alignment\',...
%     'ROI05_SO2_pMap_alpha.nii'), anatVol.dim, [64 0], anatVol.pinfo, anatVol.mat, 1, pMapAlpha);
% hdr = ioi_create_vol(fullfile('F:\Edgar\Data\PAT_Results_20130517\alignment\',...
%     'ROI05_SO2_pMap_alpha_FDR.nii'), anatVol.dim, [64 0], anatVol.pinfo, anatVol.mat, 1, pMapFDRalpha);

% EOF
