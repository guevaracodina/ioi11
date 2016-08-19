%% Load data
clear; close all; clc
dataFolder = 'D:\Edgar\OIS_Results\extent';
figFolder = 'D:\Edgar\OIS_Results\indStatMaps';

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
alphaVal = 0.05;
% Set the Min/Max values for hue coding
absmax = 1; % Pearson's r
H_range = [-absmax absmax]; % The colormap is symmetric around zero
% Set the Min/Max T-values for alpha coding
A_range = [0 1.5];
% Remove objects with fewer pixels than 5% of the total number of
% suprathreshold voxels (extent threshold)
minPixelsPerc = 5/100;
% Brain Mask
vol = spm_vol('D:\Edgar\OIS_Results\averaged_maps\16_02_25,NC01_anat_brainmask.nii');
brainMaskAll = logical(fix(ioi_MYimresize(spm_read_vols(vol), [512, 512])));
vol = spm_vol('D:\Edgar\OIS_Results\averaged_maps\AVG_Atlas.img');
Underlay = rot90(ioi_MYimresize(spm_read_vols(vol), [512, 512]),3);
        
%%
for c1 = 5:6,                       % Contrast Loop
    %% Group loop
%     for iLPS = [3, 5:7],
    for iNaCl = 9:11,
%         load(LPSIOImat{iLPS})
        load(NaClIOImat{iNaCl})
        load(IOI.fcIOS.corr.fname)
        groupToPrintString = IOI.subj_name;
        vol = spm_vol(IOI.fcIOS.mask.fname);
        brainMaskInd = logical(fix(ioi_MYimresize(spm_read_vols(vol), [512, 512])));
        brainMask = brainMaskAll | brainMaskInd;
%         myIdx = 1;
        for iR = 3:numel(seed_based_fcIOS_map),
%             pMap = seed_based_fcIOS_map{iR}{c1}.pValue;
            zCorrMap = seed_based_fcIOS_map{iR}{c1}.fisher;
            corrMap = seed_based_fcIOS_map{iR}{c1}.pearson;
            % FDR-correction
%             pMask = ~isnan(pMap) & brainMask;
%             pMapFDRtmp = ioi_fdr(pMap(pMask));
%             pMapFDR = nan(size(pMap));
%             pMapFDR(pMask) = pMapFDRtmp;
            
            % Apply threshold
            % Pmap_N_S: 'Binary map indicating significance at P<0.05 (fdr corrected)'
%             Pmap_N_S = (pMapFDR <= alphaVal)  & brainMask;
            Pmap_N_S = 1*brainMask;
            nSignifPixels = nnz(Pmap_N_S);
            ratSignifPixels = nSignifPixels / nnz(brainMask);
%             switch(c1)
%                 case 5
%                     LPSextent.HbO(myIdx, iLPS) = ratSignifPixels;
%                 case 6
%                     LPSextent.HbR(myIdx, iLPS) = ratSignifPixels;
%             end
%             myIdx = myIdx + 1;
            
%             tMap = zeros([size(zCorrMap,1) size(zCorrMap,2)]);
            tMap = zCorrMap;
%             ioi_text_waitbar(0, 'Please wait...');
%             for iRows = 1:size(zCorrMap,1)
%                 for iCols = 1:size(zCorrMap,2)
%                     [~, ~, ~, stats] = ttest(zCorrMap(iRows, iCols), alphaVal);
%                     tMap(iRows, iCols) = stats.tstat;
%                 end
%                 ioi_text_waitbar(iRows/size(zCorrMap,1), sprintf('Processing t-test %d from %d', iRows, size(zCorrMap,1)));
%             end
%             ioi_text_waitbar('Clear');

            %% Prepare data for dual code maps
            %--------------------------------------------------------------------------
            % Bmap_N_S: 'Group-averaged seed-based correlation map'
            Bmap_N_S = rot90(corrMap);
            Bmap_N_S(isnan(Bmap_N_S)) = 0;
            % Tmap_N_S: 'T-statistics for the unpaired t-test of correlation values'
            Tmap_N_S = rot90(tMap);
            % Pmap_N_S: 'Binary map indicating significance at P<0.05 (fdr corrected)'
%             Pmap_N_S = rot90(pMapFDR <= alphaVal)  & brainMask;
            nSignifPixels = nnz(Pmap_N_S);
            % Remove objects with fewer pixels than 5% of the total number of
            % suprathreshold voxels
            minPixels = fix(minPixelsPerc * nSignifPixels);
            % Remove spurious pixels
            Pmap_N_S = bwareaopen(Pmap_N_S, minPixels);
            %--------------------------------------------------------------------------
            
            %% Plot maps hue & alpha-coded
            %--------------------------------------------------------------------------
            % Plot
            brainMaskAnat = 1.*brainMask;
            brainMaskAnat(~brainMask) = 0.25;
            [hFig, hBar] = ioi_dualcodeImage(brainMaskAnat.*imadjust(mat2gray(Underlay)), brainMask.*Bmap_N_S,...
                brainMask.*Tmap_N_S, Pmap_N_S, H_range, A_range);
            %--------------------------------------------------------------------------
            
            %% Add seed size & Location
            load('C:\Edgar\OIS_Results\16_02_25,NC01\ROI\LPF\FiltNDown\GLM\corrMap\IOI.mat')
            seedX = (IOI.res.ROI{iR}.center(2) + IOI.res.ROI{iR}.radius) / IOI.res.shrink_x;
            seedY = (IOI.res.ROI{iR}.center(1) - IOI.res.ROI{iR}.radius) / IOI.res.shrink_y;
            % Seed width
            seedW = 2*IOI.res.ROI{iR}.radius / IOI.res.shrink_x;
            % Seed height
            seedH = 2*IOI.res.ROI{iR}.radius / IOI.res.shrink_y;
            seedDims =  [seedY, size(Bmap_N_S,1) - seedX, seedW, seedH];
            figure(hFig);
            hold on
            % Display ROI
            rectangle('Position',seedDims,...
                'Curvature',[1,1],...
                'LineWidth',1.5,...
                'LineStyle','-',...
                'EdgeColor','w');
            
            %% Print Statistical Map
            titleString = sprintf('indStatMap_%s_C%d_R%02d',groupToPrintString, c1, iR);
            % Specify window units
            set(hFig, 'units', 'inches')
            % Change figure and paper size
            set(hFig, 'Position', [0.1 0.1 3 3])
            set(hFig, 'PaperPosition', [0.1 0.1 3 3])
            % Save as PNG at the user-defined resolution
            print(hFig, '-dpng', ...
                fullfile(figFolder, titleString),...
                sprintf('-r%d', 300));
            close(hFig)
            close(hBar)
            
        end
    end
    

end


% EOF
