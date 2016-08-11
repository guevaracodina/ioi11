%% script_LPS_stat_map
% creates images with hue (correlation map) and alpha (T-stats) color-mapping.
% Contrast loop (5=HbO, 6=HbR)
for c1 = 6,
    % ROI loop (no visual)
    for r1 = 5,
        close all;
        % Significance (height) threshold
        alphaVal = 0.05;
        % Significance (height) threshold after FDR correction (arbitrary???)
        alphaValFDR = 0.0001;
        % Set the Min/Max values for hue coding
        absmax = 1; % Pearson's r
        H_range = [-absmax absmax]; % The colormap is symmetric around zero
        % Set the Min/Max T-values for alpha coding
        A_range = [0 15];
        % Remove objects with fewer pixels than 5% of the total number of
        % suprathreshold voxels (extent threshold)
        minPixelsPerc = 5/100;
        % Function to convert a variable name to string
        getVarName=@(x) inputname(1);
        
        %% Read files
        % Load data from one seed, one contrast
        switch(c1)
            case 5
                figFolder = fullfile('C:\Edgar\Dropbox\PostDoc\Newborn\OIS_Results\averaged_maps\HbO\',num2str(r1));
            case 6
                figFolder = fullfile('C:\Edgar\Dropbox\PostDoc\Newborn\OIS_Results\averaged_maps\HbR\',num2str(r1));
        end
        % job.parent_results_dir{1} = fullfile(figFolder,'overlay');
        load(fullfile(figFolder, sprintf('stats_R%d_C%d.mat', r1, c1)))
        
        % clean up
        clear job hSpatial pSpatial LPS_spatial_extension NaCl_spatial_extension ...
            LPSimages NaClimages
        
        % Get functional images from NaCl group after loading MAT files
%         groupToPrint = NaCl;
%         groupToPrintString = getVarName(NaCl);
        
        % Get functional images from LPS group after loading MAT files
        groupToPrint = LPS;
        groupToPrintString = getVarName(LPS);
        
        % Mask out non-brain elements
        groupToPrint(~repmat(brainMask,[1 1 size(groupToPrint,3)])) =  nan;
                
        %% Perform statistical test (2nd level analysis, i.e. look at  group)
        [tMap, pMapFDR, pMapFDRalpha] = ioi_stat_map(groupToPrint, alphaVal);
        
        %% Prepare data for dual code maps
        %--------------------------------------------------------------------------
        % Bmap_N_S: 'Group-averaged seed-based correlation map'
        Bmap_N_S = rot90(squeeze(nanmean(groupToPrint, 3)));
        Bmap_N_S(isnan(Bmap_N_S)) = 0;
        % Tmap_N_S: 'T-statistics for the unpaired t-test of correlation values'
        Tmap_N_S = rot90(tMap);
        % Pmap_N_S: 'Binary map indicating significance at P<0.05 (fdr corrected)'
        Pmap_N_S = rot90(pMapFDR <= alphaVal)  & brainMask;
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
        load('C:\Edgar\Dropbox\PostDoc\Newborn\OIS_Results\16_02_25,NC01\ROI\LPF\FiltNDown\GLM\corrMap\IOI.mat')
        seedX = (IOI.res.ROI{r1}.center(2) + IOI.res.ROI{r1}.radius) / IOI.res.shrink_x;
        seedY = (IOI.res.ROI{r1}.center(1) - IOI.res.ROI{r1}.radius) / IOI.res.shrink_y;
        % Seed width
        seedW = 2*IOI.res.ROI{r1}.radius / IOI.res.shrink_x;
        % Seed height
        seedH = 2*IOI.res.ROI{r1}.radius / IOI.res.shrink_y;
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
        titleString = sprintf('statMap_%s_C%d_R%02d',groupToPrintString, c1, r1);
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
        
        %% Save data
        save(fullfile(figFolder,[titleString, '.mat']));
        disp([titleString ' done!'])
    end % ROI loop
end % contrast loop

% EOF