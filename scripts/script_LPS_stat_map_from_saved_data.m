%% script_LPS_stat_map
% creates images with hue (correlation map) and alpha (T-stats) color-mapping.
% Contrast loop (5=HbO, 6=HbR)
clear; clc;
for c1 = 5:6,
    % ROI loop (no visual)
    for r1 = 3:10,
        %% Read files
        close all;
        % Load data from one seed, one contrast
        switch(c1)
            case 5
                figFolder = fullfile('C:\Edgar\Dropbox\PostDoc\Newborn\OIS_Results\averaged_maps\HbO\',num2str(r1));
            case 6
                figFolder = fullfile('C:\Edgar\Dropbox\PostDoc\Newborn\OIS_Results\averaged_maps\HbR\',num2str(r1));
        end
        
        %% Load data
        groupToPrintString = 'LPS';
        titleString = sprintf('statMap_%s_C%d_R%02d',groupToPrintString, c1, r1);
        load(fullfile(figFolder,[titleString, '.mat']));
        
        %% Plot maps hue & alpha-coded
        brainMaskAnat = 1.*brainMask;
        brainMaskAnat(~brainMask) = 0.25;
        A_range = [0 10];
        [hFig, hBar] = ioi_dualcodeImage(brainMaskAnat.*imadjust(mat2gray(Underlay)), brainMask.*Bmap_N_S,...
            brainMask.*Tmap_N_S, Pmap_N_S, H_range, A_range);
        
        %% Add seed size & Location
        figure(hFig);
        hold on
        % Display ROI
        rectangle('Position',seedDims,...
            'Curvature',[1,1],...
            'LineWidth',1.5,...
            'LineStyle','-',...
            'EdgeColor','w');
        
        %% Print Statistical Map
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
        disp([titleString ' done!'])
        
    end % ROI loop
end % contrast loop

% EOF