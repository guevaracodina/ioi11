%% Script to plot group correlation tests
function script_group_bilateral_corr(c1, dType)
switch dType
    case 'filt'
        fileName        = 'group_corr_pair_seeds.mat';
    case 'diff'
        fileName        = 'group_corr_pair_seeds_diff.mat';
    case 'raw'
        fileName        = 'group_corr_pair_seeds_raw.mat';
    otherwise
        fileName        = 'group_corr_pair_seeds.mat';
end

% Load data
pathName        = 'D:\Edgar\Documents\Dropbox\Batch\IOI\corr_Results_01-07_CMRO2';
load(fullfile(pathName, fileName))

% Define job
job = struct('ttest1', true, 'wilcoxon1', false, 'generate_figures', true,...
    'save_figures', true, 'paired_seeds', [(1:2:12)', (2:2:12)'], 'alpha', 0.05,...
    'parent_results_dir', {{'D:\Edgar\Documents\Dropbox\Docs\Epilepsy\figs'}},...
    'stderror', true);

% Define IOI
IOI.color = struct('eng', 'RGYLODFM', 'red', 'R', 'green', 'G', 'yellow', 'Y',...
    'laser', 'L',  'HbO', 'O',  'HbR', 'D', 'flow', 'F', 'CMRO2', 'M');

switch dType
    case 'filt'
        subfunction_plot_group_corr_test(job, IOI, c1, eTotal{c1}, yTotal{c1}, statTest);
    case 'diff'
        subfunction_plot_group_corr_test_diff(job, IOI, c1, eTotalDiff{c1}, yTotalDiff{c1}, statTestDiff)
    case 'raw'
        subfunction_plot_group_corr_test_raw(job, IOI, c1, eTotalRaw{c1}, yTotalRaw{c1}, statTestRaw)
    otherwise
        subfunction_plot_group_corr_test(job, IOI, c1, eTotal{c1}, yTotal{c1}, statTest);
end


end % script_group_bilateral_corr


function subfunction_plot_group_corr_test(job, IOI, c1, e, y, statTest)
% Plots statistical analysis group results
colorNames      = fieldnames(IOI.color);
% Positioning factor for the * mark, depends on max data value at the given seed
starPosFactor   = 1.05;
% Font Sizes
titleFontSize   = 12;
axisFontSize    = 20;
labelFontSize   = 24;
legendFontSize  = 20;
starFontSize    = 36;
yLimits         = [-0.9 1.3];
yTicks          = -1:0.5:1;

if job.ttest1
    % Display a graph with ROI labels
    if job.generate_figures
        % Display plots on SPM graphics window
%         h = spm_figure('GetWin', 'Graphics');
%         spm_figure('Clear', 'Graphics');
        h = figure; set(h,'color','w')
        % Custom bar graphs with error bars (1st arg: error)
        barwitherr(e, y)
        
        % Display colormap according to the contrast
        switch(c1)
            case 5
                % HbO contrast
                colormap([1 0 0; 1 1 1]);
                legendLocation = 'northwest';
            case 6
                % HbR contrast
                colormap([0 0 1; 1 1 1]);
                legendLocation = 'southwest';
            case 7
                % Flow contrast
                colormap([0.5 0.5 0.5; 1 1 1]);
                legendLocation = 'southwest';
            case 8
                % CMRO2 contrast
                colormap([0 0 0; 1 1 1]);
                legendLocation = 'southwest';
            otherwise
                colormap(gray)
        end
%         title(sprintf('Bilateral Correlation before/after 4-AP C%d(%s). T-test (*p<%.2g)',...
%             c1,colorNames{1+c1},job.alpha),'interpreter','none','FontSize',titleFontSize)
        set(gca,'FontSize',axisFontSize)
        ylabel('Functional correlation z(r)','FontSize',labelFontSize)
        set(gca,'XTickLabel',{'F', 'M', 'C', 'S', 'R', 'V'},'FontWeight', 'b','FontSize',labelFontSize)
        legend({'Control' '4-AP'},'FontSize',legendFontSize, 'location', legendLocation)
        axis tight
        set(gca, 'ylim', yLimits)
        set(gca, 'ytick', yTicks)
        
        % Show a * when a significant difference is found.
        for iSeeds = 1:size(job.paired_seeds, 1)
            if statTest(1).t(1).H{iSeeds,c1}
                if max(y(iSeeds,:))>=0
                    yPos = starPosFactor*(max(y(iSeeds,:)) + max(e(iSeeds,:)));
                else
                    yPos = starPosFactor*(min(y(iSeeds,:)) - max(e(iSeeds,:)));
                end
                xPos = iSeeds;
                text(xPos, yPos, '*', 'FontSize', starFontSize, 'FontWeight', 'b');
            end
        end
        if job.save_figures
            newName = sprintf('groupCorr_Ttest_C%d_(%s)',c1,colorNames{1+c1});
            % Save as EPS
            spm_figure('Print', 'Graphics', fullfile(job.parent_results_dir{1}, newName));
            % Save as PNG
            % print(h, '-dpng', fullfile(job.parent_results_dir{1},newName), '-r300');
            addpath(genpath('D:\Edgar\ssoct\Matlab'))
            export_fig(fullfile(job.parent_results_dir{1},newName),'-png',gcf)
            % Save as a figure
            saveas(h, fullfile(job.parent_results_dir{1},newName), 'fig');
        end
    end % end generate figures
end

if job.wilcoxon1
    % Display a graph with ROI labels
    if job.generate_figures
         % Display plots on SPM graphics window
%         h = spm_figure('GetWin', 'Graphics');
%         spm_figure('Clear', 'Graphics');
        h = figure; set(h,'color','w')
        % Custom bar graphs with error bars (1st arg: error)
        barwitherr(e, y)
        % Display colormap according to the contrast
        switch(c1)
            case 5
                % HbO contrast
                colormap([1 0 0; 1 1 1]);
                legendLocation = 'northwest';
            case 6
                % HbR contrast
                colormap([0 0 1; 1 1 1]);
                legendLocation = 'southwest';
            case 7
                % Flow contrast
                colormap([0.5 0.5 0.5; 1 1 1]);
                legendLocation = 'southwest';
            case 8
                % CMRO2 contrast
                colormap([0 0 0; 1 1 1]);
                legendLocation = 'southwest';
            otherwise
                colormap(gray)
        end
%         title(sprintf('Bilateral Correlation before/after 4-AP C%d(%s). Wilcoxon *(p<%.2g)',...
%             c1,colorNames{1+c1},job.alpha),'interpreter','none','FontSize',titleFontSize)
        set(gca,'FontSize',axisFontSize)
        ylabel('Functional correlation z(r)','FontSize',labelFontSize)
        set(gca,'XTickLabel',{'F', 'M', 'C', 'S', 'R', 'V'},'FontWeight', 'b','FontSize',labelFontSize)
        legend({'Control' '4-AP'},'FontSize',legendFontSize, 'location', legendLocation)
        axis tight
        set(gca, 'ylim', yLimits)
        set(gca, 'ytick', yTicks)
        
        % Show a * when a significant difference is found.
        for iSeeds = 1:size(job.paired_seeds, 1)
            if statTest(1).w(1).H{iSeeds,c1}
                if max(y(iSeeds,:))>=0
                    yPos = starPosFactor*(max(y(iSeeds,:)) + max(e(iSeeds,:)));
                else
                    yPos = starPosFactor*(min(y(iSeeds,:)) - max(e(iSeeds,:)));
                end
                xPos = iSeeds;
                text(xPos, yPos,'*', 'FontSize', starFontSize, 'FontWeight', 'b');
            end
        end
        if job.save_figures
            newName = sprintf('groupCorr_Wtest_C%d_(%s)',c1,colorNames{1+c1});
            % Save as EPS
            spm_figure('Print', 'Graphics', fullfile(job.parent_results_dir{1},newName));
            % Save as PNG
%             print(h, '-dpng', fullfile(job.parent_results_dir{1},newName), '-r300');
            addpath(genpath('D:\Edgar\ssoct\Matlab'))
            export_fig(fullfile(job.parent_results_dir{1},newName),'-png',gcf)
            % Save as a figure
            saveas(h, fullfile(job.parent_results_dir{1},newName), 'fig');
        end
    end % End generate figures
end

end % subfunction_plot_group_corr_test

function subfunction_plot_group_corr_test_diff(job, IOI, c1, eDiff, yDiff, statTestDiff)
% Plots statistical analysis group results
colorNames      = fieldnames(IOI.color);
% Positioning factor for the * mark, depends on max data value at the given seed
starPosFactor   = 1.05;
% Font Sizes
titleFontSize   = 12;
axisFontSize    = 20;
labelFontSize   = 24;
legendFontSize  = 20;
starFontSize    = 36;
yLimits         = [-0.9 1.3];
yTicks          = -1:0.5:1;

if job.ttest1
    % Display a graph with ROI labels
    if job.generate_figures
        % Display plots on SPM graphics window
%         h = spm_figure('GetWin', 'Graphics');
%         spm_figure('Clear', 'Graphics');
        h = figure; set(h,'color','w')
        % Custom bar graphs with error bars (1st arg: error)
        barwitherr(eDiff, yDiff)
        
        % Display colormap according to the contrast
        switch(c1)
            case 5
                % HbO contrast
                colormap([1 0 0; 1 1 1]);
                legendLocation = 'northwest';
            case 6
                % HbR contrast
                colormap([0 0 1; 1 1 1]);
                legendLocation = 'southwest';
            case 7
                % Flow contrast
                colormap([0.5 0.5 0.5; 1 1 1]);
                legendLocation = 'southwest';
            case 8
                % CMRO2 contrast
                colormap([0 0 0; 1 1 1]);
                legendLocation = 'southwest';
            otherwise
                colormap(gray)
        end
%         title(sprintf('Bilateral Correlation before/after 4-AP C%d(%s). T-test (*p<%.2g)',...
%             c1,colorNames{1+c1},job.alpha),'interpreter','none','FontSize',titleFontSize)
        set(gca,'FontSize',axisFontSize)
        ylabel('Functional correlation z(r)','FontSize',labelFontSize)
        set(gca,'XTickLabel',{'F', 'M', 'C', 'S', 'R', 'V'},'FontWeight', 'b','FontSize',labelFontSize)
        legend({'Control' '4-AP'},'FontSize',legendFontSize, 'location', legendLocation)
        axis tight
        set(gca, 'ylim', yLimits)
        set(gca, 'ytick', yTicks)
        
        % Show a * when a significant difference is found.
        for iSeeds = 1:size(job.paired_seeds, 1)
            if statTestDiff(1).t(1).H{iSeeds,c1}
                if max(yDiff(iSeeds,:))>=0
                    yPos = starPosFactor*(max(yDiff(iSeeds,:)) + max(eDiff(iSeeds,:)));
                else
                    yPos = starPosFactor*(min(yDiff(iSeeds,:)) - max(eDiff(iSeeds,:)));
                end
                xPos = iSeeds;
                text(xPos, yPos, '*', 'FontSize', starFontSize, 'FontWeight', 'b');
            end
        end
        if job.save_figures
            newName = sprintf('groupCorr_Ttest_C%d_(%s)_Diff',c1,colorNames{1+c1});
            % Save as EPS
            spm_figure('Print', 'Graphics', fullfile(job.parent_results_dir{1}, newName));
            % Save as PNG
            % print(h, '-dpng', fullfile(job.parent_results_dir{1},newName), '-r300');
            addpath(genpath('D:\Edgar\ssoct\Matlab'))
            export_fig(fullfile(job.parent_results_dir{1},newName),'-png',gcf)
            % Save as a figure
            saveas(h, fullfile(job.parent_results_dir{1},newName), 'fig');
        end
    end % end generate figures
end

if job.wilcoxon1
    % Display a graph with ROI labels
    if job.generate_figures
         % Display plots on SPM graphics window
%         h = spm_figure('GetWin', 'Graphics');
%         spm_figure('Clear', 'Graphics');
        h = figure; set(h,'color','w')
        % Custom bar graphs with error bars (1st arg: error)
        barwitherr(eDiff, yDiff)
        % Display colormap according to the contrast
        switch(c1)
            case 5
                % HbO contrast
                colormap([1 0 0; 1 1 1]);
                legendLocation = 'northwest';
            case 6
                % HbR contrast
                colormap([0 0 1; 1 1 1]);
                legendLocation = 'southwest';
            case 7
                % Flow contrast
                colormap([0.5 0.5 0.5; 1 1 1]);
                legendLocation = 'southwest';
            case 8
                % CMRO2 contrast
                colormap([0 0 0; 1 1 1]);
                legendLocation = 'southwest';
            otherwise
                colormap(gray)
        end
%         title(sprintf('Bilateral Correlation before/after 4-AP C%d(%s). Wilcoxon *(p<%.2g)',...
%             c1,colorNames{1+c1},job.alpha),'interpreter','none','FontSize',titleFontSize)
        set(gca,'FontSize',axisFontSize)
        ylabel('Functional correlation z(r)','FontSize',labelFontSize)
        set(gca,'XTickLabel',{'F', 'M', 'C', 'S', 'R', 'V'},'FontWeight', 'b','FontSize',labelFontSize)
        legend({'Control' '4-AP'},'FontSize',legendFontSize, 'location', legendLocation)
        axis tight
        set(gca, 'ylim', yLimits)
        set(gca, 'ytick', yTicks)
        
        % Show a * when a significant difference is found.
        for iSeeds = 1:size(job.paired_seeds, 1)
            if statTestDiff(1).w(1).H{iSeeds,c1}
                if max(yDiff(iSeeds,:))>=0
                    yPos = starPosFactor*(max(yDiff(iSeeds,:)) + max(eDiff(iSeeds,:)));
                else
                    yPos = starPosFactor*(min(yDiff(iSeeds,:)) - max(eDiff(iSeeds,:)));
                end
                xPos = iSeeds;
                text(xPos, yPos,'*', 'FontSize', starFontSize, 'FontWeight', 'b');
            end
        end
        if job.save_figures
            newName = sprintf('groupCorr_Wtest_C%d_(%s)_Diff',c1,colorNames{1+c1});
            % Save as EPS
            spm_figure('Print', 'Graphics', fullfile(job.parent_results_dir{1},newName));
            % Save as PNG
%             print(h, '-dpng', fullfile(job.parent_results_dir{1},newName), '-r300');
            addpath(genpath('D:\Edgar\ssoct\Matlab'))
            export_fig(fullfile(job.parent_results_dir{1},newName),'-png',gcf)
            % Save as a figure
            saveas(h, fullfile(job.parent_results_dir{1},newName), 'fig');
        end
    end % End generate figures
end

end % subfunction_plot_group_corr_test_diff

function subfunction_plot_group_corr_test_raw(job, IOI, c1, eRaw, yRaw, statTestRaw)
% Plots statistical analysis group results
colorNames      = fieldnames(IOI.color);
% Positioning factor for the * mark, depends on max data value at the given seed
starPosFactor   = 1.05;
% Font Sizes
titleFontSize   = 12;
axisFontSize    = 20;
labelFontSize   = 24;
legendFontSize  = 20;
starFontSize    = 36;
yLimits         = [-0.2 3];
yTicks          = 0:0.5:2.5;

if job.ttest1
    % Display a graph with ROI labels
    if job.generate_figures
        % Display plots on SPM graphics window
%         h = spm_figure('GetWin', 'Graphics');
%         spm_figure('Clear', 'Graphics');
        h = figure; set(h,'color','w')
        % Custom bar graphs with error bars (1st arg: error)
        barwitherr(eRaw, yRaw)
        
        % Display colormap according to the contrast
        switch(c1)
            case 5
                % HbO contrast
                colormap([1 0 0; 1 1 1]);
                legendLocation = 'northwest';
            case 6
                % HbR contrast
                colormap([0 0 1; 1 1 1]);
                legendLocation = 'northwest';
            case 7
                % Flow contrast
                colormap([0.5 0.5 0.5; 1 1 1]);
                legendLocation = 'northwest';
            case 8
                % CMRO2 contrast
                colormap([0 0 0; 1 1 1]);
                legendLocation = 'northwest';
            otherwise
                colormap(gray)
        end
%         title(sprintf('Bilateral Correlation before/after 4-AP C%d(%s). T-test (*p<%.2g)',...
%             c1,colorNames{1+c1},job.alpha),'interpreter','none','FontSize',titleFontSize)
        set(gca,'FontSize',axisFontSize)
        ylabel('Functional correlation z(r)','FontSize',labelFontSize)
        set(gca,'XTickLabel',{'F', 'M', 'C', 'S', 'R', 'V'},'FontWeight', 'b','FontSize',labelFontSize)
        legend({'Control' '4-AP'},'FontSize',legendFontSize, 'location', legendLocation)
        axis tight
        set(gca, 'ylim', yLimits)
        set(gca, 'ytick', yTicks)
        
        % Show a * when a significant difference is found.
        for iSeeds = 1:size(job.paired_seeds, 1)
            if statTestRaw(1).t(1).H{iSeeds,c1}
                if max(yRaw(iSeeds,:))>=0
                    yPos = starPosFactor*(max(yRaw(iSeeds,:)) + max(eRaw(iSeeds,:)));
                else
                    yPos = starPosFactor*(min(yRaw(iSeeds,:)) - max(eRaw(iSeeds,:)));
                end
                xPos = iSeeds;
                text(xPos, yPos, '*', 'FontSize', starFontSize, 'FontWeight', 'b');
            end
        end
        if job.save_figures
            newName = sprintf('groupCorr_Ttest_C%d_(%s)_Raw',c1,colorNames{1+c1});
            % Save as EPS
            spm_figure('Print', 'Graphics', fullfile(job.parent_results_dir{1}, newName));
            % Save as PNG
            % print(h, '-dpng', fullfile(job.parent_results_dir{1},newName), '-r300');
            addpath(genpath('D:\Edgar\ssoct\Matlab'))
            export_fig(fullfile(job.parent_results_dir{1},newName),'-png',gcf)
            % Save as a figure
            saveas(h, fullfile(job.parent_results_dir{1},newName), 'fig');
        end
    end % end generate figures
end

if job.wilcoxon1
    % Display a graph with ROI labels
    if job.generate_figures
         % Display plots on SPM graphics window
%         h = spm_figure('GetWin', 'Graphics');
%         spm_figure('Clear', 'Graphics');
        h = figure; set(h,'color','w')
        % Custom bar graphs with error bars (1st arg: error)
        barwitherr(eRaw, yRaw)
        % Display colormap according to the contrast
        switch(c1)
            case 5
                % HbO contrast
                colormap([1 0 0; 1 1 1]);
                legendLocation = 'northwest';
            case 6
                % HbR contrast
                colormap([0 0 1; 1 1 1]);
                legendLocation = 'southwest';
            case 7
                % Flow contrast
                colormap([0.5 0.5 0.5; 1 1 1]);
                legendLocation = 'southwest';
            case 8
                % CMRO2 contrast
                colormap([0 0 0; 1 1 1]);
                legendLocation = 'southwest';
            otherwise
                colormap(gray)
        end
%         title(sprintf('Bilateral Correlation before/after 4-AP C%d(%s). Wilcoxon *(p<%.2g)',...
%             c1,colorNames{1+c1},job.alpha),'interpreter','none','FontSize',titleFontSize)
        set(gca,'FontSize',axisFontSize)
        ylabel('Functional correlation z(r)','FontSize',labelFontSize)
        set(gca,'XTickLabel',{'F', 'M', 'C', 'S', 'R', 'V'},'FontWeight', 'b','FontSize',labelFontSize)
        legend({'Control' '4-AP'},'FontSize',legendFontSize, 'location', legendLocation)
        axis tight
        set(gca, 'ylim', yLimits)
        set(gca, 'ytick', yTicks)
        
        % Show a * when a significant difference is found.
        for iSeeds = 1:size(job.paired_seeds, 1)
            if statTestRaw(1).w(1).H{iSeeds,c1}
                if max(yRaw(iSeeds,:))>=0
                    yPos = starPosFactor*(max(yRaw(iSeeds,:)) + max(eRaw(iSeeds,:)));
                else
                    yPos = starPosFactor*(min(yRaw(iSeeds,:)) - max(eRaw(iSeeds,:)));
                end
                xPos = iSeeds;
                text(xPos, yPos,'*', 'FontSize', starFontSize, 'FontWeight', 'b');
            end
        end
        if job.save_figures
            newName = sprintf('groupCorr_Wtest_C%d_(%s)_Raw',c1,colorNames{1+c1});
            % Save as EPS
            spm_figure('Print', 'Graphics', fullfile(job.parent_results_dir{1},newName));
            % Save as PNG
%             print(h, '-dpng', fullfile(job.parent_results_dir{1},newName), '-r300');
            addpath(genpath('D:\Edgar\ssoct\Matlab'))
            export_fig(fullfile(job.parent_results_dir{1},newName),'-png',gcf)
            % Save as a figure
            saveas(h, fullfile(job.parent_results_dir{1},newName), 'fig');
        end
    end % End generate figures
end
end % subfunction_plot_group_corr_test_raw
