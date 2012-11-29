%% Script to plot group correlation tests
function script_group_bilateral_corr(c1)
% Load data
pathName        = 'D:\Edgar\Documents\Dropbox\Batch\IOI\corr_Results_01-07_CMRO2';
fileName        = 'group_corr_pair_seeds.mat';
load(fullfile(pathName, fileName))

% Define contrast
% c1              = 6;
% Define job
job = struct('ttest1', true, 'wilcoxon1', false, 'generate_figures', true,...
    'save_figures', true, 'paired_seeds', [(1:2:12)', (2:2:12)'], 'alpha', 0.05,...
    'parent_results_dir', {{'D:\Edgar\Documents\Dropbox\Docs\Epilepsy\figs'}},...
    'stderror', true);

% Define IOI
IOI.color = struct('eng', 'RGYLODFM', 'red', 'R', 'green', 'G', 'yellow', 'Y',...
    'laser', 'L',  'HbO', 'O',  'HbR', 'D', 'flow', 'F', 'CMRO2', 'M');

[job, IOI, e, y, statTest] = subfunction_group_corr_test(job, IOI, c1, groupCorrData, groupCorrIdx);

subfunction_plot_group_corr_test(job, IOI, c1, e, y, statTest);

end % script_group_bilateral_corr

function [job, IOI, e, y, statTest] = subfunction_group_corr_test(job, IOI, c1, groupCorrData, ~)
% Do a separate paired t-test for each seed data
for iSeeds = 1:size(job.paired_seeds, 1)
    % Average of control group
    meanCorr{iSeeds,c1}(1) = nanmean(groupCorrData{iSeeds,c1}(:,1));
    % Average of treatment group
    meanCorr{iSeeds,c1}(2) = nanmean(groupCorrData{iSeeds,c1}(:,2));
    % Standard deviation of control group
    stdCorr{iSeeds,c1}(1) = nanstd(groupCorrData{iSeeds,c1}(:,1));
    % Standard deviation oftreatment group
    stdCorr{iSeeds,c1}(2) = nanstd(groupCorrData{iSeeds,c1}(:,2));
    
    % Used to plot bar graphs with errorbars
    y(iSeeds,:) = meanCorr{iSeeds,c1};
    e(iSeeds,:) = stdCorr{iSeeds,c1};
    
    % Paired-sample t-test
    if job.ttest1
        [statTest(1).t(1).H{iSeeds,c1}, statTest(1).t(1).P{iSeeds,c1}, ...
            statTest(1).t(1).CI{iSeeds,c1}, statTest(1).t(1).STATS{iSeeds,c1}] ...
            = ttest(...
            groupCorrData{iSeeds,c1}(:,1), groupCorrData{iSeeds,c1}(:,2),...
            job.alpha,'both');
        statTest(1).t(1).id = 'Paired-sample t-test';
    end % t-test
    
    % Wilcoxon rank sum test
    if job.wilcoxon1
        ctrlGroup = groupCorrData{iSeeds,c1}(:,1);
        % ignore NaN values
        ctrlGroup = ctrlGroup(~isnan(ctrlGroup));
        treatmentGroup = groupCorrData{iSeeds,c1}(:,2);
        % ignore NaN values
        treatmentGroup = treatmentGroup(~isnan(treatmentGroup));
        % Perform such test
        [statTest(1).w(1).P{iSeeds,c1}, statTest(1).w(1).H{iSeeds,c1},...
            statTest(1).w(1).STATS{iSeeds,c1}] = ranksum...
            (ctrlGroup, treatmentGroup, 'alpha', job.alpha);
        statTest(1).w(1).id = 'Wilcoxon rank sum test';
    end % Wilcoxon test
    
end % paired-seeds loop
% Show standard error bars instead of standard deviation
if job.stderror
    % std error bars: sigma/sqrt(N)
    e = e / sqrt(size(groupCorrData{iSeeds,c1}, 1));
end

end % subfunction_group_corr_test

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
