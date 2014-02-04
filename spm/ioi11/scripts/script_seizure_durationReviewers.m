%% script_seizure_duration
load('F:\Edgar\Data\IOS_Results\corr_Results_seizuresReviewers\group_corr_pair_seeds.mat')
load('F:\Edgar\Data\IOS_Results\corr_Results_seizuresReviewers\group_corr_pair_seeds_diff.mat')

%% Fix correlation data on EG03 and EG04
% CBF
c1 = 7;
load('F:\Edgar\Data\IOS_Results\12_08_24,EG03\FiltNDown\GLMfcIOS\corrMap\IOI.mat')
load(IOI.fcIOS.filtNdown.fname)
CBFrS02 = [];
% Loop over sessions
for s1=1:length(IOI.sess_res)
    pair = 1;
    % Loop over ROIs
    for r1 = 1:2:numel(IOI.ROIname)
        CBFrS02(s1, pair) = fisherz(corr(filtNdownROI{r1}{s1,c1}', filtNdownROI{r1+1}{s1,c1}'));
        pair = pair + 1;
    end
end

load('F:\Edgar\Data\IOS_Results\12_08_24,EG04\FiltNDown\IOI.mat')
load(IOI.fcIOS.filtNdown.fname)
CBFrS03 = [];
% Loop over sessions
for s1=1:length(IOI.sess_res)
    pair = 1;
    % Loop over ROIs
    for r1 = 1:2:numel(IOI.ROIname)
        CBFrS03(s1, pair) = fisherz(corr(filtNdownROI{r1}{s1,c1}', filtNdownROI{r1+1}{s1,c1}'));
        pair = pair + 1;
    end
end

% Fix correlation values
for pair = 1:size(groupCorrData,1)-1
    groupCorrData{pair,c1}(7,1) = CBFrS02(1,pair);
    groupCorrData{pair,c1}(7,2) = CBFrS02(3,pair);
    groupCorrData{pair,c1}(8,2) = CBFrS02(3,pair);
    groupCorrData{pair,c1}(9,1) = CBFrS03(1,pair);
    groupCorrData{pair,c1}(9,2) = CBFrS03(3,pair);
    groupCorrData{pair,c1}(10,1) = CBFrS03(2,pair);
    groupCorrData{pair,c1}(10,2) = CBFrS03(3,pair);
    groupCorrData{pair,c1}(11,1) = CBFrS03(1,pair);
    groupCorrData{pair,c1}(11,2) = CBFrS03(4,pair);
    groupCorrData{pair,c1}(12,1) = CBFrS03(2,pair);
    groupCorrData{pair,c1}(12,2) = CBFrS03(4,pair);
end
% Fix correlation values for 6th seeds pair
groupCorrData{6,c1}(9,1) = CBFrS03(1,pair);
groupCorrData{6,c1}(10,1) = CBFrS03(2,pair);
groupCorrData{6,c1}(11,1) = CBFrS03(1,pair);
groupCorrData{6,c1}(11,2) = CBFrS03(4,pair);
groupCorrData{6,c1}(12,1) = CBFrS03(2,pair);
groupCorrData{6,c1}(12,2) = CBFrS03(4,pair);

%% Average seizure duration 
% Preallocate and gather data for the sessions after the 4AP injection
seizureDuration = num2cell(nan(max(groupCorrIdx{1,5}(:,3)), 5));
% seizureDuration{iSubject, iSession}
seizureDuration{1, 4} = mean([203-158; 375-332; 660-627; 838-806]);
seizureDuration{1, 5} = mean([244-209; 433-401; 729-702]);

% Only 1 4AP session
seizureDuration{3-1, 3} = mean([81-54; 332-299; 852-818]);

seizureDuration{4-1, 3} = mean([334-220]);
seizureDuration{4-1, 4} = mean([392-321; 812-700]);


seizureDuration{7-3, 3} = mean([200-162; 384-354; 748-718]);
seizureDuration{7-3, 4} = mean([382-349; 566-534]);

seizureDuration{9-4, 3} = mean([459-329]);
seizureDuration{9-4, 4} = mean([432-271]);

%% Seizure percentage
sessionTotalDuration = 863;

% seizurePercent{iSubject, iSession}

seizurePercent{1, 4} = 100*([203-158; 375-332; 660-627; 838-806])/sessionTotalDuration;
seizurePercent{1, 5} = 100*([244-209; 433-401; 729-702])/sessionTotalDuration;

% Only 1 4AP session
seizurePercent{3-1, 3} = 100*([81-54; 332-299; 852-818])/sessionTotalDuration;

seizurePercent{4-1, 3} = 100*([334-220])/sessionTotalDuration;
seizurePercent{4-1, 4} = 100*([392-321; 812-700])/sessionTotalDuration;

seizurePercent{7-3, 3} = 100*([200-162; 384-354; 748-718])/sessionTotalDuration;
seizurePercent{7-3, 4} = 100*([382-349; 566-534])/sessionTotalDuration;

seizurePercent{9-4, 3} = 100*([459-329])/sessionTotalDuration;
seizurePercent{9-4, 4} = 100*([432-271])/sessionTotalDuration;

%% Average number of seizures and duration
% nSeizureAvg{iSubject, iSession}

nSeizureAvg{1, 4} = mean([numel([203-158; 375-332; 660-627; 838-806]) numel([244-209; 433-401; 729-702])]);

% Only 1 4AP session
nSeizureAvg{3-1, 3} = mean(numel([81-54; 332-299; 852-818]));

nSeizureAvg{4-1, 3} = mean([numel([334-220]) numel([392-321; 812-700])]);

nSeizureAvg{7-3, 3} = mean([numel([200-162; 384-354; 748-718]) numel([382-349; 566-534])]);

nSeizureAvg{9-4, 3} = mean([numel([459-329]) numel([432-271])]);


%% Only subjects where seizures where observed
seizureDurationTrans = seizureDuration';
seizDurationVector = cell2mat(seizureDurationTrans(~isnan(cell2mat(seizureDurationTrans))));

% Time percent of the session, marked as seizure
seizurePercentVector = 100*seizDurationVector/sessionTotalDuration;

%% Get correlation values
groupCorrDataSubtr = cell(size(groupCorrData));
groupCorrDataSubtrAvg = cell(size(groupCorrData));
for r1 = 1:size(groupCorrData,1)
    for c1 = 5:size(groupCorrData,2)
        groupCorrDataSubtr{r1,c1}(:,1) = groupCorrData{r1,c1}(:,2)-groupCorrData{r1,c1}(:,1);
        groupCorrDataSubtr{r1,c1}(:,2) = groupCorrIdx{r1,c1}(:,2);
        groupCorrDataSubtr{r1,c1}(:,3) = groupCorrIdx{r1,c1}(:,3);
        groupCorrDataSubtrAvg{r1,c1}(1,1) = nanmean(groupCorrDataSubtr{r1,c1}(1:3,1));
        groupCorrDataSubtrAvg{r1,c1}(2,1) = nanmean(groupCorrDataSubtr{r1,c1}(4:6,1));
        groupCorrDataSubtrAvg{r1,c1}(3,1) = nanmean(groupCorrDataSubtr{r1,c1}(7:8,1));
        groupCorrDataSubtrAvg{r1,c1}(4,1) = nanmean(groupCorrDataSubtr{r1,c1}(9:10,1));
        groupCorrDataSubtrAvg{r1,c1}(5,1) = nanmean(groupCorrDataSubtr{r1,c1}(11:12,1));
        groupCorrDataSubtrAvg{r1,c1}(6,1) = nanmean(groupCorrDataSubtr{r1,c1}(13:14,1));
        groupCorrDataSubtrAvg{r1,c1}(7,1) = nanmean(groupCorrDataSubtr{r1,c1}(15:16,1));
        groupCorrDataSubtrAvg{r1,c1}(8,1) = nanmean(groupCorrDataSubtr{r1,c1}(17:18,1));
        groupCorrDataSubtrAvg{r1,c1}(9,1) = nanmean(groupCorrDataSubtr{r1,c1}(19:20,1));
    end
end

%% Plot results change in correlation vs, seizures
plotTitles = {{'F'}; {'M'}; {'C'}; {'S'}; {'R'}; {'V'}};
% Font sizes
axisFont = 22;
textFont = 18;
axisLabelFont = 28;
dottedLineWidth = 2;
markSize = 12;

% Preallocate cells
p = cell(size(groupCorrData));
R = cell(size(groupCorrData));
% p-values of correlation
pVals = cell(size(groupCorrData));
f = cell(size(groupCorrData));
yLimits = [-2.2 2];
xLimits = [0.8*min(seizDurationVector) 1.05*max(seizDurationVector)];
% # points to interpolate fitted line
nPoints = 100;

figure; set (gcf,'color','w')
for r1=1:6,
    for c1=5:7,
        if c1 == 5 
            plotType = 'ro';
            lineType = 'r-';
        elseif c1 == 6
            plotType = 'bx';
            lineType = 'b-.';
        else
            plotType = 'k^';
            lineType = 'k--';
        end
        subplot(2,3,r1)
        x = seizDurationVector(~isnan(groupCorrDataSubtrAvg{r1,c1}));
        y = groupCorrDataSubtrAvg{r1,c1}(~isnan(groupCorrDataSubtrAvg{r1,c1}));
        plot(x, y, plotType, 'LineWidth', dottedLineWidth, 'MarkerSize', markSize)
        p{r1,c1} = polyfit(x,y,1);
        % Measure of correlation r^2
        [R{r1,c1} pVals{r1,c1}] = corr(x,y);
        R{r1,c1} = R{r1,c1} .^ 2;
        if c1 == 5 
            text(40, -0.8, sprintf('r^2(HbO_2)=%0.2f,p=%0.2f', R{r1,c1}, pVals{r1,c1}),...
                'FontSize', textFont, 'Color', 'r')
        elseif c1 == 6
            text(40, -1.25, sprintf('r^2(HbR)=%0.2f,p=%0.2f', R{r1,c1}, pVals{r1,c1}),...
                'FontSize', textFont, 'Color', 'b')
        else
            text(40, -1.70, sprintf('r^2(CBF)=%0.2f,p=%0.2f', R{r1,c1}, pVals{r1,c1}),...
                'FontSize', textFont, 'Color', 'k')
        end
        hold on
        f{r1,c1} = polyval(p{r1,c1},x);
        xinterp = linspace(x(1), x(end), nPoints);
        plot(x, f{r1,c1}, lineType, 'LineWidth', dottedLineWidth)
        title(plotTitles{r1}, 'FontSize', axisLabelFont);
        set(gca,'FontSize',axisFont);
        ylim(yLimits)
        xlim(xLimits)
        if r1 == 1
            ylabel('z_{4AP}(r) - z_0(r)','FontSize',axisLabelFont);
        end
        if r1 == 1 || r1 == 2 || r1 == 3
            set(gca, 'XTick', []);
        end
        if r1 == 2 || r1 == 3 || r1 == 5 || r1 == 6
            set(gca, 'YTick', []);
        end
        if r1 == 5
            xlabel('Seizure duration [s]','FontSize',axisLabelFont);
        end
    end
    if r1 == 6
        legend({'HbO_2';'HbO_2 fit';'HbR';'HbR fit';'Flow'; 'Flow fit'})
    end
end

%% Print graphics
% addpath(genpath('D:\Edgar\ssoct\Matlab'))
% export_fig(fullfile('D:\Edgar\Documents\Dropbox\Docs\Epilepsy\figs','z_vs_seizureReviewers'),'-png',gcf)

%% Save results
save (fullfile('F:\Edgar\Data\IOS_Results\corr_Results_seizuresReviewers','seizDur.mat'),...
    'seizureDuration', 'seizDurationVector','groupCorrData', 'groupCorrDataSubtr',...
    'groupCorrDataSubtrAvg', 'seizurePercentVector', 'R', 'p', 'pVals')


%% Plot results change in correlation vs, seizures duration % (FIG4 in JBO paper)
clear
% Load saved data
load('F:\Edgar\Data\IOS_Results\corr_Results_seizuresReviewers\seizDur.mat');
plotTitles = {{'F'}; {'M'}; {'C'}; {'S'}; {'R'}; {'V'}};
% Font sizes
axisFont = 22;
textFont = 20;
axisLabelFont = 30;
dottedLineWidth = 3;
markSize = 14;
% significance threshold
alpha = 0.05;

% Preallocate cells
f = cell(size(groupCorrData));
finterp = cell(size(groupCorrData));
yLimits = [-2.2 2];
xLimits = [0.8*min(seizurePercentVector) 1.05*max(seizurePercentVector)];
% # points to interpolate fitted line
nPoints = 100;

% New figure
h = figure; set (gcf,'color','w')
for r1=1:6,
    for c1=5:7,
        if c1 == 5 
            plotType = 'o';
            lineType = '-';
            % Color red
            LineColor = [161 0 64]/255;
        elseif c1 == 6
            plotType = 'x';
            lineType = '-.';
            % Color blue
            LineColor = [0 138 255]/255;
        else
            plotType = '^';
            lineType = '--';
            % Color black
            LineColor = 'k';
        end
        subplot(2,3,r1)
        x = seizurePercentVector(~isnan(groupCorrDataSubtrAvg{r1,c1}));
        y = groupCorrDataSubtrAvg{r1,c1}(~isnan(groupCorrDataSubtrAvg{r1,c1}));
        plot(x, y, plotType, 'LineWidth', dottedLineWidth, 'MarkerSize', markSize,...
            'Color', LineColor)
        p{r1,c1} = polyfit(x,y,1);
        if pVals{r1,c1} <= alpha
            % Plot only significant results labels
            if c1 == 5
                text(5, -0.8, sprintf('r^2(HbO_2)=%0.2f*', R{r1,c1}),...
                    'FontSize', textFont, 'FontWeight', 'b', 'Color', LineColor)
            elseif c1 == 6
                text(5, -1.25, sprintf('r^2(HbR)=%0.2f*', R{r1,c1}),...
                    'FontSize', textFont, 'FontWeight', 'b',  'Color', LineColor)
            else
                text(5, -1.70, sprintf('r^2(CBF)=%0.2f*', R{r1,c1}),...
                    'FontSize', textFont, 'FontWeight', 'b', 'Color', LineColor)
            end
        end
        hold on
        % Fitted line (1st degree polynomial)
        f{r1,c1} = polyval(p{r1,c1},x);
        plot(x, f{r1,c1}, lineType, 'LineWidth', dottedLineWidth, 'Color', LineColor)
%         title(plotTitles{r1}, 'FontSize', axisLabelFont);
        set(gca,'FontSize',axisFont);
        ylim(yLimits)
        xlim(xLimits)
        if r1 == 1
            ylabel('z_{4AP}(r) - z_0(r)','FontSize',axisLabelFont);
        end
        if r1 == 1 || r1 == 2 || r1 == 3
            set(gca, 'XTick', []);
        end
        if r1 == 2 || r1 == 3 || r1 == 5 || r1 == 6
            set(gca, 'YTick', []);
        end
        if r1 == 5
            xlabel('Seizure duration [%]','FontSize',axisLabelFont);
        end
    end
    if r1 == 6
%         legend({'HbO_2';'HbO_2 fit';'HbR';'HbR fit';'CBF'; 'CBF fit'})
    end
end

job.figSize = [20 10];
job.figRes = 1200;

% Specify window units
set(h, 'units', 'inches')
% Change figure and paper size
set(h, 'Position', [0.1 0.1 job.figSize(1) job.figSize(2)])
set(h, 'PaperPosition', [0.1 0.1 job.figSize(1) job.figSize(2)])

%% Print graphics

% addpath(genpath('D:\Edgar\ssoct\Matlab'))
% export_fig(fullfile('D:\Edgar\Documents\Dropbox\Docs\Epilepsy\figs','z_vs_seizurePercReviewers'),'-png',gcf)

clc
% Save as PNG at the user-defined resolution
print(h, '-dpng', ...
    fullfile('D:\Edgar\Documents\Dropbox\Docs\Epilepsy\figs', 'z_vs_seizurePercReviewers'),...
    sprintf('-r%d',job.figRes));
% Return the property to its default
set(h, 'units', 'pixels')
disp('Printing done!')


% EOF
