%% script_peak_threshold
load('E:\Edgar\Data\IOS_Results\corr_Results_seizures\group_corr_pair_seeds.mat')
load('E:\Edgar\Data\IOS_Results\corr_Results_seizures\group_corr_pair_seeds_diff.mat')

%% Preallocate and gather data for the sessions after the 4AP injection
seizureDuration = num2cell(nan(max(groupCorrIdx{1,5}(:,3)), 5));
% seizureDuration{iSubject, iSession}
seizureDuration{1, 4} = mean([203-158; 375-332; 660-627; 838-806]);
seizureDuration{1, 5} = mean([244-209; 433-401; 729-702]);

% seizureDuration{2, 3} = NaN;
% seizureDuration{2, 4} = NaN;

% Only 1 4AP session
seizureDuration{3-1, 3} = mean([81-54; 332-299; 852-818]);

seizureDuration{4-1, 3} = mean([334-220]);
seizureDuration{4-1, 4} = mean([392-321; 812-700]);

% seizureDuration{5, 3} = NaN;
% seizureDuration{5, 4} = NaN;
% 
% seizureDuration{6, 3} = NaN;
% seizureDuration{6, 4} = NaN;

seizureDuration{7-3, 3} = mean([200-162; 384-354; 748-718]);
seizureDuration{7-3, 4} = mean([382-349; 566-534]);

seizureDuration{9-4, 3} = mean([459-329]);
seizureDuration{9-4, 4} = mean([432-271]);

% seizureDuration{10-1, 3} = NaN;
% seizureDuration{10-1, 4} = NaN;


%% Only subjects where seizures where observed
seizureDurationTrans = seizureDuration';
seizDurationVector = cell2mat(seizureDurationTrans(~isnan(cell2mat(seizureDurationTrans))));

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

% %% Arrange data
% clear seizureDurationComb
% for idxSub = 1:size(seizureDuration,1)
%     for s1 = 3:size(seizureDuration,2)
%         for r1 = 1:size(groupCorrData,1)
%             for c1 = 5:size(groupCorrData,2)
%                 % s2 goes from 1 to 36 , it is the umber of combined sessions
%                 for s2 = 1:size(groupCorrIdx{r1,c1},1)
%                     if (groupCorrIdx{r1,c1}(s2,2) == s1) && (groupCorrIdx{r1,c1}(s2,3) == idxSub)
%                         % 1st column: seizure durations
%                         seizureDurationComb{s2, 1}{r1,c1} = seizureDuration{idxSub, s1};
%                         % 2nd column: 4-APsessions
%                         seizureDurationComb{s2, 2}{r1,c1} = groupCorrIdx{r1,c1}(s2, 2);
%                         % 3rd column: subject indices
%                         seizureDurationComb{s2, 3}{r1,c1} = groupCorrIdx{r1,c1}(s2, 3);
%                         % 4th column: Absolute loss of correlation
%                         seizureDurationComb{s2, 4}{r1,c1} = (groupCorrData{r1,c1}(s2,2) - groupCorrData{r1,c1}(s2,1));
%                         % 5th column: combined vectors: 1st row duration, 2nd
%                         % row loss of correlation
%                         seizureDurationComb{s2, 5}{r1,c1} = CombVec (seizureDurationComb{s2, 1}{r1,c1}', seizureDurationComb{s2, 4}{r1,c1})';
%                         % 6th column, 
%                         seizureDurationComb{s2, 6}{r1,c1} = [mean(seizureDurationComb{s2, 1}{r1,c1}) mean(seizureDurationComb{s2, 4}{r1,c1})];
%                     end
%                 end
%             end
%         end
%     end
% end
% 
% % Arrange final results
% lossVSseiz = cell(size(groupCorrData));
% for s2 = 1:size(groupCorrIdx{r1,c1},1)
%     for r1 = 1:size(groupCorrData,1)
%         for c1 = 5:size(groupCorrData,2)
%             lossVSseiz{r1, c1} = [lossVSseiz{r1, c1};  seizureDurationComb{s2, 6}{r1,c1}];
%         end
%     end
% end

%% Plot results
% Load saved data
load('E:\Edgar\Data\IOS_Results\corr_Results_seizures\seizDur.mat');
plotTitles = {{'F'}; {'M'}; {'C'}; {'S'}; {'R'}; {'V'}};
% Font sizes
axisFont = 22;
textFont = 18;
axisLabelFont = 26;
dottedLineWidth = 4;
markSize = 12;

%
p = cell(size(groupCorrData));
R = cell(size(groupCorrData));
f = cell(size(groupCorrData));
yLimits = [-2 2];
xLimits = [0.8*min(seizDurationVector) 1.05*max(seizDurationVector)];

figure; set (gcf,'color','w')
for r1=1:6,
    for c1=5:7,
        if c1 == 5 
            plotType = 'ro';
            lineType = 'r-';
        elseif c1 == 6
            plotType = 'bx';
            lineType = 'b--';
        else
            plotType = 'ks';
            lineType = 'k:';
        end
        subplot(2,3,r1)
        x = seizDurationVector(~isnan(groupCorrDataSubtrAvg{r1,c1}));
        y = groupCorrDataSubtrAvg{r1,c1}(~isnan(groupCorrDataSubtrAvg{r1,c1}));
        plot(x, y, plotType, 'LineWidth', dottedLineWidth, 'MarkerSize', markSize)
        p{r1,c1} = polyfit(x,y,1);
        % Measure of correlation r^2
        R{r1,c1} = corr(x,y).^2;
        if c1 == 5 
            text(40, -0.8, ['r^2(HbO)=', sprintf('%0.2f',R{r1,c1})],...
                'FontSize', textFont, 'Color', 'r')
        elseif c1 == 6
            text(40, -1.25, ['r^2(HbR)=', sprintf('%0.2f',R{r1,c1})],...
                'FontSize', textFont, 'Color', 'b')
        else
            text(40, -1.70, ['r^2(F)=', sprintf('%0.2f',R{r1,c1})],...
                'FontSize', textFont, 'Color', 'k')
        end
        hold on
        f{r1,c1} = polyval(p{r1,c1},x);
        plot(x, f{r1,c1}, lineType, 'LineWidth', dottedLineWidth)
        title(plotTitles{r1}, 'FontSize', axisLabelFont);
        set(gca,'FontSize',axisFont);
        ylim(yLimits)
        xlim(xLimits)
        if r1 == 1 || r1 == 4
            ylabel('z_{4AP}(r) - z_0(r)','FontSize',axisLabelFont+6);
        end
        if r1 == 5
            xlabel('Seizure duration [s]','FontSize',axisLabelFont);
        end
    end
    if r1 == 6
        legend({'HbO_2';'HbO_2 fit';'HbR';'HbR fit';'Flow'; 'Flow fit'})
    end
end

%% Save results
save (fullfile('E:\Edgar\Data\IOS_Results\corr_Results_seizures','seizDur.mat'),...
    'seizureDuration', 'seizDurationVector','groupCorrData', 'groupCorrDataSubtr',...
    'groupCorrDataSubtrAvg', 'R', 'p', 'f')

%% Print graphics
addpath(genpath('D:\Edgar\ssoct\Matlab'))
export_fig(fullfile('D:\Edgar\Documents\Dropbox\Docs\Epilepsy\figs','z_vs_seizure'),'-png',gcf)

%% Plot resting state somatosensory seeds time traces.
% Font sizes
axisFont = 18;
axisLabelFont = 20;
dottedLineWidth = 2;
markSize = 12;
xLimits = [200 500];

load('E:\Edgar\Data\IOS_Results\12_08_21,EG01\FiltNDownButter4ff\GLMfcIOS\ROIregress.mat')
r1 = 3;
s1 = 5;


figure; set(gcf,'color','w')
hold on

% HbO
c1 = 5;
plot(ROIregress{r1}{s1,c1}, 'r-', 'LineWidth', dottedLineWidth)
plot(ROIregress{r1+1}{s1,c1}, 'r--', 'LineWidth', dottedLineWidth)
corr(ROIregress{r1}{s1,c1}',ROIregress{r1+1}{s1,c1}')

% HbR
c1 = 6;
plot(ROIregress{r1}{s1,c1}, 'b-', 'LineWidth', dottedLineWidth)
plot(ROIregress{r1+1}{s1,c1}, 'b--', 'LineWidth', dottedLineWidth)
corr(ROIregress{r1}{s1,c1}',ROIregress{r1+1}{s1,c1}')

xlim(xLimits)
legend({'HbO_2 Left';'HbO_2 Right';'HbR Left';'HbR Right'})
set(gca,'FontSize',axisFont);
xlabel('time [s]','FontSize',axisLabelFont);
ylabel('\DeltaHb [mM]','FontSize',axisLabelFont);

%% Print graphics
addpath(genpath('D:\Edgar\ssoct\Matlab'))
export_fig(fullfile('D:\Edgar\Documents\Dropbox\Docs\Epilepsy\figs','resting_state_somato'),'-png',gcf)

