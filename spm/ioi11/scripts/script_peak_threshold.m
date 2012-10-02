%% script_peak_threshold
load('E:\Edgar\Data\IOS_Results\corr_Results_01-07\group_corr_pair_seeds.mat')


%% Preallocate and gather data for the sessions after the 4AP injection
seizureDuration = num2cell(nan(max(groupCorrIdx{1,5}(:,3)), 5));
% seizureDuration{iSubject}{iSession}
seizureDuration{1}{4} = [244-209; 433-401; 729-702];
seizureDuration{1}{5} = [244-209; 433-401; 729-702];

seizureDuration{2}{3} = NaN;
seizureDuration{2}{4} = NaN;

% Only 1 4AP session
seizureDuration{3}{3} = [81-54; 332-299; 852-818];

seizureDuration{4}{3} = [392-321; 812-700];
seizureDuration{4}{4} = [392-321; 812-700];

seizureDuration{5}{3} = NaN;
seizureDuration{5}{4} = NaN;

seizureDuration{6}{3} = NaN;
seizureDuration{6}{4} = NaN;

seizureDuration{7}{3} = [200-162; 384-354; 748-718];
seizureDuration{7}{4} = [200-162; 384-354; 748-718];


%% Plot results

%% Save results
save (fullfile('E:\Edgar\Data\IOS_Results\corr_seizure_01-07','seizDur.mat'),...
    'seizureDuration')
