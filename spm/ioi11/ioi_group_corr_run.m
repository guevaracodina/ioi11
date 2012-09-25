function out = ioi_group_corr_run(job)
% Performs a group comparison based on bilateral correlation between seeds
% time-courses, though a paired t-test data on the correlation before the 4-AP
% injection and its value after the epileptogenic injection. Performs the t-test
% for each pair of seeds.
%_______________________________________________________________________________
% Copyright (C) 2012 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________

% For each subject:
% - Display resumeexp.txt in figures windows (to check which sessions are
%   control/treatment)
% - choose the pairs of ROIs (1-2, 3-4, ..., 11-12 by default) (batch)
% - choose the control sessions
% - choose the post-injection sessions
% for each color
%   - get the seed-to-seed correlation matrix
%   - transform Pearson's r to Fisher's z
%   - choose only 6 values for each mouse, each color
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% For all the subjects
% - Choose parent folder (batch)
% - Do a separate statistical for each seed data (paired t-test, wil)
% - Display a graph with ROI labels
% - Save a .csv file in parent folder


% Initialize cell to hold group data
groupCorrData = cell([size(job.paired_seeds, 1) 7]);
groupCorrIdx  = cell([size(job.paired_seeds, 1) 7]);

if isfield (job,'derivative')
    groupCorrDataDiff = cell([size(job.paired_seeds, 1) 7]);
end

%Big loop over subjects
for SubjIdx = 1:size(job.IOImat, 1)
    try
        tic
        %Load IOI.mat information
        [IOI IOImat dir_ioimat]= ioi_get_IOI(job,SubjIdx);
        
        if ~isfield(IOI.fcIOS.corr,'corrMatrixOK') % correlation matrix OK
            disp(['No seed-to-seed correlation matrix available for subject ' int2str(SubjIdx) ' ... skipping correlation map']);
        else
            if ~isfield(IOI.fcIOS.corr,'corrGroupOK') || job.force_redo
                % Get colors to include information
                IC = job.IC;
                colorNames = fieldnames(IOI.color);
                
                % File name where correlation data is saved
                IOI.fcIOS.corr(1).fnameGroup = fullfile(job.parent_results_dir{1},'group_corr_pair_seeds.mat');
                
                % File name where correlation data of 1st derivative is saved
                if isfield (job,'derivative')
                    IOI.fcIOS.corr(1).fnameGroupDiff = fullfile(job.parent_results_dir{1},'group_corr_pair_seeds_diff.mat');
                end
                
                % - Display resumeexp.txt in command window (to check which sessions are
                %   control/treatment)
                % fileread(fullfile(IOI.dir.dir_subj_raw,'resumeExp.txt'))
                
                if isfield(job.AutoSessionChoice, 'ManualSessions')
                    % - choose the control sessions
                    [ctrlSessList, sts] = cfg_getfile(Inf,'dir','Select control sessions',[], IOI.dir.dir_subj_res, 'S[0-9]+');
                    
                    % - choose the post-injection sessions
                    [treatSessList, sts] = cfg_getfile(Inf,'dir','Select 4-AP sessions',[], IOI.dir.dir_subj_res, 'S[0-9]+');
                    
                    ctrlSessVec = [];
                    % Get control session numbers
                    for iSess = 1:numel(ctrlSessList)
                        sessString = regexp(ctrlSessList{iSess},'\\S[0-9]+\\','match','once');
                        sessString = regexp(sessString,'[0-9]+','match','once');
                        ctrlSessVec(iSess) = str2double(sessString);
                    end
                    
                    treatSessVec = [];
                    % Get treatment session numbers
                    for iSess = 1:numel(treatSessList)
                        sessString = regexp(treatSessList{iSess},'\\S[0-9]+\\','match','once');
                        sessString = regexp(sessString,'[0-9]+','match','once');
                        treatSessVec(iSess) = str2double(sessString);
                    end
                elseif isfield(job.AutoSessionChoice, 'AutoSessions')
                    % Retrieve sessions
                    ctrlSessVec = job.AutoSessionChoice.AutoSessions.control_sessions{SubjIdx};
                    treatSessVec = job.AutoSessionChoice.AutoSessions.treatment_sessions{SubjIdx};
                else
                    % Do nothing for the moment
                end
                
                % Get indices for the paired sessions
                % 1st column contains the indices to the control sessions
                % 2nd column contains the indices to the 4-AP sessions
                idxSess = CombVec(ctrlSessVec, treatSessVec)';
                
                % Additional 3rd column with subject index
                idxSess(:,3) = SubjIdx;
                
                % Get the seed-to-seed correlation matrix
                load(IOI.fcIOS.corr.corrMatrixFname)
                
                if isfield (job,'derivative')
                    % Get the seed-to-seed correlation of 1st derivative matrix
                    load(IOI.fcIOS.corr.corrMatrixDiffFname)
                end
                
                % Initialize cell to hold paired seeds correlation data
                pairedSeeds = cell([size(job.paired_seeds, 1) 1]);
                for iCell = 1:size(job.paired_seeds, 1),
                    pairedSeeds{iCell} = cell([size(ctrlSessVec,1) + size(treatSessVec, 1) 7]);
                    if isfield (job,'derivative')
                        pairedSeedsDiff{iCell} = cell([size(ctrlSessVec,1) + size(treatSessVec, 1) 7]);
                    end
                end

                % Loop over available colors
                for c1 = 1:size(IOI.fcIOS.corr.corrMapName{1}, 2)
                    doColor = ioi_doColor(IOI,c1,IC);
                    if doColor
                        %skip laser - only extract for flow
                        if ~(IOI.color.eng(c1)==IOI.color.laser)
                            % Loop over sessions
                            for s1 = [ctrlSessVec treatSessVec],
                                % Get current correlation matrix
                                currCorrMat = seed2seedCorrMat{1}{s1,c1};
                                % transform Pearson's r to Fisher's z
                                currCorrMat = fisherz(currCorrMat);
                                if isnan(currCorrMat)
                                    for iROI = 1:size(job.paired_seeds, 1)
                                        pairedSeeds{iROI}{s1,c1} = NaN;
                                    end
                                else
                                    % choose only 6 values for each mouse, each color
                                    for iROI = 1:size(job.paired_seeds, 1)
                                        pairedSeeds{iROI}{s1,c1} = currCorrMat(job.paired_seeds(iROI,1), job.paired_seeds(iROI,2));
                                    end
                                end
                                if isfield (job,'derivative')
                                    % Get current correlation matrix
                                    currCorrMatDiff = seed2seedCorrMatDiff{1}{s1,c1};
                                    % transform Pearson's r to Fisher's z
                                    currCorrMatDiff = fisherz(currCorrMatDiff);
                                    if isnan(currCorrMatDiff)
                                        for iROI = 1:size(job.paired_seeds, 1)
                                            pairedSeedsDiff{iROI}{s1,c1} = NaN;
                                        end
                                    else
                                        % choose only 6 values for each mouse, each color
                                        for iROI = 1:size(job.paired_seeds, 1)
                                            pairedSeedsDiff{iROI}{s1,c1} = currCorrMatDiff(job.paired_seeds(iROI,1), job.paired_seeds(iROI,2));
                                        end
                                    end
                                    
                                end % derivative
                            end % sessions loop
                        end
                        % Arrange the paired seeds according to idxSessions
                        fprintf('Retrieving data %s C%d (%s)...\n',IOI.subj_name, c1,colorNames{1+c1});
                        tmpArray = zeros([size(job.paired_seeds, 1), size(idxSess)]);
                        if isfield (job,'derivative')
                            tmpArrayDiff = zeros([size(job.paired_seeds, 1), size(idxSess)]);
                        end
                        for iROI = 1:size(job.paired_seeds, 1)
                            for iSess = 1:size(idxSess, 1)
                                tmpArray(iROI,iSess,1) = pairedSeeds{iROI}{idxSess(iSess,1), c1};
                                tmpArray(iROI,iSess,2) = pairedSeeds{iROI}{idxSess(iSess,2), c1};
                                if isfield (job,'derivative')
                                    tmpArrayDiff(iROI,iSess,1) = pairedSeedsDiff{iROI}{idxSess(iSess,1), c1};
                                    tmpArrayDiff(iROI,iSess,2) = pairedSeedsDiff{iROI}{idxSess(iSess,2), c1};
                                end
                            end
                        end
                        for iROI = 1:size(job.paired_seeds, 1)
                            groupCorrData{iROI,c1} = [groupCorrData{iROI,c1}; squeeze(tmpArray(iROI,:,:))];
                            groupCorrIdx{iROI,c1} = [groupCorrIdx{iROI,c1}; idxSess];
                            if isfield (job,'derivative')
                                groupCorrDataDiff{iROI,c1} = [groupCorrDataDiff{iROI,c1}; squeeze(tmpArrayDiff(iROI,:,:))];
                            end
                        end
                    % else
                        % correlation map failed!
                        % IOI.fcIOS.corr(1).corrGroupOK{1}{1, c1} = false;
                        % fprintf('Paired seeds correlation analysis failed! C%d (%s)...\n',c1,colorNames{1+c1});
                    end
                end % colors loop
                % Save group correlation data data (data is appended for more subjects)
                save(IOI.fcIOS.corr(1).fnameGroup,'groupCorrData','groupCorrIdx','groupCorrDataDiff')
                if isfield (job,'derivative')
                    save(IOI.fcIOS.corr(1).fnameGroupDiff,'groupCorrDataDiff');
                end
                % Save IOI matrix
                save(IOImat,'IOI');
            end % correlation OK or redo job
        end % corrMap OK
        disp(['Elapsed time: ' datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS')]);
        disp(['Subject ' int2str(SubjIdx) ' (' IOI.subj_name ')' ' complete']);
        out.IOImat{SubjIdx} = IOImat;
    catch exception
        out.IOImat{SubjIdx} = IOImat;
        disp(exception.identifier)
        disp(exception.stack(1))
    end
end % Big loop over subjects

% For all the subjects
% Choose parent folder (batch)
% Loop over available colors and paired seeds
% Do group statistical tests
% Plot/print graphics
for c1 = 1:size(IOI.fcIOS.corr.corrMapName{1}, 2)
    doColor = ioi_doColor(IOI,c1,IC);
    if doColor
        colorOK = true;
        %skip laser - only extract for flow
        if ~(IOI.color.eng(c1)==IOI.color.laser)
            
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
                    treatmentlGroup = groupCorrData{iSeeds,c1}(:,2);
                    % ignore NaN values
                    treatmentlGroup = treatmentlGroup(~isnan(treatmentlGroup));
                    % Perform such test
                    [statTest(1).w(1).P{iSeeds,c1}, statTest(1).w(1).H{iSeeds,c1},...
                        statTest(1).w(1).STATS{iSeeds,c1}] = ranksum...
                        (ctrlGroup, treatmentlGroup, 'alpha', job.alpha);
                    statTest(1).w(1).id = 'Wilcoxon rank sum test';
                 end % Wilcoxon test
                 
            end % paired-seeds loop
            
            % Show standard error bars instead of standard deviation
            if job.stderror
                % std error bars: sigma/sqrt(N)
                e = e / sqrt(size(groupCorrData{iSeeds,c1}, 1));
            end
            
            if isfield (job,'derivative')
                for iSeeds = 1:size(job.paired_seeds, 1)
                    % Average of control group
                    meanCorrDiff{iSeeds,c1}(1) = nanmean(groupCorrDataDiff{iSeeds,c1}(:,1));
                    % Average of treatment group
                    meanCorrDiff{iSeeds,c1}(2) = nanmean(groupCorrDataDiff{iSeeds,c1}(:,2));
                    % Standard deviation of control group
                    stdCorrDiff{iSeeds,c1}(1) = nanstd(groupCorrDataDiff{iSeeds,c1}(:,1));
                    % Standard deviation oftreatment group
                    stdCorrDiff{iSeeds,c1}(2) = nanstd(groupCorrDataDiff{iSeeds,c1}(:,2));
                    
                    % Used to plot bar graphs with errorbars
                    yDiff(iSeeds,:) = meanCorrDiff{iSeeds,c1};
                    eDiff(iSeeds,:) = stdCorrDiff{iSeeds,c1};
                    
                    % Paired-sample t-test
                    if job.ttest1
                        [statTestDiff(1).t(1).H{iSeeds,c1}, statTestDiff(1).t(1).P{iSeeds,c1}, ...
                            statTestDiff(1).t(1).CI{iSeeds,c1}, statTestDiff(1).t(1).STATS{iSeeds,c1}] ...
                            = ttest(...
                            groupCorrDataDiff{iSeeds,c1}(:,1), groupCorrDataDiff{iSeeds,c1}(:,2),...
                            job.alpha,'both');
                        statTestDiff(1).t(1).id = 'Paired-sample t-test(1st derivative)';
                    end % t-test
                    
                    % Wilcoxon rank sum test
                    if job.wilcoxon1
                        ctrlGroupDiff = groupCorrDataDiff{iSeeds,c1}(:,1);
                        % ignore NaN values
                        ctrlGroupDiff = ctrlGroupDiff(~isnan(ctrlGroupDiff));
                        treatmentlGroupDiff = groupCorrDataDiff{iSeeds,c1}(:,2);
                        % ignore NaN values
                        treatmentlGroupDiff = treatmentlGroupDiff(~isnan(treatmentlGroupDiff));
                        % Perform such test
                        [statTestDiff(1).w(1).P{iSeeds,c1}, statTestDiff(1).w(1).H{iSeeds,c1},...
                            statTestDiff(1).w(1).STATS{iSeeds,c1}] = ranksum...
                            (ctrlGroupDiff, treatmentlGroupDiff, 'alpha', job.alpha);
                        statTestDiff(1).w(1).id = 'Wilcoxon rank sum test(1st derivative)';
                    end % Wilcoxon test
                    
                end % paired-seeds loop
                
                % Show standard error bars instead of standard deviation
                if job.stderror
                    % std error bars: sigma/sqrt(N)
                    eDiff = eDiff / sqrt(size(groupCorrDataDiff{iSeeds,c1}, 1));
                end
                % Save results in .mat file
                save(IOI.fcIOS.corr(1).fnameGroupDiff,'groupCorrDataDiff','groupCorrIdx',...
                'meanCorrDiff','stdCorrDiff','statTestDiff');
            end
            
            % Plot results
            private_plot_group_corr_test(job, IOI, c1, e, y, statTest);
            
            % Plot results based on 1st derivative
            private_plot_group_corr_test_diff(job, IOI, c1, eDiff, yDiff, statTestDiff);
            
            % Save results in .mat file
            save(IOI.fcIOS.corr(1).fnameGroup,'groupCorrData','groupCorrIdx',...
                'meanCorr','stdCorr','statTest');
        end
    end
end % loop over colors
% Group comparison of bilateral correlation succesful!
fprintf('Group comparison of bilateral correlation succesful!\n');

function private_plot_group_corr_test(job, IOI, c1, e, y, statTest)
% Plots statistical analysis group results
colorNames = fieldnames(IOI.color);
% Positioning factor for the * mark, depends on max data value at the given seed
starPosFactor = 1.05;

if job.ttest1
    % Display a graph with ROI labels
    if job.generate_figures
        % Display plots on SPM graphics window
        h = spm_figure('GetWin', 'Graphics');
        spm_figure('Clear', 'Graphics');
        % Custom bar graphs with error bars (1st arg: error)
        barwitherr(e, y)
        % Display colormap according to the contrast
        switch(c1)
            case 5
                % HbO contrast
                colormap([1 0 0; 1 1 1]);
            case 6
                % HbR contrast
                colormap([0 0 1; 1 1 1]);
            case 7
                % Flow contrast
                colormap([0.5 0.5 0.5; 1 1 1])
            otherwise
                colormap(gray)
        end
        title(sprintf('Bilateral Correlation before/after 4-AP C%d(%s). T-test (*p<%.2g)',...
            c1,colorNames{1+c1},job.alpha),'interpreter','none','FontSize',12)
        ylabel('Functional correlation z(r)','FontSize',18)
        set(gca,'XTickLabel',{'F', 'M', 'C', 'S', 'R', 'V'},'FontWeight', 'b')
        set(gca,'FontSize',14)
        legend({'Control' '4-AP'},'FontSize',12)
        % Show a * when a significant difference is found.
        for iSeeds = 1:size(job.paired_seeds, 1)
            if statTest(1).t(1).H{iSeeds,c1}
                if max(y(iSeeds,:))>=0
                    yPos = starPosFactor*(max(y(iSeeds,:)) + max(e(iSeeds,:)));
                else
                    yPos = starPosFactor*(min(y(iSeeds,:)) - max(e(iSeeds,:)));
                end
                xPos = iSeeds;
                text(xPos, yPos, '*', 'FontSize', 24, 'FontWeight', 'b');
            end
        end
        if job.save_figures
            newName = sprintf('groupCorr_Ttest_C%d_(%s)',c1,colorNames{1+c1});
            % Save as EPS
            spm_figure('Print', 'Graphics', fullfile(job.parent_results_dir{1}, newName));
            % Save as PNG
            print(h, '-dpng', fullfile(job.parent_results_dir{1},newName), '-r300');
        end
    end % end generate figures
end

if job.wilcoxon1
    % Display a graph with ROI labels
    if job.generate_figures
        % Display plots on SPM graphics window
        h = spm_figure('GetWin', 'Graphics');
        spm_figure('Clear', 'Graphics');
        % Custom bar graphs with error bars (1st arg: error)
        barwitherr(e, y)
        % Display colormap according to the contrast
        switch(c1)
            case 5
                % HbO contrast
                colormap([1 0 0; 1 1 1]);
            case 6
                % HbR contrast
                colormap([0 0 1; 1 1 1]);
            case 7
                % Flow contrast
                colormap([0.5 0.5 0.5; 1 1 1])
            otherwise
                colormap(gray)
        end
        title(sprintf('Bilateral Correlation before/after 4-AP C%d(%s). Wilcoxon *(p<%.2g)',...
            c1,colorNames{1+c1},job.alpha),'interpreter','none','FontSize',12)
        ylabel('Functional correlation z(r)','FontSize',18)
        set(gca,'XTickLabel',{'F', 'M', 'C', 'S', 'R', 'V'},'FontWeight', 'b')
        set(gca,'FontSize',14)
        legend({'Control' '4-AP'},'FontSize',12)
        % Show a * when a significant difference is found.
        for iSeeds = 1:size(job.paired_seeds, 1)
            if statTest(1).w(1).H{iSeeds,c1}
                if max(y(iSeeds,:))>=0
                    yPos = starPosFactor*(max(y(iSeeds,:)) + max(e(iSeeds,:)));
                else
                    yPos = starPosFactor*(min(y(iSeeds,:)) - max(e(iSeeds,:)));
                end
                xPos = iSeeds;
                text(xPos,yPos,'*', 'FontSize', 24, 'FontWeight', 'b');
            end
        end
        if job.save_figures
            newName = sprintf('groupCorr_Wtest_C%d_(%s)',c1,colorNames{1+c1});
            % Save as EPS
            spm_figure('Print', 'Graphics', fullfile(job.parent_results_dir{1},newName));
            % Save as PNG
            print(h, '-dpng', fullfile(job.parent_results_dir{1},newName), '-r300');
        end
    end % End generate figures
end

function private_plot_group_corr_test_diff(job, IOI, c1, eDiff, yDiff, statTestDiff)
if isfield (job,'derivative')
    % Plots statistical analysis group results
    colorNames = fieldnames(IOI.color);
    % Positioning factor for the * mark, depends on max data value at the given seed
    starPosFactor = 1.05;
    
    if job.ttest1
        % Display a graph with ROI labels
        if job.generate_figures
            % Display plots on SPM graphics window
            h = spm_figure('GetWin', 'Graphics');
            spm_figure('Clear', 'Graphics');
            % Custom bar graphs with error bars (1st arg: error)
            barwitherr(eDiff, yDiff)
            % Display colormap according to the contrast
            switch(c1)
                case 5
                    % HbO contrast
                    colormap([1 0 0; 1 1 1]);
                case 6
                    % HbR contrast
                    colormap([0 0 1; 1 1 1]);
                case 7
                    % Flow contrast
                    colormap([0.5 0.5 0.5; 1 1 1])
                otherwise
                    colormap(gray)
            end
            title(sprintf('Bilateral Correlation before/after 4-AP C%d(%s) Diff. T-test (*p<%.2g)',...
                c1,colorNames{1+c1},job.alpha),'interpreter','none','FontSize',12)
            ylabel('Functional correlation z(r)','FontSize',18)
            set(gca,'XTickLabel',{'F', 'M', 'C', 'S', 'R', 'V'},'FontWeight', 'b')
            set(gca,'FontSize',14)
            legend({'Control' '4-AP'},'FontSize',12)
            % Show a * when a significant difference is found.
            for iSeeds = 1:size(job.paired_seeds, 1)
                if statTestDiff(1).t(1).H{iSeeds,c1}
                    if max(yDiff(iSeeds,:))>=0
                        yPos = starPosFactor*(max(yDiff(iSeeds,:)) + max(eDiff(iSeeds,:)));
                    else
                        yPos = starPosFactor*(min(yDiff(iSeeds,:)) - max(eDiff(iSeeds,:)));
                    end
                    xPos = iSeeds;
                    text(xPos, yPos, '*', 'FontSize', 24, 'FontWeight', 'b');
                end
            end
            if job.save_figures
                newName = sprintf('groupCorr_Ttest_C%d_(%s)_diff',c1,colorNames{1+c1});
                % Save as EPS
                spm_figure('Print', 'Graphics', fullfile(job.parent_results_dir{1}, newName));
                % Save as PNG
                print(h, '-dpng', fullfile(job.parent_results_dir{1},newName), '-r300');
            end
        end % end generate figures
    end
    
    if job.wilcoxon1
        % Display a graph with ROI labels
        if job.generate_figures
            % Display plots on SPM graphics window
            h = spm_figure('GetWin', 'Graphics');
            spm_figure('Clear', 'Graphics');
            % Custom bar graphs with error bars (1st arg: error)
            barwitherr(eDiff, yDiff)
            % Display colormap according to the contrast
            switch(c1)
                case 5
                    % HbO contrast
                    colormap([1 0 0; 1 1 1]);
                case 6
                    % HbR contrast
                    colormap([0 0 1; 1 1 1]);
                case 7
                    % Flow contrast
                    colormap([0.5 0.5 0.5; 1 1 1])
                otherwise
                    colormap(gray)
            end
            title(sprintf('Bilateral Correlation before/after 4-AP C%d(%s) Diff. Wilcoxon (*p<%.2g)',...
                c1,colorNames{1+c1},job.alpha),'interpreter','none','FontSize',12)
            ylabel('Functional correlation z(r)','FontSize',18)
            set(gca,'XTickLabel',{'F', 'M', 'C', 'S', 'R', 'V'},'FontWeight', 'b')
            set(gca,'FontSize',14)
            legend({'Control' '4-AP'},'FontSize',12)
            % Show a * when a significant difference is found.
            for iSeeds = 1:size(job.paired_seeds, 1)
                if statTestDiff(1).w(1).H{iSeeds,c1}
                    if max(yDiff(iSeeds,:))>=0
                        yPos = starPosFactor*(max(yDiff(iSeeds,:)) + max(eDiff(iSeeds,:)));
                    else
                        yPos = starPosFactor*(min(yDiff(iSeeds,:)) - max(eDiff(iSeeds,:)));
                    end
                    xPos = iSeeds;
                    text(xPos,yPos,'*', 'FontSize', 24, 'FontWeight', 'b');
                end
            end
            if job.save_figures
                newName = sprintf('groupCorr_Wtest_C%d_(%s)_diff',c1,colorNames{1+c1});
                % Save as EPS
                spm_figure('Print', 'Graphics', fullfile(job.parent_results_dir{1},newName));
                % Save as PNG
                print(h, '-dpng', fullfile(job.parent_results_dir{1},newName), '-r300');
            end
        end % End generate figures
    end % Wilcoxon
end % derivative

% EOF
