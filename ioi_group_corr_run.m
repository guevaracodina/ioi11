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
% - Do a separate paired t-test for each seed data
% - Display a graph with ROI labels
% - Save a .csv file in parent folder


% Initialize cell to hold group data
groupCorrData = cell([length(job.paired_seeds) 7]);
groupCorrIdx  = cell([length(job.paired_seeds) 7]);


%Big loop over subjects
for SubjIdx=1:length(job.IOImat)
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
                
                % - Display resumeexp.txt in command window (to check which sessions are
                %   control/treatment)
                % fileread(fullfile(IOI.dir.dir_subj_raw,'resumeExp.txt'))
                
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
                
                % Get indices for the paired sessions
                idxSess = CombVec(ctrlSessVec, treatSessVec)';
                
                %   - get the seed-to-seed correlation matrix
                load(IOI.fcIOS.corr.corrMatrixFname)
                
                % Initialize cell to hold paired seeds correlation data
                pairedSeeds = cell([length(job.paired_seeds) 1]);
                for iCell = 1:length(job.paired_seeds),
                    pairedSeeds{iCell} = cell([length(ctrlSessVec)+length(treatSessVec) 7]);
                end

                % Loop over available colors
                for c1=1:length(IOI.fcIOS.corr.corrMapName{1})
                    doColor = ioi_doColor(IOI,c1,IC);
                    if doColor
                        colorOK = true;
                        %skip laser - only extract for flow
                        if ~(IOI.color.eng(c1)==IOI.color.laser)
                            % Loop over sessions
                            for s1 = [ctrlSessVec treatSessVec],
                                %% Do my stuff here!
                                % Get current correlation matrix
                                currCorrMat = seed2seedCorrMat{1}{s1,c1};
                                %   - transform Pearson's r to Fisher's z
                                currCorrMat = fisherz(currCorrMat);
                                if isnan(currCorrMat)
                                    for iROI = 1:length(job.paired_seeds)
                                        pairedSeeds{iROI}{s1,c1} = NaN;
                                    end
                                else
                                    %   - choose only 6 values for each mouse, each color
                                    for iROI = 1:length(job.paired_seeds)
                                        pairedSeeds{iROI}{s1,c1} = currCorrMat(job.paired_seeds(iROI,1), job.paired_seeds(iROI,2));
                                    end
                                end
                                % Group comparison of bilateral correlation succesful!
                                % fprintf('Paired seeds correlation analysis succesful! C%d (%s)...\n', c1,colorNames{1+c1});
                            end % sessions loop
                        end
                        % Arrange the paired seeds according to idxSessions
                        fprintf('Arranging data C%d (%s)...\n',c1,colorNames{1+c1});
                        tmpArray = zeros([length(job.paired_seeds), size(idxSess)]);
                        for iROI = 1:length(job.paired_seeds)
                            for iSess = 1:length(idxSess)
                                tmpArray(iROI,iSess,1) = pairedSeeds{iROI}{idxSess(iSess,1), c1};
                                tmpArray(iROI,iSess,2) = pairedSeeds{iROI}{idxSess(iSess,2), c1};
                            end
                        end
                        for iROI = 1:length(job.paired_seeds)
                            groupCorrData{iROI,c1} = [groupCorrData{iROI,c1}; squeeze(tmpArray(iROI,:,:))];
                            groupCorrIdx{iROI,c1} = [groupCorrIdx{iROI,c1}; idxSess];
                        end
                    % else
                        % correlation map failed!
                        % IOI.fcIOS.corr(1).corrGroupOK{1}{1, c1} = false;
                        % fprintf('Paired seeds correlation analysis failed! C%d (%s)...\n',c1,colorNames{1+c1});
                    end
                end % colors loop
                % Save group correlation data data
                save(IOI.fcIOS.corr(1).fnameGroup,'groupCorrData','groupCorrIdx')
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
% - Choose parent folder (batch)
% Loop over available colors
for c1=1:length(IOI.fcIOS.corr.corrMapName{1})
    doColor = ioi_doColor(IOI,c1,IC);
    if doColor
        colorOK = true;
        %skip laser - only extract for flow
        if ~(IOI.color.eng(c1)==IOI.color.laser)
            % - Do a separate paired t-test for each seed data
            for iSeeds = 1:length(job.paired_seeds)
                meanCorr{iSeeds,c1}(1) = nanmean(groupCorrData{iSeeds,c1}(:,1));
                meanCorr{iSeeds,c1}(2) = nanmean(groupCorrData{iSeeds,c1}(:,2));
                stdCorr{iSeeds,c1}(1) = nanstd(groupCorrData{iSeeds,c1}(:,1));
                stdCorr{iSeeds,c1}(2) = nanstd(groupCorrData{iSeeds,c1}(:,2));
                y(iSeeds,:) = meanCorr{iSeeds,c1};
                e(iSeeds,:) = stdCorr{iSeeds,c1};
                 
                [H{iSeeds,c1},P{iSeeds,c1},CI{iSeeds,c1},STATS{iSeeds,c1}] = ttest(...
                    groupCorrData{iSeeds,c1}(:,1), groupCorrData{iSeeds,c1}(:,2),...
                    job.alpha,'both');
            end
            % std error bars: sigma/sqrt(N)
            e=e/sqrt(length(groupCorrData{iSeeds,c1}));
            % - Display a graph with ROI labels
            if job.generate_figures
                % Display plots on SPM graphics window
                h = spm_figure('GetWin', 'Graphics');
                set(gcf,'color','w')
                barwitherr(e,y)
                % Display colormap according to the contrast
                switch(c1)
                    case 5
                        colormap([1 0 0; 1 1 1]);
                    case 6
                        colormap([0 0 1; 1 1 1]);
                    case 7 
                        colormap(gray)
                    otherwise
                        colormap(gray)
                end
                title(sprintf('%s Bilateral correlation before/after 4-AP C%d(%s)',IOI.subj_name,c1,colorNames{1+c1}),'interpreter','none','FontSize',12)
                ylabel('functional correlation z(r)','FontSize',12)
                set(gca,'XTickLabel',{'F', 'M', 'C', 'S', 'R', 'V'})
                set(gca,'FontSize',12)
                legend({'Control' '4-AP'},'FontSize',12)
                if job.save_figures
                    newName = sprintf('%s_groupCorr_C%d_(%s)',IOI.subj_name,c1,colorNames{1+c1});
                    % Save as EPS
                    spm_figure('Print', 'Graphics', fullfile(job.parent_results_dir{1},newName));
                    % Save as PNG
                    print(h, '-dpng', fullfile(job.parent_results_dir{1},newName), '-r300');
                end
            end
            % - Save a .csv file in parent folder
            save(IOI.fcIOS.corr(1).fnameGroup,'groupCorrData','groupCorrIdx',...
                'meanCorr','H','P','CI','STATS');
        end
    end
end % loop over colors
% correlation group analysis succesful!
% IOI.fcIOS.corr(1).corrGroupOK = true;
% Group comparison of bilateral correlation succesful!
fprintf('Group comparison of bilateral correlation succesful!\n');

% EOF
