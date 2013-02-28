function out = ioi_network_analyses_run(job)
% network analyses
%_______________________________________________________________________________
% Copyright (C) 2012 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________

% ------------------------------------------------------------------------------
% REMOVE AFTER FINISHING THE FUNCTION //EGC
% ------------------------------------------------------------------------------
% fprintf('Work in progress...\nEGC\n')
% out.IOImat = job.IOImat;
% return
% ------------------------------------------------------------------------------

tic
% CONN toolbox v.13
addpath(genpath('D:\Edgar\conn'))
% current folder
currFolder = pwd;
% Load first IOI matrix
[IOI IOImat dir_ioimat] = ioi_get_IOI(job, 1);
% Get sessions info
% [all_sessions selected_sessions] = ioi_get_sessions(job);
% Get colors to include information
IC = job.IC;
colorNames = fieldnames(IOI.color);
% Flag to indicate replacement of seeds names
REPLACE_ROI_NAME = true;
% Seeds names
% newROIname = {  'Frontal Right'
%                 'Frontal Left'
%                 'Motor Right'
%                 'Motor Left'
%                 'Cingulate Right'
%                 'Cingulate Left'
%                 'Somatosensory Right'
%                 'Somatosensory Left'
%                 'Retrosplenial Right'
%                 'Retrosplenial Left'
%                 'Visual Right'
%                 'Visual Left' };
newROIname = {  'Fr_R'
                'Fr_L'
                'M_R'
                'M_L'
                'Cg_R'
                'Cg_L'
                'S_R'
                'S_L'
                'RS_R'
                'RS_L'
                'V_R'
                'V_L' };
% Get ROIs
[all_ROIs selected_ROIs] = ioi_get_ROIs(job);
if all_ROIs
    % All the ROIs
    nROI = length(IOI.res.ROI);
    % Index vector of ROIs
    roiIndex = 1:nROI;
else
    nROI = numel(selected_ROIs);
    roiIndex = selected_ROIs;
end
% Variable with the source rois names {1 x nROI}
names = newROIname';
% Variable with the target rois names {1 x nROI}
names2 = names;
% Preallocate connectivity values [nROI x nROI x length(job.IOImat)]
Z = zeros([nROI nROI length(job.IOImat)]);
% Go to network results folder
cd(job.results_dir{1});
% Loop over sessions
s1 = 1;
% Initialize progress bar
spm_progress_bar('Init', length(IOI.sess_res{s1}.fname), sprintf('Network analysis'), 'Colors');
% ioi_text_waitbar(0, sprintf('Network analysis\n'));
% Loop over available colors
for c1=1:length(IOI.sess_res{s1}.fname)
    doColor = ioi_doColor(IOI,c1,IC);
    if doColor
        %skip laser - only extract for flow
        if ~(IOI.color.eng(c1)==IOI.color.laser)
            %% Main processing loop
            % Create filename coherent with conn_network.m
            currentFname = fullfile(job.results_dir{1}, sprintf('resultsROI_Condition%02d_%s', s1, colorNames{1+c1}));
            % loop over subjects
            for SubjIdx = 1:length(job.IOImat)
                try
                    %Load IOI.mat information
                    [IOI IOImat dir_ioimat] = ioi_get_IOI(job,SubjIdx);
                    if REPLACE_ROI_NAME
                        IOI.ROIname = newROIname;
                    end
                    if ~isfield(IOI.fcIOS.corr, 'corrOK') % correlation analysis OK
                        fprintf('No correlation analysis available for subject %d of %d... skipping network analysis\n', SubjIdx, length(job.IOImat));
                    else
                        if ~isfield(IOI.fcIOS.corr,'networkOK') || job.force_redo
                            % Check if seed to seed correlation matrix was computed
                            if ~IOI.fcIOS.corr.corrMatrixOK
                                [seed2seedCorrMat seed2seedCorrMatDiff IOI.fcIOS.corr(1).corrMatrixFname IOI.fcIOS.corr(1).corrMatrixDiffFname] = ioi_roi_corr(job, SubjIdx);
                                % Save seed-to-seed correlation data
                                save(IOI.fcIOS.corr(1).corrMatrixFname,'seed2seedCorrMat')
                                % Save seed-to-seed derivatives correlation data
                                save(IOI.fcIOS.corr(1).corrMatrixDiffFname,'seed2seedCorrMatDiff')
                            end
                            % Check if mouse is tratment (1) or control (0)
                            isTreatment(SubjIdx,1) = ~isempty(regexp(IOI.subj_name, [job.treatmentString '[0-9]+'], 'once'));
                            % Load seed to seed correlation matrix
                            load(IOI.fcIOS.corr(1).corrMatrixFname,'seed2seedCorrMat')
                            % Fisher transform
                            seed2seedCorrMat{1}{s1, c1} = fisherz(seed2seedCorrMat{1}{s1, c1});
                            % Replace Inf by NaN
                            seed2seedCorrMat{1}{s1, c1}(~isfinite(seed2seedCorrMat{1}{s1, c1})) = NaN;
                            % Update data in Z
                            Z(:,:,SubjIdx) = seed2seedCorrMat{1}{s1, c1};
                            % Create filename coherent with conn_network.m
                            IOI.fcIOS.corr.networkDataFname{s1, c1} = currentFname;
                            % Network analysis succesful!
                            IOI.fcIOS.corr(1).networkOK = true;
                            % Save IOI matrix
                            save(IOImat,'IOI');
                        end % network OK or redo job
                    end % correlation maps OK
                    out.IOImat{SubjIdx} = IOImat;
                catch exception
                    out.IOImat{SubjIdx} = IOImat;
                    disp(exception.identifier)
                    disp(exception.stack(1))
                end
            end % loop over subjects
            % Save resultsROI_Condition*_Color
            save(IOI.fcIOS.corr.networkDataFname{s1, c1}, ...
                'Z', 'names', 'names2');
            % Process network analyses here
            if job.threshold ~= 0,
            results = conn_network(IOI.fcIOS.corr.networkDataFname{s1, c1}, ...
                roiIndex, job.measures, job.normalType, job.threshold);
            else
                results = conn_network(IOI.fcIOS.corr.networkDataFname{s1, c1}, ...
                roiIndex, job.measures, job.normalType);
            end
            fclose('all')
            varName = sprintf('results_S%02d_%s', s1, colorNames{1+c1});
            controlGroupIdx = find(~isTreatment);
            treatmentGroupIdx = find(isTreatment);
            % Save results in .mat file
            save(fullfile(job.results_dir{1}, [varName '.mat']), 'results', ...
                'controlGroupIdx', 'treatmentGroupIdx');
            % Update progress bar
            spm_progress_bar('Set', c1);
            % ioi_text_waitbar(c1/length(IOI.sess_res{s1}.fname), sprintf('Processing color %d from %d', c1, length(IOI.sess_res{s1}.fname)));
        end
    end
end % colors loop
% Return to working folder
cd(currFolder);
% Clear progress bar
spm_progress_bar('Clear');
% ioi_text_waitbar('Clear');
fprintf('Elapsed time: %s', datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS\n') );

% EOF
