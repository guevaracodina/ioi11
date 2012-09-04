function [corrMatrix corrMatrixFname] = ioi_roi_corr(job,SubjIdx)
% Gets the correlation matrix for every seed/ROI time trace.
%_______________________________________________________________________________
% Copyright (C) 2012 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________

% Get sessions info
[all_sessions selected_sessions] = ioi_get_sessions(job);
%Load IOI.mat information
[IOI IOImat dir_ioimat]= ioi_get_IOI(job,SubjIdx);

if ~isfield(IOI.fcIOS.SPM, 'GLMOK') % GLM OK
    disp(['No GLM regression available for subject ' int2str(SubjIdx) ' ... skipping seed to seed correlation matrix']);
else
    % Get colors to include information
    IC = job.IC;
    [all_ROIs selected_ROIs] = ioi_get_ROIs(job);
    nROI = 1:length(IOI.res.ROI); % All the ROIs
    
    % Load regressed ROI data in cell ROIregress
    load(IOI.fcIOS.SPM.fnameROIregress)
    
    % Loop over sessions
    for s1=1:length(IOI.sess_res)
        if all_sessions || sum(s1==selected_sessions)
            % Loop over available colors
            for c1=1:length(IOI.sess_res{s1}.fname)
                doColor = ioi_doColor(IOI,c1,IC);
                if doColor
                    % Loop over ROI/seeds
                    for r1 = nROI,
                        if IOI.fcIOS.SPM.ROIregressOK{r1}{s1,c1}
                            % Preallocate map for the seed-to-seed correlation matrix
                            tVector = numel(ROIregress{r1}{s1,c1});
                            fprintf('time vector size found for seed %d\n',r1)
                            break; % end loop for as soon as a good ROI is found
                        end
                    end
                    roiMatrix = zeros([tVector numel(nROI)]);
                    % Loop over ROI/seeds
                    for r1 = nROI,
                        if all_ROIs || sum(r1==selected_ROIs)
                            roiMatrix(:, r1) = ROIregress{r1}{s1,c1};
                        end
                    end % loop over sessions
                    % Compute seed-to-seed correlation matrix
                    corrMatrix = corrcoef(roiMatrix);
                    if job.generate_figures
                        figure;
                        if job.save_figures
                            [~, oldName, oldExt] = fileparts(IOI.fcIOS.SPM.fnameROInifti{r1}{s1, c1});
                            newName = [oldName '_fcIOS_map'];
                            % Save as EPS
                            spm_figure('Print', 'Graphics', fullfile(dir_ioimat,newName));
                            % Save as PNG
                            print(h, '-dpng', fullfile(dir_ioimat,newName), '-r300');
                            % Save as nifti
                            ioi_save_nifti(tempCorrMap, fullfile(dir_ioimat,[newName oldExt]), vx);
                            IOI.fcIOS.corr(1).corrMapName{r1}{s1, c1} = fullfile(dir_ioimat,[newName oldExt]);
                        end
                    end
                end
            end % loop over colors
        end
    end % loop over sessions
end % GLM regression ok

% EOF
