function out = ioi_correlation_map_run(job)
% A functional connectivity (fcIOS) map is made by correlating the seed/ROI with
% all other brain (non-masked) pixels
%_______________________________________________________________________________
% Copyright (C) 2012 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________

% Get sessions info
[all_sessions selected_sessions] = ioi_get_sessions(job);

%Big loop over subjects
for SubjIdx=1:length(job.IOImat)
    try
        tic
        clear IOI
        %Load IOI.mat information
        IOImat = job.IOImat{SubjIdx};
        [dir_ioimat dummy] = fileparts(job.IOImat{SubjIdx});
        if isfield(job.IOImatCopyChoice,'IOImatCopy')
            newDir = job.IOImatCopyChoice.IOImatCopy.NewIOIdir;
            newDir = fullfile(dir_ioimat,newDir);
            if ~exist(newDir,'dir'),mkdir(newDir); end
            IOImat = fullfile(newDir,'IOI.mat');
        else
            newDir = dir_ioimat;
        end
        try
            load(IOImat);
        catch
            load(job.IOImat{SubjIdx});
        end
        if ~isfield(IOI.fcIOS.SPM, 'GLMOK') % GLM OK
            disp(['No GLM regression available for subject ' int2str(SubjIdx) ' ... skipping correlation map']);
        else
            if ~isfield(IOI.fcIOS.corr,'corrOK') || job.force_redo
                % Get colors to include information
                IC = job.IC;
                colorNames = fieldnames(IOI.color);
                
                if IOI.res.shrinkageOn
                    vx = [IOI.res.shrink_x IOI.res.shrink_y 1];
                else
                    vx = [1 1 1];
                end
                
                % File name where correlation data is saved
                IOI.fcIOS.corr(1).fname = fullfile(newDir,'seed_based_fcIOS_map.mat');
                
                % Loop over sessions
                for s1=1:length(IOI.sess_res)
                    if all_sessions || sum(s1==selected_sessions)
                        % Loop over available colors
                        for c1=1:length(IOI.sess_res{s1}.fname)
                            doColor = ioi_doColor(IOI,c1,IC);
                            if doColor
                                colorOK = true;
                                %skip laser - only extract for flow
                                if ~(IOI.color.eng(c1)==IOI.color.laser)
                                    %% Do my stuff here!
                                    [all_ROIs selected_ROIs] = ioi_get_ROIs(job);
                                    nROI = 1:length(IOI.res.ROI); % All the ROIs
                                    
                                    % Initialize progress bar
                                    spm_progress_bar('Init', numel(nROI), sprintf('fcIOS map S%d C%d (%s)\n',s1,c1,colorNames{1+c1}), 'Seeds');
                                    
                                    % Loop over ROI/seeds
                                    for r1 = nROI,
                                        if all_ROIs || sum(r1==selected_ROIs)
                                            
                                            % Checking if regression was
                                            % succesful for both the seeds and
                                            % the brain pixels time-courses
                                            if IOI.fcIOS.SPM.wholeImageRegressOK{s1, c1} && IOI.fcIOS.SPM.ROIregressOK{r1}{s1, c1}
                                                fprintf('Loading data, seed %d (%s) session %d C%d (%s)...\n',r1,IOI.ROIname{r1},s1,c1,colorNames{1+c1});
                                                % Load brain pixels time-course
                                                vol = spm_vol(IOI.fcIOS.SPM.fname{s1, c1});
                                                y = spm_read_vols(vol);
                                                % Load ROI time-course
                                                ROIvol = spm_vol(IOI.fcIOS.SPM.fnameROInifti{r1}{s1, c1});
                                                ROI = spm_read_vols(ROIvol);
                                                % Load brain mask
                                                maskVol = spm_vol(IOI.fcIOS.mask.fname);
                                                mask = spm_read_vols(maskVol);
                                                % Load anatomical image
                                                anatVol = spm_vol(IOI.res.file_anat);
                                                anat = spm_read_vols(anatVol);
                                                % Load ROI mask
                                                % ROImaskVol = spm_vol(IOI.res.ROI{r1}.fname);
                                                % ROImask = spm_read_vols(ROImaskVol);
                                                % Preallocate
                                                tempCorrMap = zeros([size(y,1) size(y,2)]);
                                                pValuesMap =  zeros([size(y,1) size(y,2)]);
                                                % Find Pearson's correlation coefficient
                                                fprintf('Computing Pearson''s correlation map...\n');
                                                for iX = 1:size(y,1),
                                                    for iY = 1:size(y,2),
                                                        if mask(iX, iY)
                                                            [tempCorrMap(iX, iY) pValuesMap(iX, iY)]= corr(squeeze(ROI), squeeze(y(iX, iY, 1, :)));
                                                        end
                                                    end
                                                end
                                                % Assign data to be saved to
                                                % .mat file
                                                seed_based_fcIOS_map{r1}{s1,c1}.pearson = tempCorrMap;
                                                seed_based_fcIOS_map{r1}{s1,c1}.pValue = pValuesMap;
                                                
                                                if job.generate_figures
                                                    % Display plots on SPM graphics window
                                                    h = spm_figure('GetWin', 'Graphics');
                                                    spm_figure('Clear', 'Graphics');
                                                    
%                                                     % --------------------------
%                                                     % Compute image ranges to
%                                                     % use a split colormap
%                                                     newMin = -max(tempCorrMap(:))+min(tempCorrMap(:)) - 0.001;
%                                                     newMax = min(tempCorrMap(:));
%                                                     oldMin = min(anat(:));
%                                                     oldMax = max(anat(:));
%                                                     % Value in new range
%                                                     anat = (anat - oldMin)*(newMax-newMin)/(oldMax-oldMin)+newMin;
%                                                     tempCorrMap(mask==0) = anat(mask==0);
%                                                     spm_figure('ColorMap','gray-jet')
%                                                     % -------------------------

                                                    % Improve display
                                                    tempCorrMap(mask==0) = median(tempCorrMap(:));
                                                    spm_figure('ColorMap','jet')
                                                    subplot(211)
                                                    % Correlation map
                                                    imagesc(tempCorrMap); colorbar; axis image
                                                    title(sprintf('fcIOS map Seed %d (%s) S%d C%d (%s)\n',r1,IOI.ROIname{r1},s1,c1,colorNames{1+c1}))
                                                    subplot(212)
                                                    % Show only significant
                                                    % pixels
                                                    imagesc(tempCorrMap .* (pValuesMap <= job.pValue), [-1 1]); colorbar; axis image
                                                    title(sprintf('Significant pixels (p<%f) Seed %d (%s) S%d C%d (%s)\n',job.pValue,r1,IOI.ROIname{r1},s1,c1,colorNames{1+c1}))
                                                    
                                                    if job.save_figures
                                                        [~, oldName, oldExt] = fileparts(IOI.fcIOS.SPM.fnameROInifti{r1}{s1, c1});
                                                        newName = [oldName '_fcIOS_map'];
                                                        % Save as EPS
                                                        spm_figure('Print', 'Graphics', fullfile(newDir,newName));
                                                        % Save as PNG
                                                        print(h, '-dpng', fullfile(newDir,newName), '-r300');
                                                        % Save as nifti
                                                        ioi_save_nifti(tempCorrMap, fullfile(newDir,[newName oldExt]), vx);
                                                        IOI.fcIOS.corr(1).corrMapName{r1}{s1, c1} = fullfile(newDir,[newName oldExt]);
                                                    end
                                                end
                                                
                                                % correlation map succesful!
                                                fprintf('Pearson''s correlation coefficient computed. Seed %d (%s) session %d C%d (%s)\n',r1,IOI.ROIname{r1},s1,c1,colorNames{1+c1});
                                                IOI.fcIOS.corr(1).corrMapOK{r1}{s1, c1} = true;
                                            else
                                                % correlation map failed!
                                                IOI.fcIOS.corr(1).corrMapOK{r1}{s1, c1} = false;
                                                fprintf('Pearson''s correlation coefficient failed! Seed %d (%s) S%d C%d (%s)\n',r1,IOI.ROIname{r1},s1,c1,colorNames{1+c1});
                                            end
                                        end
                                        spm_progress_bar('Set', r1);
                                    end % ROI/seeds loop
                                    % Clear progress bar
                                    spm_progress_bar('Clear');
                                    
                                end
                            end
                        end % colors loop
                    end
                end % Sessions Loop
                % LPF succesful!
                IOI.fcIOS.corr(1).corrOK = true;
                % Save fcIOS data
                save(IOI.fcIOS.corr(1).fname,'seed_based_fcIOS_map')
                % Save IOI matrix
                save(IOImat,'IOI');
            end % LPF OK or redo job
        end % Concentrations OK
        disp(['Elapsed time: ' datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS')]);
        disp(['Subject ' int2str(SubjIdx) ' (' IOI.subj_name ')' ' complete']);
        out.IOImat{SubjIdx} = IOImat;
    catch exception
        out.IOImat{SubjIdx} = IOImat;
        disp(exception.identifier)
        disp(exception.stack(1))
    end
end

% EOF
