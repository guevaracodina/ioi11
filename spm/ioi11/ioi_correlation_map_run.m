function out = ioi_correlation_map_run(job)
% A functional connectivity (fcIOS) map is made by correlating the seed/ROI with
% all other brain (non-masked) pixels
%_______________________________________________________________________________
% Copyright (C) 2012 LIOM Laboratoire d'Imagerie Optique et Mol�culaire
%                    �cole Polytechnique de Montr�al
%_______________________________________________________________________________

% Get sessions info
[all_sessions selected_sessions] = ioi_get_sessions(job);

%Big loop over subjects
for SubjIdx=1:length(job.IOImat)
    try
        tic
        %Load IOI.mat information
        [IOI IOImat dir_ioimat]= ioi_get_IOI(job,SubjIdx);  

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
                IOI.fcIOS.corr(1).fname = fullfile(dir_ioimat,'seed_based_fcIOS_map.mat');
                
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
                                    % Load regressed ROI data in cell ROIregress
                                    load(IOI.fcIOS.SPM.fnameROIregress)
                                    % Loop over ROI/seeds
                                    for r1 = nROI,
                                        if all_ROIs || sum(r1==selected_ROIs)
                                            
                                            % Checking if regression was
                                            % succesful for both the seeds and
                                            % the brain pixels time-courses
                                            if IOI.fcIOS.SPM.wholeImageRegressOK{s1, c1} && IOI.fcIOS.SPM.ROIregressOK{r1}{s1, c1}
                                                fprintf('Loading data, seed %d (%s) session %d C%d (%s)...\n',r1,IOI.ROIname{r1},s1,c1,colorNames{1+c1});
                                                % Load brain pixels time-course
                                                % already filtered/downsampled &
                                                % regressed
                                                vol = spm_vol(IOI.fcIOS.SPM.fname{s1, c1});
                                                y = spm_read_vols(vol);
                                                % Load ROI time-course already
                                                % filtered/downsampled &
                                                % regressed (column vector)
                                                % ROIvol = spm_vol(IOI.fcIOS.SPM.fnameROInifti{r1}{s1, c1});
                                                % ROI = spm_read_vols(ROIvol);
                                                ROI = ROIregress{r1}{s1, c1}';
                                                % Load brain mask
                                                brainMaskVol = spm_vol(IOI.fcIOS.mask.fname);
                                                brainMask = logical(spm_read_vols(brainMaskVol));
                                                if size(brainMask,1)~= size(y,1)|| size(brainMask,2)~= size(y,2)
                                                    brainMask = ioi_MYimresize(brainMask, [size(y,1) size(y,2)]);
                                                end
                                                % Load anatomical image
                                                % anatVol = spm_vol(IOI.res.file_anat);
                                                % anat = spm_read_vols(anatVol);
                                                % if size(anat,1)~= size(y,1)|| size(anat,2)~= size(y,2)
                                                %     anat = ioi_MYimresize(anat, [size(y,1) size(y,2)]);
                                                % end
                                                % Load ROI mask
                                                % ROImaskVol = spm_vol(IOI.res.ROI{r1}.fname);
                                                % ROImask = spm_read_vols(ROImaskVol);
                                                % if size(ROImask,1)~= size(y,1)|| size(ROImask,2)~= size(y,2)
                                                %     ROImask = ioi_MYimresize(ROImask, [size(y,1) size(y,2)]);
                                                % end
                                                % Preallocate
                                                tempCorrMap = zeros([size(y,1) size(y,2)]);
                                                pValuesMap =  zeros([size(y,1) size(y,2)]);
                                                % Find Pearson's correlation coefficient
                                                fprintf('Computing Pearson''s correlation map...\n');
                                                for iX = 1:size(y,1),
                                                    for iY = 1:size(y,2),
                                                        if brainMask(iX, iY)
                                                            [tempCorrMap(iX, iY) pValuesMap(iX, iY)]= corr(squeeze(ROI), squeeze(y(iX, iY, 1, :)));
                                                        end
                                                    end
                                                end
                                                % Assign data to be saved to
                                                % .mat file
                                                seed_based_fcIOS_map{r1}{s1,c1}.pearson = tempCorrMap;
                                                seed_based_fcIOS_map{r1}{s1,c1}.pValue = pValuesMap;
                                                
                                                if job.generate_figures
                                                    
%                                                     % --------------------------
%                                                     % Compute image ranges to
%                                                     % use a split colormap
%                                                     newMin = -max(tempCorrMap(:))+min(tempCorrMap(:)) - 0.001;
%                                                     newMax = min(tempCorrMap(:));
%                                                     oldMin = min(anat(:));
%                                                     oldMax = max(anat(:));
%                                                     % Value in new range
%                                                     anat = (anat - oldMin)*(newMax-newMin)/(oldMax-oldMin)+newMin;
%                                                     tempCorrMap(brainMask==0) = anat(brainMask==0);
%                                                     spm_figure('ColorMap','gray-jet')
%                                                     % -------------------------

                                                    % Improve display
                                                    tempCorrMap(~brainMask) = median(tempCorrMap(:));
                                                    % Seed annotation dimensions
                                                    % the lower left corner of
                                                    % the bounding rectangle at
                                                    % the point seedX, seedY
                                                    seedX = IOI.res.ROI{r1}.center(2) - IOI.res.ROI{r1}.radius;
                                                    seedY = IOI.res.ROI{r1}.center(1) - IOI.res.ROI{r1}.radius;
                                                    seedW = 2*IOI.res.ROI{r1}.radius;
                                                    seedH = 2*IOI.res.ROI{r1}.radius;
                                                    if isfield(IOI.res,'shrinkageOn')
                                                        if IOI.res.shrinkageOn == 1
                                                            seedX = seedX / IOI.res.shrink_x;
                                                            seedY = seedY / IOI.res.shrink_y;
                                                            seedW = seedW / IOI.res.shrink_x;
                                                            seedH = seedH / IOI.res.shrink_y;
                                                        end
                                                    end
                                                    
                                                    % Display plots on SPM graphics window
                                                    h = spm_figure('GetWin', 'Graphics');
                                                    spm_figure('Clear', 'Graphics');
                                                    spm_figure('ColorMap','jet')
                                                    
                                                    % Correlation map
                                                    subplot(211)
                                                    imagesc(tempCorrMap); colorbar; axis image; 
                                                    % Display ROI
                                                    rectangle('Position',[seedX seedY seedW seedH],...
                                                        'Curvature',[1,1],...
                                                        'LineWidth',2,'LineStyle','-');
                                                    set(gca,'Xtick',[]); set(gca,'Ytick',[]);
                                                    xlabel('Left', 'FontSize', 14); ylabel('Rostral', 'FontSize', 14);
                                                    title(sprintf('%s fcIOS map Seed %d (%s) S%d C%d (%s)\n',IOI.subj_name,r1,IOI.ROIname{r1},s1,c1,colorNames{1+c1}),'interpreter', 'none', 'FontSize', 14)
                                                    
                                                    % Show only significant
                                                    % pixels
                                                    subplot(212)
                                                    imagesc(tempCorrMap .* (pValuesMap <= job.pValue), [-1 1]); colorbar; axis image;
                                                    % Display ROI
                                                    rectangle('Position',[seedX seedY seedW seedH],...
                                                        'Curvature',[1,1],...
                                                        'LineWidth',2,'LineStyle','-');
                                                    set(gca,'Xtick',[]); set(gca,'Ytick',[]);
                                                    xlabel('Left', 'FontSize', 14); ylabel('Rostral', 'FontSize', 14);
                                                    title(sprintf('%s significant pixels (p<%.2f) Seed %d (%s) S%d C%d (%s)\n',IOI.subj_name,job.pValue,r1,IOI.ROIname{r1},s1,c1,colorNames{1+c1}),'interpreter', 'none', 'FontSize', 14)
                                                    
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
                                                
                                                % correlation map succesful!
                                                fprintf('Pearson''s correlation coefficient computed. Seed %d (%s) session %d C%d (%s)\n',r1,IOI.ROIname{r1},s1,c1,colorNames{1+c1});
                                                IOI.fcIOS.corr(1).corrMapOK{r1}{s1, c1} = true;
                                            else
                                                % correlation map failed!
                                                IOI.fcIOS.corr(1).corrMapOK{r1}{s1, c1} = false;
                                                fprintf('Pearson''s correlation coefficient failed! Seed %d (%s) S%d C%d (%s)\n',r1,IOI.ROIname{r1},s1,c1,colorNames{1+c1});
                                            end
                                        end
                                        % Update progress bar
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