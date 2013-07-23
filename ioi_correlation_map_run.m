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
                                    %% Main processing loop
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
                                                
                                                % Preallocate
                                                tempCorrMap = zeros([size(y,1) size(y,2)]);
                                                pValuesMap =  zeros([size(y,1) size(y,2)]);
                                                if isfield (job,'derivative')
                                                    tempCorrMapDiff = zeros([size(y,1) size(y,2)]);
                                                    pValuesMapDiff =  zeros([size(y,1) size(y,2)]);
                                                end
                                                % Find Pearson's correlation coefficient
                                                fprintf('Computing Pearson''s correlation map...\n');
                                                for iX = 1:size(y,1),
                                                    for iY = 1:size(y,2),
                                                        if brainMask(iX, iY)
                                                            [tempCorrMap(iX, iY) pValuesMap(iX, iY)]= corr(squeeze(ROI), squeeze(y(iX, iY, 1, :)));
                                                            if isfield (job,'derivative')
                                                                [tempCorrMapDiff(iX, iY) pValuesMapDiff(iX, iY)]= corr(diff(squeeze(ROI)), diff(squeeze(y(iX, iY, 1, :))));
                                                            end
                                                        end
                                                    end
                                                end
                                                % Assign data to be saved to
                                                % .mat file
                                                seed_based_fcIOS_map{r1}{s1,c1}.pearson = tempCorrMap;
                                                seed_based_fcIOS_map{r1}{s1,c1}.pValue = pValuesMap;
                                                
                                                if isfield (job,'derivative')
                                                    seed_based_fcIOS_map{r1}{s1,c1}.pearsonDiff = tempCorrMapDiff;
                                                    seed_based_fcIOS_map{r1}{s1,c1}.pValueDiff = pValuesMapDiff;
                                                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                    if job.generate_figures
                                                        % Improve display
                                                        tempCorrMapDiff(~brainMask) = median(tempCorrMapDiff(:));
                                                        % Seed annotation dimensions
                                                        % the lower left corner of
                                                        % the bounding rectangle at
                                                        % the point seedX, seedY
                                                        seedX = IOI.res.ROI{r1}.center(2) - IOI.res.ROI{r1}.radius;
                                                        seedY = IOI.res.ROI{r1}.center(1) - IOI.res.ROI{r1}.radius;
                                                        % Seed width
                                                        seedW = 2*IOI.res.ROI{r1}.radius;
                                                        % Seed height
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
                                                        imagesc(tempCorrMapDiff); colorbar; axis image;
                                                        % Display ROI
                                                        rectangle('Position',[seedX seedY seedW seedH],...
                                                            'Curvature',[1,1],...
                                                            'LineWidth',2,'LineStyle','-');
                                                        set(gca,'Xtick',[]); set(gca,'Ytick',[]);
                                                        xlabel('Left', 'FontSize', 14); ylabel('Rostral', 'FontSize', 14);
                                                        title(sprintf('%s fcIOS map Seed %d (%s) S%d C%d (%s) Diff\n',IOI.subj_name,r1,IOI.ROIname{r1},s1,c1,colorNames{1+c1}),'interpreter', 'none', 'FontSize', 14)
                                                        
                                                        % Show only significant
                                                        % pixels
                                                        subplot(212)
                                                        imagesc(tempCorrMapDiff .* (pValuesMapDiff <= job.pValue), [-1 1]); colorbar; axis image;
                                                        % Display ROI
                                                        rectangle('Position',[seedX seedY seedW seedH],...
                                                            'Curvature',[1,1],...
                                                            'LineWidth',2,'LineStyle','-');
                                                        set(gca,'Xtick',[]); set(gca,'Ytick',[]);
                                                        xlabel('Left', 'FontSize', 14); ylabel('Rostral', 'FontSize', 14);
                                                        title(sprintf('%s significant pixels (p<%.2f) Seed %d (%s) S%d C%d (%s) Diff\n',IOI.subj_name,job.pValue,r1,IOI.ROIname{r1},s1,c1,colorNames{1+c1}),'interpreter', 'none', 'FontSize', 14)
                                                        
                                                        if job.save_figures
                                                            try
                                                                if ~isempty(IOI.fcIOS.SPM.fnameROInifti{r1}{s1, c1})
                                                                    [~, oldName, oldExt] = fileparts(IOI.fcIOS.SPM.fnameROInifti{r1}{s1, c1});
                                                                else
                                                                    oldExt = '.nii';
                                                                end
                                                            catch
                                                                oldExt = '.nii';
                                                            end
                                                            % newName = [oldName '_fcIOS_map'];
                                                            newName = [sprintf('%s_R%02d_S%02d_C%d',IOI.subj_name,r1,s1,c1) '_fcIOS_map_diff'];
                                                            dir_corrfigDiff = fullfile(dir_ioimat,'fig_corrMapDiff');
                                                            if ~exist(dir_corrfigDiff,'dir'), mkdir(dir_corrfigDiff); end
                                                            % Save as PNG
                                                            print(h, '-dpng', fullfile(dir_corrfigDiff,newName), '-r300');
                                                            % Save as a figure
                                                            saveas(h, fullfile(dir_corrfigDiff,newName), 'fig');
                                                            % Save as EPS
%                                                             spm_figure('Print', 'Graphics', fullfile(dir_corrfigDiff,newName));
                                                            % Save as nifti
                                                            ioi_save_nifti(tempCorrMap, fullfile(dir_corrfigDiff,[newName oldExt]), vx);
                                                            IOI.fcIOS.corr(1).corrMapNameDiff{r1}{s1, c1} = fullfile(dir_corrfigDiff,[newName oldExt]);
                                                        end
                                                    end
                                                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                    % Derivative processing
                                                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                end
                                                
                                                if job.generate_figures
                                                    % Improve display
                                                    tempCorrMap(~brainMask) = median(tempCorrMap(:));
                                                    % Seed annotation dimensions
                                                    % the lower left corner of
                                                    % the bounding rectangle at
                                                    % the point seedX, seedY
                                                    seedX = IOI.res.ROI{r1}.center(2) - IOI.res.ROI{r1}.radius;
                                                    seedY = IOI.res.ROI{r1}.center(1) - IOI.res.ROI{r1}.radius;
                                                    % Seed width
                                                    seedW = 2*IOI.res.ROI{r1}.radius;
                                                    % Seed height
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
                                                        try
                                                            if ~isempty(IOI.fcIOS.SPM.fnameROInifti{r1}{s1, c1})
                                                                [~, oldName, oldExt] = fileparts(IOI.fcIOS.SPM.fnameROInifti{r1}{s1, c1});
                                                            else
                                                                oldExt = '.nii';
                                                            end
                                                        catch
                                                            oldExt = '.nii';
                                                        end
                                                        % newName = [oldName '_fcIOS_map'];
                                                        newName = [sprintf('%s_R%02d_S%02d_C%d',IOI.subj_name,r1,s1,c1) '_fcIOS_map'];
                                                        if isfield(job.IOImatCopyChoice,'IOImatCopy')
                                                            dir_corrfig = fullfile(dir_ioimat,strcat('fig_',job.IOImatCopyChoice.IOImatCopy.NewIOIdir));
                                                        else
                                                            dir_corrfig = fullfile(dir_ioimat,'fig_corrMap');
                                                        end
                                                        if ~exist(dir_corrfig,'dir'), mkdir(dir_corrfig); end
                                                        % Save as PNG
                                                        print(h, '-dpng', fullfile(dir_corrfig,newName), '-r300');
                                                        % Save as a figure
                                                        saveas(h, fullfile(dir_corrfig,newName), 'fig');
                                                        % Save as EPS
%                                                         spm_figure('Print', 'Graphics', fullfile(dir_corrfig,newName));
                                                        % Save as nifti
                                                        ioi_save_nifti(tempCorrMap, fullfile(dir_corrfig,[newName oldExt]), vx);
                                                        IOI.fcIOS.corr(1).corrMapName{r1}{s1, c1} = fullfile(dir_corrfig,[newName oldExt]);
                                                    end
                                                end
                                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                
                                                % Convert Pearson's correlation
                                                % coeff. r to Fisher's z value.
                                                if job.fisherZ
                                                    % Find Pearson's correlation coefficient
                                                    fprintf('Computing Fisher''s Z statistic...\n');
                                                    try
                                                        if ~isempty(IOI.fcIOS.SPM.fnameROInifti{r1}{s1, c1})
                                                            [~, oldName, oldExt] = fileparts(IOI.fcIOS.SPM.fnameROInifti{r1}{s1, c1});
                                                        else
                                                            oldExt = '.nii';
                                                        end
                                                    catch
                                                        oldExt = '.nii';
                                                    end
                                                    % newName = [oldName '_fcIOS_Zmap'];
                                                    newName = [sprintf('%s_R%02d_S%02d_C%d',IOI.subj_name,r1,s1,c1) '_fcIOS_Zmap'];
                                                    dir_fisherZfig = fullfile(dir_ioimat,'fig_fisherZ');
                                                    if ~exist(dir_fisherZfig,'dir'), mkdir(dir_fisherZfig); end
                                                    % Convert r to z, masking
                                                    % out non-brain voxels
                                                    seed_based_fcIOS_map{r1}{s1,c1}.fisher = fisherz(tempCorrMap) .* brainMask;
                                                    
                                                    
                                                    if job.generate_figures
                                                        % Save as nifti
                                                        ioi_save_nifti(seed_based_fcIOS_map{r1}{s1,c1}.fisher, fullfile(dir_fisherZfig,[newName oldExt]), vx);
                                                        % Save masked file as nifti
                                                        ioi_save_nifti(seed_based_fcIOS_map{r1}{s1,c1}.fisher .* (pValuesMap < job.pValue), fullfile(dir_fisherZfig,[newName '_Masked' oldExt]), vx);
                                                        IOI.fcIOS.corr(1).zMapName{r1}{s1, c1} = fullfile(dir_fisherZfig,[newName oldExt]);
                                                        IOI.fcIOS.corr(1).zMapNameMask{r1}{s1, c1} = fullfile(dir_fisherZfig,[newName '_Masked' oldExt]);
                                                        % Find 1st positive value
                                                        z1 = sort(seed_based_fcIOS_map{r1}{s1,c1}.fisher(:),1,'ascend');
                                                        idx1  = find(z1>0, 1, 'first');
                                                        minPosVal = z1(idx1);
                                                        
                                                        % Find 1st negative value
                                                        z2 = sort(seed_based_fcIOS_map{r1}{s1,c1}.fisher(:),1,'descend');
                                                        idx2  = find(z2<0, 1, 'first');
                                                        minNegVal = z2(idx2);
                                                        
                                                        % Get parameters for overlay
                                                        anatomical      = IOI.res.file_anat;
                                                        % positiveMap     = IOI.fcIOS.corr.zMapName{r1}{s1, c1};
                                                        positiveMap     = IOI.fcIOS.corr.zMapNameMask{r1}{s1, c1};
                                                        % negativeMap     = IOI.fcIOS.corr.zMapName{r1}{s1, c1};
                                                        negativeMap     = IOI.fcIOS.corr.zMapNameMask{r1}{s1, c1};
                                                        colorNames      = fieldnames(IOI.color);
                                                        mapRange        = {[minPosVal max(z1)], [minNegVal min(z2)]};
                                                        titleString     = sprintf('%s seed%d S%d(%s)Z-map',IOI.subj_name,r1,s1,colorNames{1+c1});
                                                        % Display plots on SPM graphics window
                                                        h = spm_figure('GetWin', 'Graphics');
                                                        spm_figure('Clear', 'Graphics');
                                                        h = ioi_overlay_map(anatomical, positiveMap, negativeMap, mapRange, titleString);
                                                        if job.save_figures
                                                            % Save as PNG
                                                            print(h, '-dpng', fullfile(dir_fisherZfig,newName), '-r150');
                                                            % Save as a figure
                                                            saveas(h, fullfile(dir_fisherZfig,newName), 'fig');
                                                            % Save as EPS
%                                                             spm_figure('Print', 'Graphics', fullfile(dir_fisherZfig,newName));
                                                        end
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
                % correlation succesful!
                IOI.fcIOS.corr(1).corrOK = true;
                % Compute seed to seed correlation matrix
                if job.seed2seedCorrMat
                    [seed2seedCorrMat seed2seedCorrMatDiff ...
                        IOI.fcIOS.corr(1).corrMatrixFname ...
                        IOI.fcIOS.corr(1).corrMatrixDiffFname ...
                        pMatrixMat pMatrixMatDiff] = ioi_roi_corr(job, SubjIdx);
                    % seed-to-seed correlation succesful!
                    IOI.fcIOS.corr(1).corrMatrixOK = true;
                    % Save seed-to-seed correlation data
                    save(IOI.fcIOS.corr(1).corrMatrixFname,'seed2seedCorrMat', 'pMatrixMat')
                    % Save seed-to-seed derivatives correlation data
                    save(IOI.fcIOS.corr(1).corrMatrixDiffFname,'seed2seedCorrMatDiff', 'pMatrixMatDiff')
                end
                % Save fcIOS data
                save(IOI.fcIOS.corr(1).fname,'seed_based_fcIOS_map')
                % Save IOI matrix
                save(IOImat,'IOI');
            end % correlation OK or redo job
        end % GLM OK
        disp(['Elapsed time: ' datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS')]);
        disp(['Subject ' int2str(SubjIdx) ' (' IOI.subj_name ')' ' complete']);
        out.IOImat{SubjIdx} = IOImat;
    catch exception
        out.IOImat{SubjIdx} = IOImat;
        disp(exception.identifier)
        disp(exception.stack(1))
    end
end % Big loop over subjects

% EOF
