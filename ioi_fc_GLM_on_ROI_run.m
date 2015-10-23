function out = ioi_fc_GLM_on_ROI_run(job)
% GLM regression of global brain signal in resting-state from ROI/seeds time
% trace in order to remove global source of variance.
%_______________________________________________________________________________
% Copyright (C) 2012 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________

% Get sessions and ROI info
[all_sessions selected_sessions] = ioi_get_sessions(job);
[all_ROIs selected_ROIs] = ioi_get_ROIs(job);

%Big loop over subjects
for SubjIdx=1:length(job.IOImat)
    try
        eTime = tic;
        clear ROI SPM
        %Load IOI.mat information
        [IOI IOImat dir_ioimat]= ioi_get_IOI(job,SubjIdx);
        
        if ~isfield(IOI.res,'ROIOK')
            disp(['No ROI available for subject ' int2str(SubjIdx) ' ... skipping series extraction']);
        else
            if ~isfield(IOI.res,'seriesOK') || ~isfield(IOI.fcIOS.mask,'seriesOK')
                disp(['Extracted series not available for subject ' int2str(SubjIdx) ' ... skipping GLM']);
            else
                if ~isfield(IOI.fcIOS.filtNdown,'filtNdownOK')
                    disp(['Filtered/downsampled time-courses not available for subject ' int2str(SubjIdx) ' ... skipping GLM']);
                else
                    if ~isfield(IOI.fcIOS,'SPM')
                        % It is the first time to run SPM
                        IOI.fcIOS.SPM = struct([]);
                    end
                    if ~isfield(IOI.fcIOS.SPM,'GLMOK') || job.force_redo
                        % Get colors to include information
                        IC = job.IC;
                        colorNames = fieldnames(IOI.color);
                        % Load filtered downsampled signals
                        filtNdownData = load(IOI.fcIOS.filtNdown.fname);
                        fnameROIregress = fullfile(dir_ioimat,'ROIregress.mat');
                        % Loop over sessions
                        for s1=1:length(IOI.sess_res)
                            if all_sessions || sum(s1==selected_sessions)
                                sessionDir = [dir_ioimat filesep 'S' sprintf('%02d',s1)];
                                if ~exist(sessionDir,'dir'),mkdir(sessionDir); end
                                % Loop over available colors
                                for c1=1:length(IOI.sess_res{s1}.fname)
                                    doColor = ioi_doColor(IOI,c1,IC);
                                    if doColor
                                        %skip laser - only extract for flow
                                        if ~(IOI.color.eng(c1)==IOI.color.laser)
                                            colorDir = [sessionDir filesep 'C' sprintf('%d',c1)];
                                            if ~exist(colorDir,'dir'),mkdir(colorDir); end
                                            % Get filtered downsampled signals
                                            brainSignal = filtNdownData.filtNdownBrain{1}{s1, c1};
                                            % Initialize single voxel 4-D series
                                            brainSignalRep = zeros([1 1 1 numel(brainSignal)]);
                                            dim = [1 1 1/IOI.fcIOS.filtNdown.downFreq];
                                            if job.wholeImage
                                                %% GLM on images here!
                                                % y = ioi_get_images(IOI,1:IOI.sess_res{s1}.n_frames,c1,s1,dir_ioimat,shrinkage_choice);
                                                
                                                % Constructing inputs required for
                                                % GLM analysis within the SPM
                                                % framework
                                                clear SPM
                                                SPM.xY.VY = spm_vol(IOI.fcIOS.filtNdown.fnameWholeImage{s1, c1});
                                                y = spm_read_vols(SPM.xY.VY);
                                                % Preallocating output images
                                                yRegress = zeros(size(y));
                                                
                                                % All regressors are identified
                                                % here, lets take the mean global
                                                % signal just to test
                                                SPM.xX.name = cellstr(['Global Brain Signal']);
                                                SPM.xX.X = brainSignal';        % Regression is along first dimension. For one regressor it is a column vector.
                                                if job.heartRate
                                                    load(IOI.fcIOS.SPM.physioHeartFile{1})
                                                    heartRate = physio.heartRate(1:numel(brainSignal));
                                                    SPM.xX.name = {'Global Brain Signal' 'Heart Rate'};
                                                    selectedRegressorsArray = [brainSignal' heartRate];
                                                    SPM.xX.X = selectedRegressorsArray;        % Regression is along first dimension. For one regressor it is a column vector.
                                                end
                                                
                                                % A revoir
                                                SPM.xX.iG = [];
                                                SPM.xX.iH = [];
                                                SPM.xX.iC = 1:size(SPM.xX,2);   % Indices of regressors of interest
                                                SPM.xX.iB = [];                 % Indices of confound regressors
                                                SPM.xGX.rg = [];                % Raw globals, need to check this //EGC
                                                
                                                SPM.xVi.Vi = {speye(size(SPM.xX.X,1))}; % Time correlation
                                                
                                                % Directory to save SPM and img/hdr
                                                % files
                                                SPM.swd = colorDir;
                                                if ~exist(SPM.swd,'dir'),mkdir(SPM.swd); end
                                                fprintf('\nPerforming GLM on %s whole images, Session %d Color %d (%s)...\n',IOI.subj_name,s1,c1,colorNames{1+c1})
                                                try
                                                    % GLM is performed here
                                                    % at whole image level.
                                                    % Orlando!
                                                    SPM = spm_spm(SPM);
                                                    
                                                    if job.regressBrainSignal == 1,
                                                        % EGC: Remember to
                                                        % uncomment this
                                                        % section to make a
                                                        % cleaner code
                                                        % Subtract global brain signal
                                                        % from every pixel time course
%                                                         betaVol = spm_vol(fullfile(SPM.swd,SPM.Vbeta.fname));
%                                                         beta = spm_read_vols(betaVol);
                                                        % Create a single voxel 4-D
                                                        % series
%                                                         brainSignalRep(1,1,1,:) = brainSignal;
%                                                         yRegress = y - repmat(beta,[1 1 1 size(y,4)]) .* repmat(brainSignalRep,[size(beta,1) size(beta,2) 1 1]);
                                                        if job.heartRate
                                                            beta = [];
                                                            for iRegressors = 1:size(SPM.xX.X,2)
                                                                betaVol{iRegressors} = spm_vol(fullfile(SPM.swd,SPM.Vbeta(iRegressors).fname));
                                                                betaCell{iRegressors} = spm_read_vols(betaVol{iRegressors});
                                                                beta = [beta, betaCell{iRegressors}(:)];
                                                            end
                                                            signaltoRegress = beta * selectedRegressorsArray';
                                                            % Create a single voxel 4-D series
                                                            signaltoRegress = reshape(signaltoRegress,[size(y,1) size(y,2) 1 size(y,4)]);
                                                            % Subtract ROIs signal from every pixel time course
                                                            yRegress = y - signaltoRegress;
%                                                             filtNdownfnameRegress = fullfile(sessionDir,[IOI.subj_name '_OD_' IOI.color.eng(c1) '_regress_' sprintf('%05d',1) 'to' sprintf('%05d',IOI.sess_res{s1}.n_frames) '.nii']);
% 
%                                                             % Save NIFTI file
%                                                             ioi_save_nifti(yRegress, filtNdownfnameRegress, dim);
%                                                             % Brain signal regression succesful!
%                                                             IOI.fcIOS.SPM(1).wholeImageRegressOK{s1, c1} = true;
                                                            fprintf('\nHeart rate regressed from %s whole images in Session %d Color %d (%s) done!\n',IOI.subj_name,s1,c1,colorNames{1+c1})
                                                        end
                                                        filtNdownfnameRegress = fullfile(sessionDir,[IOI.subj_name '_OD_' IOI.color.eng(c1) '_regress_' sprintf('%05d',1) 'to' sprintf('%05d',IOI.sess_res{s1}.n_frames) '.nii']);
                                                        % Save NIFTI file
                                                        ioi_save_nifti(yRegress, filtNdownfnameRegress, dim);
                                                        % Brain signal regression succesful!
                                                        IOI.fcIOS.SPM(1).wholeImageRegressOK{s1, c1} = true;
                                                        fprintf('\nGlobal brain signal regressed from %s whole images in Session %d Color %d (%s) done!\n',IOI.subj_name,s1,c1,colorNames{1+c1})

                                                        
                                                    else
                                                        %% Just copy the filtered signal
                                                        yRegress = y;
                                                        filtNdownfnameRegress = fullfile(sessionDir,[IOI.subj_name '_OD_' IOI.color.eng(c1) '_regress_' sprintf('%05d',1) 'to' sprintf('%05d',IOI.sess_res{s1}.n_frames) '.nii']);
                                                        % Save NIFTI file
                                                        ioi_save_nifti(yRegress, filtNdownfnameRegress, dim);
                                                        % Brain signal regression succesful!
                                                        IOI.fcIOS.SPM(1).wholeImageRegressOK{s1, c1} = true;
                                                        fprintf('\nGlobal brain signal NOT regressed from %s whole images in Session %d Color %d (%s) done!\n',IOI.subj_name,s1,c1,colorNames{1+c1})
                                                    end
                                                    
                                                    % Update SPM matrix info
                                                    IOI.fcIOS.SPM(1).fnameSPM{s1, c1} = SPM.swd;
                                                    IOI.fcIOS.SPM(1).fname{s1, c1} = filtNdownfnameRegress;
                                                catch exception
                                                    % Brain signal regression failed!
                                                    IOI.fcIOS.SPM(1).wholeImageRegressOK{s1, c1} = false;
                                                    fprintf('\nGlobal brain signal regressed from %s whole images in Session %d Color %d (%s) failed!\n',IOI.subj_name,s1,c1,colorNames{1+c1})
                                                    disp(exception.identifier)
                                                    disp(exception.stack(1))
                                                end
                                            end % end on GLM on images
                                            
                                            % Loop over ROIs
                                            for r1=1:length(IOI.res.ROI)
                                                if all_ROIs || sum(r1==selected_ROIs)
                                                    ROIdir = [colorDir filesep 'ROI' sprintf('%02d',r1)];
                                                    if ~exist(ROIdir,'dir'),mkdir(ROIdir); end
                                                    % Initialize y tilde (ROIregress)
                                                    ROIregress{r1}{s1,c1} = [];
                                                    %% GLM on ROI code
                                                    y = filtNdownData.filtNdownROI{r1}{s1, c1};
                                                    % Initialize single voxel
                                                    % 4-D series
                                                    y2 = zeros([1 1 1 numel(y)]);
                                                    
                                                    if job.generate_figures
                                                        % Display plots on SPM graphics window
                                                        spm_figure('GetWin', 'Graphics');
                                                        spm_figure('Clear', 'Graphics');
                                                        subplot(311); plot(y);
                                                        title(sprintf('Seed %d time-course, S%d, C%d (%s)',r1,s1,c1,colorNames{1+c1}),'FontSize',14);
                                                        subplot(312); plot(brainSignal);
                                                        title(sprintf('Mean global signal time-course, S%d, C%d (%s)',s1,c1,colorNames{1+c1}),'FontSize',14);
                                                    end
                                                    
                                                    % Creating nifti files to be able to use SPM later
                                                    % --------------------------
                                                    fnameNIFTI = fullfile(ROIdir,['ROI' sprintf('%02d',r1) '_S' sprintf('%02d',s1) '_C' num2str(c1),'.nii']);
                                                    dim = [1 1 1];
                                                    % Create a single voxel 4-D
                                                    % series
                                                    y2(1,1,1,:) = y;
                                                    ioi_save_nifti(y2, fnameNIFTI, dim);
                                                    % --------------------------
                                                    % end of nifti processing
                                                    
                                                    % Constructing inputs
                                                    % required for GLM analysis
                                                    % within the SPM framework
                                                    clear SPM
                                                    SPM.xY.VY = spm_vol(fnameNIFTI);
                                                    
                                                    % All regressors are
                                                    % identified here, lets take
                                                    % the mean global signal
                                                    % just to test
                                                    SPM.xX.name = cellstr(['Global Brain Signal']);
                                                    % Orlando: Add heart
                                                    % rate signal
                                                    SPM.xX.X = brainSignal';        % Regression is along first dimension. For one regressor it is a column vector.
                                                    
                                                    % A revoir
                                                    SPM.xX.iG = [];
                                                    SPM.xX.iH = [];
                                                    SPM.xX.iC = 1:size(SPM.xX,2);   % Indices of regressors of interest
                                                    SPM.xX.iB = [];                 % Indices of confound regressors
                                                    SPM.xGX.rg = [];                % Raw globals, need to check this //EGC
                                                    
                                                    SPM.xVi.Vi = {speye(size(SPM.xX.X,1))}; % Time correlation
                                                    
                                                    % Directory to save SPM and
                                                    % img/hdr files
                                                    SPM.swd = ROIdir;
                                                    if ~exist(SPM.swd,'dir'),mkdir(SPM.swd); end
                                                    fprintf('\nPerforming GLM for %s ROI %d Session %d Color %d (%s)...\n',IOI.subj_name,r1,s1,c1,colorNames{1+c1})
                                                    try
                                                        % GLM is performed here
                                                        SPM = spm_spm(SPM);
                                                        
                                                        if job.regressBrainSignal == 1,
                                                            % Subtract global brain
                                                            % signal from ROI time
                                                            % courses
                                                            betaVol = spm_vol(fullfile(SPM.swd,SPM.Vbeta.fname));
                                                            beta = spm_read_vols(betaVol);
                                                            % Orlando: Add
                                                            % heart rate
                                                            % signal
                                                            ROIregress{r1}{s1, c1} = y - beta * brainSignal;
                                                            
                                                            % Brain signal regression succesful!
                                                            IOI.fcIOS.SPM(1).ROIregressOK{r1}{s1, c1} = true;
                                                            
                                                            % Identify in IOI the file name of the time series
                                                            IOI.fcIOS.SPM(1).fnameROIregress = fnameROIregress;
                                                            fprintf('\nGlobal brain signal regressed from %s ROI %d (%s) Session %d Color %d (%s) done!\n',IOI.subj_name,r1,IOI.ROIname{r1},s1,c1,colorNames{1+c1})
                                                            if job.generate_figures
                                                                spm_figure('GetWin', 'Graphics');
                                                                subplot(313); plot(ROIregress{r1}{s1, c1});
                                                                title(sprintf('Global signal regressed from ROI time-course %d, S%d, C%d (%s)',r1,s1,c1,colorNames{1+c1}),'FontSize',14);
                                                            end
                                                        else
                                                            %% Just copy the filtered signal
                                                            ROIregress{r1}{s1, c1} = y;
                                                            % Brain signal regression succesful!
                                                            IOI.fcIOS.SPM(1).ROIregressOK{r1}{s1, c1} = true;
                                                            
                                                            % Identify in IOI the file name of the time series
                                                            IOI.fcIOS.SPM(1).fnameROIregress = fnameROIregress;
                                                            fprintf('\nGlobal brain signal NOT regressed from %s ROI %d (%s) Session %d Color %d (%s) done!\n',IOI.subj_name,r1,IOI.ROIname{r1},s1,c1,colorNames{1+c1})
                                                            if job.generate_figures
                                                                spm_figure('GetWin', 'Graphics');
                                                                subplot(313); plot(ROIregress{r1}{s1, c1});
                                                                title(sprintf('Global signal NOT regressed from ROI time-course %d, S%d, C%d (%s)',r1,s1,c1,colorNames{1+c1}),'FontSize',14);
                                                            end
                                                        end
                                                        
                                                        % Update SPM matrix info
                                                        IOI.fcIOS.SPM(1).fnameROISPM{r1}{s1, c1} = SPM.swd;
                                                        IOI.fcIOS.SPM(1).fnameROInifti{r1}{s1, c1} = fnameNIFTI;
                                                        % Contrasts
                                                        % --------------------------
                                                        % [Ic, xCon] = spm_conman(SPM, 'T|F', Inf, 'Select contrasts...', 'Contrasts amongst spectral regressors', 1);
                                                        % SPM.xCon = xCon;
                                                        % [SPM,xSPM]=spm_getSPM(SPM,[]);
                                                        % --------------------------
                                                        fprintf('\nGLM for %s ROI %d Session %d Color %d (%s) done!\n',IOI.subj_name,r1,s1,c1,colorNames{1+c1})
                                                    catch exception
                                                        % Brain signal regression on ROI failed!
                                                        IOI.fcIOS.SPM(1).ROIregressOK{r1}{s1, c1} = false;
                                                        fprintf('\nGLM for %s ROI %d Session %d Color %d (%s) failed!\n',IOI.subj_name,r1,s1,c1,colorNames{1+c1})
                                                        disp(exception.identifier)
                                                        disp(exception.stack(1))
                                                    end
                                                end
                                            end % ROI loop
                                        end
                                    end 
                                end % colors loop
                            end
                        end % sessions loop
                        
                        %%
                        % ------------------------------------------------------
                        % If ROI regression failed, then get ROI time course
                        % from the whole image regressed series.
                        % ------------------------------------------------------
                        % Identify in IOI the file name of the time series
                        IOI.fcIOS.SPM(1).fnameROIregress = fnameROIregress;
                        % Get mask for each ROI
                        [~, mask] = ioi_get_ROImask(IOI,job);
                        Amask = []; % Initialize activation mask
                        % We are not extracting brain mask here
                        job.extractingBrainMask = false;
                        job.extractBrainMask = false;
                        % Extract ROI from regressed whole image series.
                        [ROIregressTmp IOI] = ...
                            ioi_extract_core(IOI,job,mask,Amask,'regressData');
                        % Loop over sessions
                        for s1=1:length(IOI.sess_res)
                            if all_sessions || sum(s1==selected_sessions)
                                % Loop over available colors
                                for c1=1:length(IOI.sess_res{s1}.fname)
                                    doColor = ioi_doColor(IOI,c1,IC);
                                    if doColor
                                        %skip laser - only extract for flow
                                        if ~(IOI.color.eng(c1)==IOI.color.laser)
                                            % Loop over ROIs
                                            for r1=1:length(IOI.res.ROI)
                                                if all_ROIs || sum(r1==selected_ROIs)
                                                    if ~IOI.fcIOS.SPM(1).ROIregressOK{r1}{s1, c1}
                                                        % GLM on ROI not succesful, so
                                                        % we write the ROI extracted frm
                                                        % the regressed whole image.
                                                        ROIregress{r1}{s1, c1} =  ROIregressTmp{r1}{s1, c1};
                                                        % Brain signal regression succesful!
                                                        IOI.fcIOS.SPM(1).ROIregressOK{r1}{s1, c1} = true;
                                                    end
                                                end
                                            end % ROI loop
                                        end
                                    end
                                end % colors loop
                            end
                        end % sessions loop
                        % ------------------------------------------------------
                        
                        %% GLM regression succesful!
                        IOI.fcIOS.SPM(1).GLMOK = true;
                        % Save IOI.fcIOS.SPM.fnameROIregress anyway!
                        %if job.regressBrainSignal == 1,
                            save(fnameROIregress,'ROIregress');
                        %end
                        if job.cleanupGLM
                            % Keeps only NIfTI files of succesfully regressed ROIs
                            IOI.fcIOS.SPM.cleanupOK = ioi_fc_GLM_on_ROI_cleanup(IOI, job);
                        end
                        % Save IOI matrix
                        save(IOImat,'IOI');
                    end % GLM OK or redo job
                end % Filtering&Downsampling OK
            end % Time-series OK
        end % ROI OK
        disp(['Elapsed time: ' datestr(datenum(0,0,0,0,0,toc(eTime)),'HH:MM:SS')]);
        disp(['Subject ' int2str(SubjIdx) ' (' IOI.subj_name ')' ' complete']);
        out.IOImat{SubjIdx} = IOImat;
        cd(spm('Dir'));     % Return to SPM working directory
    catch exception
        out.IOImat{SubjIdx} = IOImat;
        disp(exception.identifier)
        disp(exception.stack(1))
    end
end % Big loop over subjects

% EOF

