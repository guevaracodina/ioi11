function out = ioi_fc_GLM_on_ROI_run(job)
% GLM regression of global brain signal in resting-state from ROI/seeds time
% trace in order to remove global source of variance.
%_______________________________________________________________________________
% Copyright (C) 2012 LIOM Laboratoire d'Imagerie Optique et Mol�culaire
%                    �cole Polytechnique de Montr�al
%_______________________________________________________________________________

% Get sessions and ROI info
[all_sessions selected_sessions] = ioi_get_sessions(job);
[all_ROIs selected_ROIs] = ioi_get_ROIs(job);

%Big loop over subjects
for SubjIdx=1:length(job.IOImat)
    try
        tic
        clear IOI ROI SPM
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
                        % Load filtered downsampled signals
                        filtNdownData = load(IOI.fcIOS.filtNdown.fname);
                        fnameROIregress = fullfile(newDir,'ROIregress.mat');
                        % Loop over sessions
                        for s1=1:length(IOI.sess_res)
                            if all_sessions || sum(s1==selected_sessions)
                                sessionDir = [newDir filesep 'S' sprintf('%02d',s1)];
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
                                            dim = [1 1 1/IOI.fcIOS.filtNdown.downFreq];
                                            if job.wholeImage
                                                %% GLM on images here!
                                                % y = ioi_get_images(IOI,1:IOI.sess_res{s1}.n_frames,c1,s1,dir_ioimat,shrinkage_choice);
                                                
                                                % Constructing inputs required for
                                                % GLM analysis within the SPM
                                                % framework
                                                SPM.xY.VY = spm_vol(IOI.fcIOS.filtNdown.fnameWholeImage{s1, c1});
                                                y = spm_read_vols(SPM.xY.VY);
                                                % Preallocating output images
                                                yRegress = zeros(size(y));
                                                
                                                % All regressors are identified
                                                % here, lets take the mean global
                                                % signal just to test
                                                SPM.xX.name = cellstr(['Global Brain Signal']);
                                                SPM.xX.X = brainSignal';        % Regression is along first dimension. For one regressor it is a column vector.
                                                
                                                % A revoir
                                                SPM.xX.iG = [];
                                                SPM.xX.iH = [];
                                                SPM.xX.iC = 1:size(SPM.xX,2);   % Indices of regressors of interest
                                                SPM.xX.iB = [];                 % Indices of confound regressors
                                                
                                                SPM.xVi.Vi = {speye(size(SPM.xX.X,1))}; % Time correlation
                                                
                                                % Directory to save SPM and img/hdr
                                                % files
                                                SPM.swd = colorDir;
                                                if ~exist(SPM.swd,'dir'),mkdir(SPM.swd); end
                                                
                                                % GLM is performed here
                                                SPM = spm_spm(SPM);
                                                
                                                if job.regressBrainSignal == 1,
                                                    % Subtract global brain signal
                                                    % from every pixel time course
                                                    betaVol = spm_vol(fullfile(SPM.swd,SPM.Vbeta.fname));
                                                    beta = spm_read_vols(betaVol);
                                                    brainSignalRep(1,1,1,:) = brainSignal;
                                                    yRegress = y - repmat(beta,[1 1 1 size(y,4)]) .* repmat(brainSignalRep,[size(beta,1) size(beta,2) 1 1]);
                                                    
                                                    filtNdownfnameRegress = fullfile(sessionDir,[IOI.subj_name '_OD_' IOI.color.eng(c1) '_regress_' sprintf('%05d',1) 'to' sprintf('%05d',IOI.sess_res{s1}.n_frames) '.nii']);
                                                    % Save NIFTI file
                                                    ioi_save_nifti(yRegress, filtNdownfnameRegress, dim);
                                                    % Brain signal regression succesful!
                                                    IOI.fcIOS.SPM(1).wholeImageRegressOK{s1, c1} = true;
                                                end
                                                % Update SPM matrix info
                                                IOI.fcIOS.SPM(1).fnameSPM{s1, c1} = SPM.swd;
                                                IOI.fcIOS.SPM(1).fname{s1, c1} = filtNdownfnameRegress;
                                            end % end on GLM on images
                                            
                                            % Loop over ROIs
                                            for r1=1:length(IOI.res.ROI)
                                                if all_ROIs || sum(r1==selected_ROIs)
                                                    ROIdir = [colorDir filesep 'ROI' sprintf('%02d',r1)];
                                                    if ~exist(ROIdir,'dir'),mkdir(ROIdir); end
                                                    % Initialize y tilde (ROIregress)
                                                    ROIregress{r1}{s1,c1} = [];
                                                    %% Do my GLM on ROI code here
                                                    y = filtNdownData.filtNdownROI{r1}{s1, c1};
                                                    
                                                    % Display plots on SPM graphics window
                                                    spm_figure('GetWin', 'Graphics');
                                                    spm_figure('Clear', 'Graphics');
                                                    subplot(311); plot(y);
                                                    title(sprintf('Seed %d time-course, S%d, C%d',r1,s1,c1));
                                                    subplot(312); plot(brainSignal);
                                                    title(sprintf('Mean global signal time-course, S%d, C%d',s1,c1));
                                                    
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
                                                    SPM.xY.VY = spm_vol(fnameNIFTI);
                                                    
                                                    % All regressors are
                                                    % identified here, lets take
                                                    % the mean global signal
                                                    % just to test
                                                    SPM.xX.name = cellstr(['Global Brain Signal']);
                                                    SPM.xX.X = brainSignal';        % Regression is along first dimension. For one regressor it is a column vector.
                                                    
                                                    % A revoir
                                                    SPM.xX.iG = [];
                                                    SPM.xX.iH = [];
                                                    SPM.xX.iC = 1:size(SPM.xX,2);   % Indices of regressors of interest
                                                    SPM.xX.iB = [];                 % Indices of confound regressors
                                                    
                                                    SPM.xVi.Vi = {speye(size(SPM.xX.X,1))}; % Time correlation
                                                    
                                                    % Directory to save SPM and
                                                    % img/hdr files
                                                    SPM.swd = ROIdir;
                                                    if ~exist(SPM.swd,'dir'),mkdir(SPM.swd); end
                                                    % GLM is performed here
                                                    SPM = spm_spm(SPM);
                                                    
                                                    if job.regressBrainSignal == 1,
                                                        % Subtract global brain
                                                        % signal from ROI time
                                                        % courses
                                                        betaVol = spm_vol(fullfile(SPM.swd,SPM.Vbeta.fname));
                                                        beta = spm_read_vols(betaVol);
                                                        ROIregress{r1}{s1, c1} = y - beta * brainSignal;
                                                        
                                                        % Brain signal regression succesful!
                                                        IOI.fcIOS.SPM(1).ROIregressOK{r1}{s1, c1} = true;
                                                        
                                                        % Identify in IOI the file name of the time series
                                                        IOI.fcIOS.SPM(1).fnameROIregress = fnameROIregress;
                                                        colorNames = fieldnames(IOI.color);
                                                        fprintf('Global brain signal regressed from ROI %d (%s) Session %d Color %d (%s) done!\n',r1,IOI.ROIname{r1},s1,c1,colorNames{1+c1})
                                                        subplot(313); plot(ROIregress{r1}{s1, c1});
                                                        title(sprintf('Global signal regressed from ROI time-course %d, S%d, C%d',r1,s1,c1));
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
                                                    fprintf('GLM for ROI %d Session %d Color %d done!\n',r1,s1,c1)
                                                end
                                            end
                                        end
                                    end % ROI loop
                                end % colors loop
                            end
                        end % sessions loop
                        % GLM regression succesful!
                        IOI.fcIOS.SPM(1).GLMOK = true;
                        if job.regressBrainSignal == 1,
                            save(fnameROIregress,'ROIregress');
                        end
                        save(IOImat,'IOI');
                    end % GLM OK or redo job
                end % Filtering&Downsampling OK
            end % Time-series OK
        end % ROI OK
        toc
        disp(['Subject ' int2str(SubjIdx) ' (' IOI.subj_name ')' ' complete']);
        out.IOImat{SubjIdx} = IOImat;
    catch exception
        out.IOImat{SubjIdx} = IOImat;
        disp(exception.identifier)
        disp(exception.stack(1))
    end
end % Big loop over subjects
end

%spm get spm
