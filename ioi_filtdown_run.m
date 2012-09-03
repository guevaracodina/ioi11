function out = ioi_filtdown_run(job)
% Band-pass filters (usually [0.009 - 0.8]Hz) and downsamples (usually to 1 Hz)
% the time series of an ROI (seed) and the whole brain mask.
%_______________________________________________________________________________
% Copyright (C) 2012 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________

for SubjIdx=1:length(job.IOImat)
    try
        tic
        %Load IOI.mat information
        [IOI IOImat dir_ioimat]= ioi_get_IOI(job,SubjIdx);
        
        if ~isfield(IOI.fcIOS.filtNdown,'filtNdownOK') || job.force_redo
            % Get shrinkage options
            shrinkage_choice = IOI.res.shrinkageOn;
            % Filter & Downsample time series
            % ------------------------------------------------------------------
            % Original Sampling Frequency (5 Hz per color, data is sampled at 20
            % Hz for 4 colors RGYL)
            fs = 1/IOI.dev.TR;
            
            % Desired downsampling frequency
            fprintf('Desired downsampling frequency: %0.1f Hz \n',job.downFreq);
            
            % Real downsampling frequency
            samples2skip = round(fs/job.downFreq);
            fprintf('Real downsampling frequency: %0.1f Hz \n',fs/samples2skip);
            IOI.fcIOS.filtNdown(1).fs = fs/samples2skip;
            
            % Filter order
            filterOrder = job.bpf.bpf_On.bpf_order;
            
            % Band-pass cut-off frequencies
            BPFfreq = job.bpf.bpf_On.bpf_freq;
            
            % Filter type
            fType = job.bpf.bpf_On.bpf_type;
            
            % Retrieve data
            if IOI.res.seriesOK
                ROIdata = load(IOI.ROI.ROIfname);
                ROIdata = ROIdata.ROI;
            end
            if job.wholeImage && IOI.fcIOS.mask.seriesOK
                brainMaskData = load(IOI.fcIOS.mask.fnameSeries);
                brainMaskData = brainMaskData.brainMaskSeries;
            end
            
            [all_sessions selected_sessions] = ioi_get_sessions(job);
            
            % Include colors
            IC = job.IC;
            % Loop over sessions
            for s1=1:length(IOI.sess_res)
                if all_sessions || sum(s1==selected_sessions)
                    colorNames = fieldnames(IOI.color);
                    % Loop over available colors
                    for c1=1:length(IOI.sess_res{s1}.fname)
                        doColor = ioi_doColor(IOI,c1,IC);
                        if doColor
                            colorOK = 1;
                            if ~(IOI.color.eng(c1)==IOI.color.laser)
                                if job.wholeImage
                                    %% Filtering & Downsampling whole images (y)
                                    y = ioi_get_images(IOI,1:IOI.sess_res{s1}.n_frames,c1,s1,dir_ioimat,shrinkage_choice);
                                    % Preallocating output images
                                    filtY = zeros([size(y,1) size(y,2) 1 size(y,3)]);
                                    % Computing lenght of time vector
                                    nT = ceil(size(y,3) / samples2skip);
                                    filtNdownY = zeros([size(y,1) size(y,2) 1 nT]);
                                    % Read brain mask file
                                    vol = spm_vol(IOI.fcIOS.mask.fname);
                                    brainMask = logical(spm_read_vols(vol));
                                    % Test if there was shrinkage
                                    if size(brainMask,1)~= size(y,1)|| size(brainMask,2)~= size(y,2)
                                        brainMask = ioi_MYimresize(brainMask, [size(y,1) size(y,2)]);
                                    end
                                    % Color names
                                    colorNames = fieldnames(IOI.color);
                                    % Initialize progress bar
                                    spm_progress_bar('Init', size(filtY,1), sprintf('Filtering & Downsampling session %d, color %d (%s)\n',s1,c1,colorNames{1+c1}), 'Pixels along X');
                                    for iX = 1:size(filtY,1)
                                        spm_progress_bar('Set', iX);
                                        for iY = 1:size(filtY,2)
                                            if brainMask(iX,iY) == 1
                                                % Only non-masked pixels are
                                                % band-passs filtered
                                                filtY(iX,iY,1,:) = temporalBPF(fType,fs,BPFfreq,filterOrder,squeeze(y(iX,iY,:)));
                                                % Downsampling
                                                filtNdownY(iX,iY,1,:) = downsample(squeeze(filtY(iX,iY,1,:)), samples2skip);
                                            end
                                        end
                                    end
                                    spm_progress_bar('Clear');
                                    % Saving images
                                    sessionDir = [dir_ioimat filesep 'S' sprintf('%02d',s1)];
                                    if ~exist(sessionDir,'dir'),mkdir(sessionDir); end
                                    filtNdownfnameWholeImage = fullfile(sessionDir,[IOI.subj_name '_OD_' IOI.color.eng(c1) '_filtNdown_' sprintf('%05d',1) 'to' sprintf('%05d',IOI.sess_res{s1}.n_frames) '.nii']);
                                    ioi_save_nifti(filtNdownY, filtNdownfnameWholeImage, [1 1 samples2skip/fs]);
                                    IOI.fcIOS.filtNdown.fnameWholeImage{s1, c1} = filtNdownfnameWholeImage;
                                    if IOI.fcIOS.mask.seriesOK
                                        % Retrieve time-course
                                        % signal for brain mask
                                        brainSignal = brainMaskData{1}{s1, c1};
                                        % Band-passs filtering
                                        brainSignal = temporalBPF(fType,fs,BPFfreq,filterOrder,brainSignal);
                                        % Downsampling
                                        brainSignal = downsample(brainSignal, samples2skip);
                                        % Update data cell
                                        filtNdownBrain{1}{s1,c1} = brainSignal;
                                    end
                                    fprintf('Filtering and downsampling whole images for session %d and color %d (%s) completed\n',s1,c1,colorNames{1+c1})
                                end % End of filtering & downsampling whole images
                                %skip laser - only extract for flow
                                [all_ROIs selected_ROIs] = ioi_get_ROIs(job);
                                nROI = 1:length(IOI.res.ROI); % All the ROIs
                                msg_ColorNotOK = 1;
                                % Initialize output filtNdownROI
                                for r1 = nROI;
                                    if all_ROIs || sum(r1==selected_ROIs)
                                        filtNdownROI{r1}{s1,c1} = [];
                                    end
                                end
                                % Loop over ROIs
                                for r1 = nROI; % All the ROIs
                                    if all_ROIs || sum(r1==selected_ROIs)
                                        try
                                            % Retrieve time-series signal for
                                            % given ROI, session and color
                                            ROIsignal = ROIdata{r1}{s1, c1};
                                            
                                            % Band-passs filtering
                                            ROIsignal = temporalBPF(fType,fs,BPFfreq,filterOrder,ROIsignal);
                                            
                                            % Downsampling
                                            ROIsignal = downsample(ROIsignal, samples2skip);
                                            
                                            if job.generate_figures
                                                % ---- Plotting results ----
                                                % Display plots on SPM graphics window
                                                h = spm_figure('GetWin', 'Graphics');
                                                spm_figure('Clear', 'Graphics');
                                                % Positive FFT
                                                [X, freq] = ioi_positiveFFT(ROIsignal, fs);
                                                % Time vector
                                                t = 0:1/fs:(1/fs)*(numel(ROIsignal)-1);
                                                
                                                subplot(221)
                                                plot(t, ROIsignal,'k-','LineWidth',2)
                                                title(sprintf('%s_R%02d(%s)_S%02d_C%d(%s)\n',IOI.subj_name,r1,IOI.ROIname{r1},s1,c1,colorNames{1+c1}),'interpreter', 'none','FontSize',14);
                                                xlabel('t [s]','FontSize',14)
                                                set(gca,'FontSize',12)
                                                
                                                subplot(222)
                                                set(gca,'FontSize',12)
                                                semilogx(freq, abs(X),'k-','LineWidth',2);
                                                title(sprintf('Unfiltered spectrum'),'interpreter', 'none','FontSize',14);
                                                xlabel('f [Hz]','FontSize',14)
                                                % --------------------------
                                                % ---- Plotting results ----
                                                % Positive FFT
                                                [X, freq] = ioi_positiveFFT(ROIsignal, fs);
                                                subplot(223)
                                                plot(t, ROIsignal,'k-','LineWidth',2)
                                                title(sprintf('Filtered time-course\n'),'interpreter', 'none','FontSize',14);
                                                xlabel('t [s]','FontSize',14)
                                                set(gca,'FontSize',12)
                                                hold on
                                                
                                                subplot(224)
                                                set(gca,'FontSize',12)
                                                semilogx(freq, abs(X),'k-','LineWidth',2);
                                                % title(sprintf('%s_R%02d filtered spectrum S%d C%d (%s)\n',IOI.subj_name,r1,IOI.ROIname{r1},s1,c1,colorNames{1+c1}),'interpreter', 'none','FontSize',14);
                                                title('Filtered spectrum','interpreter', 'none','FontSize',14);
                                                xlabel('f [Hz]','FontSize',14)
                                                hold on
                                                % --------------------------
                                                % ---- Plotting results ----
                                                % Downsampled Time vector
                                                t = 0:1/IOI.fcIOS.filtNdown(1).fs:(1/IOI.fcIOS.filtNdown(1).fs)*(numel(ROIsignal)-1);
                                                % Positive FFT
                                                [X, freq] = ioi_positiveFFT(ROIsignal, IOI.fcIOS.filtNdown(1).fs);
                                                
                                                subplot(223)
                                                plot(t, ROIsignal,'r--','LineWidth',2)
                                                legend({'Filt' 'Filt&Down'},'FontSize',14)
                                                
                                                subplot(224)
                                                semilogx(freq, abs(X),'r--','LineWidth',2);
                                                legend({'Filt' 'Filt&Down'},'FontSize',14)
                                            end
                                            
                                            [oldDir, oldName, oldExt] = fileparts(IOI.res.ROI{1,1}.fname);
                                            newName = [sprintf('%s_R%02d_S%02d_C%d',IOI.subj_name,r1,s1,c1) '_filtNdown'];
                                            
                                            if job.save_figures
                                                % Save as PNG
                                                print(h, '-dpng', fullfile(dir_ioimat,newName), '-r300');
                                                % --------------------------
                                            end
                                        catch
                                            if msg_ColorNotOK
                                                msg = ['Problem filtering/downsampling for color ' int2str(c1) ', session ' int2str(s1) ...
                                                    ',region ' int2str(r1) ': size ROIsignal= ' int2str(size(ROIsignal,1)) 'x' ...
                                                    int2str(size(brainSignal,2)) ', but brainSignal= ' int2str(size(brainSignal,1)) 'x' ...
                                                    int2str(size(brainSignal,2))];
                                                IOI = disp_msg(IOI,msg);
                                                msg_ColorNotOK = 0;
                                            end
                                            if colorOK
                                                try
                                                    % Retrieve time-series signal for
                                                    % given ROI, session and color
                                                    ROIsignal = ROIdata{r1}{s1, c1};
                                                    
                                                    % Band-passs filtering
                                                    ROIsignal = temporalBPF(fType,fs,BPFfreq,filterOrder,ROIsignal);
                                                        
                                                    % Downsampling
                                                    ROIsignal = downsample(ROIsignal, samples2skip);
                                                    
                                                    if job.generate_figures
                                                        % ---- Plotting results ----
                                                        % Display plots on SPM graphics window
                                                        h = spm_figure('GetWin', 'Graphics');
                                                        spm_figure('Clear', 'Graphics');
                                                        % Positive FFT
                                                        [X, freq] = ioi_positiveFFT(ROIsignal, fs);
                                                        % Time vector
                                                        t = 0:1/fs:(1/fs)*(numel(ROIsignal)-1);
                                                        
                                                        subplot(221)
                                                        plot(t, ROIsignal,'k-','LineWidth',2)
                                                        title(sprintf('%s_R%02d(%s)_S%02d_C%d(%s)\n',IOI.subj_name,r1,IOI.ROIname{r1},s1,c1,colorNames{1+c1}),'interpreter', 'none','FontSize',14);
                                                        xlabel('t [s]','FontSize',14)
                                                        set(gca,'FontSize',12)
                                                        
                                                        subplot(222)
                                                        set(gca,'FontSize',12)
                                                        semilogx(freq, abs(X),'k-','LineWidth',2);
                                                        title(sprintf('Unfiltered spectrum'),'interpreter', 'none','FontSize',14);
                                                        xlabel('f [Hz]','FontSize',14)
                                                        % --------------------------
                                                        % ---- Plotting results ----
                                                        % Positive FFT
                                                        [X, freq] = ioi_positiveFFT(ROIsignal, fs);
                                                        subplot(223)
                                                        plot(t, ROIsignal,'k-','LineWidth',2)
                                                        title(sprintf('Filtered time-course\n'),'interpreter', 'none','FontSize',14);
                                                        xlabel('t [s]','FontSize',14)
                                                        set(gca,'FontSize',12)
                                                        hold on
                                                        
                                                        subplot(224)
                                                        set(gca,'FontSize',12)
                                                        semilogx(freq, abs(X),'k-','LineWidth',2);
                                                        % title(sprintf('%s_R%02d filtered spectrum S%d C%d (%s)\n',IOI.subj_name,r1,IOI.ROIname{r1},s1,c1,colorNames{1+c1}),'interpreter', 'none','FontSize',14);
                                                        title('Filtered spectrum','interpreter', 'none','FontSize',14);
                                                        xlabel('f [Hz]','FontSize',14)
                                                        hold on
                                                        % --------------------------
                                                        % ---- Plotting results ----
                                                        % Downsampled Time vector
                                                        t = 0:1/IOI.fcIOS.filtNdown(1).fs:(1/IOI.fcIOS.filtNdown(1).fs)*(numel(ROIsignal)-1);
                                                        % Positive FFT
                                                        [X, freq] = ioi_positiveFFT(ROIsignal, IOI.fcIOS.filtNdown(1).fs);
                                                        
                                                        subplot(223)
                                                        plot(t, ROIsignal,'r--','LineWidth',2)
                                                        legend({'Filt' 'Filt&Down'},'FontSize',14)
                                                        
                                                        subplot(224)
                                                        semilogx(freq, abs(X),'r--','LineWidth',2);
                                                        legend({'Filt' 'Filt&Down'},'FontSize',14)
                                                    end
                                                    
                                                    [oldDir, oldName, oldExt] = fileparts(IOI.res.ROI{1,1}.fname);
                                                    newName = [sprintf('%s_R%02d_S%02d_C%d',IOI.subj_name,r1,s1,c1) '_filtNdown'];
                                                    
                                                    if job.save_figures
                                                        % Save as PNG
                                                        print(h, '-dpng', fullfile(dir_ioimat,newName), '-r300');
                                                        % --------------------------
                                                    end
                                                    
                                                catch
                                                    msg = ['Unable to extract color ' int2str(c1) ', session ' int2str(s1)];
                                                    IOI = disp_msg(IOI,msg);
                                                    colorOK = 0;
                                                end
                                            end
                                        end
                                        if colorOK
                                            filtNdownROI{r1}{s1,c1} = ROIsignal;
                                            filtNdownBrain{1}{s1,c1} = brainSignal;
                                        end
                                    end
                                end % ROI loop
                                if colorOK
                                    filtNdownROI{r1}{s1,c1} = ROIsignal;
                                    filtNdownBrain{1}{s1,c1} = brainSignal;
                                    fprintf('Filtering and downsampling ROIs/seeds for session %d and color %d (%s) completed\n',s1,c1,colorNames{1+c1})
                                end
                            end
                        end
                    end % Colors loop
                end
            end % Sessions loop
            % ------------------------------------------------------------------
            
            % Filter and Downsampling succesful!
            IOI.fcIOS.filtNdown(1).filtNdownOK = true;
            
            filtNdownfname = fullfile(dir_ioimat,'filtNdown.mat');
            save(filtNdownfname,'filtNdownROI','filtNdownBrain');
            IOI.fcIOS.filtNdown.fname = filtNdownfname;
            % Desired downsampling frequency, it could be different to real
            % downsampling frequency (filtNdown.fs)
            IOI.fcIOS.filtNdown.downFreq = job.downFreq;
            IOI.fcIOS.filtNdown.BPFfreq = job.bpf.bpf_On.bpf_freq;
            save(IOImat,'IOI');
        end
        out.IOImat{SubjIdx} = IOImat;
        disp(['Elapsed time: ' datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS')]);
        disp(['Subject ' int2str(SubjIdx) ' (' IOI.subj_name ')' ' complete']);
    catch exception
        out.IOImat{SubjIdx} = job.IOImat{SubjIdx};
        disp(exception.identifier)
        disp(exception.stack(1))
    end % End of try
end % End of main for
end % End of function

% EOF
