function out = ioi_filtdown_run(job)
% Band-pass filters (usually [0.009 - 0.8]Hz) and downsamples (usually to 1 Hz)
% a time series of an ROI (seed) and the whole brain mask.
%_______________________________________________________________________________
% Copyright (C) 2012 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________

for SubjIdx=1:length(job.IOImat)
    try
        tic
        clear IOI
        %Load IOI.mat information
        IOImat = job.IOImat{SubjIdx};
        load(IOImat);
        
        % IOI copy/overwrite
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
        
        if ~isfield(IOI.fcIOS.filtNdown,'filtNdownOK') || job.force_redo
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
            
            % Filter order
            filterOrder = 4;
            
            % Retrieve data
            ROIdata = load(IOI.ROI.ROIfname);
            ROIdata = ROIdata.ROI;
            brainMaskData = load(IOI.fcIOS.mask.fnameSeries);
            brainMaskData = brainMaskData.brainMaskSeries;
            
            [all_sessions selected_sessions] = ioi_get_sessions(job);
            
            IC = job.IC;
            % Loop over sessions
            for s1=1:length(IOI.sess_res)
                if all_sessions || sum(s1==selected_sessions)
                    % Loop over available colors
                    for c1=1:length(IOI.sess_res{s1}.fname)
                        doColor = ioi_doColor(IOI,c1,IC);
                        if doColor
                            colorOK = 1;
                            if ~(IOI.color.eng(c1)==IOI.color.laser)
                                %skip laser - only extract for flow
                                [all_ROIs selected_ROIs] = ioi_get_ROIs(job);
                                msg_ColorNotOK = 1;
                                % Initialize output filtNdownROI
                                for r1 = length(ROIdata);
                                    if all_ROIs || sum(r1==selected_ROIs)
                                        filtNdownROI{r1}{s1,c1} = [];
                                        filtNdownBrain{1}{s1,c1} = [];
                                    end
                                end
                                % Loop over ROIs
                                for r1 = 1:length(ROIdata); % All the ROIs
                                    if all_ROIs || sum(r1==selected_ROIs)
                                        try
                                             % Retrieve time-series signal for
                                             % given ROI, session and color
                                             ROIsignal = ROIdata{r1}{s1, c1};
                                             % Retrieve time-course
                                             % signal for brain mask
                                             brainSignal = brainMaskData{1}{s1, c1};
                                             % Band-passs filtering
                                             ROIsignal = ButterLPF(fs,job.BPFfreq(2),filterOrder,ROIsignal);
                                             ROIsignal = ButterHPF(fs,job.BPFfreq(1),filterOrder,ROIsignal);
                                             brainSignal = ButterLPF(fs,job.BPFfreq(2),filterOrder,brainSignal);
                                             brainSignal = ButterHPF(fs,job.BPFfreq(1),filterOrder,brainSignal);
                                             % Downsampling
                                             ROIsignal = downsample(ROIsignal, samples2skip);
                                             brainSignal = downsample(brainSignal, samples2skip);
                                        catch
                                            if msg_ColorNotOK
                                                msg = ['Problem extracting for color ' int2str(c1) ', session ' int2str(s1) ...
                                                    ',region ' int2str(r1) ': size ROIcell= ' int2str(size(ROIcell{r1},1)) 'x' ...
                                                    int2str(size(ROIcell{r1},2)) ', but size image= ' int2str(size(tmp_d,1)) 'x' ...
                                                    int2str(size(tmp_d,2))];
                                                IOI = disp_msg(IOI,msg);
                                                msg_ColorNotOK = 0;
                                            end
                                            if colorOK
                                                try
                                                    % Retrieve time-series signal for
                                                    % given ROI, session and color
                                                    ROIsignal = ROIdata{r1}{s1, c1};
                                                    % Retrieve time-course
                                                    % signal for brain mask
                                                    brainSignal = brainMaskData{1}{s1, c1};
                                                    % Band-passs filtering
                                                    ROIsignal = ButterLPF(fs,job.BPFfreq(2),filterOrder,ROIsignal);
                                                    ROIsignal = ButterHPF(fs,job.BPFfreq(1),filterOrder,ROIsignal);
                                                    brainSignal = ButterLPF(fs,job.BPFfreq(2),filterOrder,brainSignal);
                                                    brainSignal = ButterHPF(fs,job.BPFfreq(1),filterOrder,brainSignal);
                                                    % Downsampling
                                                    ROIsignal = downsample(ROIsignal, samples2skip);
                                                    brainSignal = downsample(brainSignal, samples2skip);
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
                                    disp(['Filtering and Downsampling for session ' int2str(s1) ' and color ' IOI.color.eng(c1) ' completed']);
                                end
                            end
                        end
                    end % Colors loop
                end
            end % Sessions loop
            % ------------------------------------------------------------------
            
            % Filter and Downsampling succesful!
            IOI.fcIOS.filtNdown(1).filtNdownOK = true;
            
            % IOI copy/overwrite
            if isfield(job.IOImatCopyChoice,'IOImatCopy')
                newDir = job.IOImatCopyChoice.IOImatCopy.NewIOIdir;
                newDir = fullfile(dir_ioimat,newDir);
                if ~exist(newDir,'dir'),mkdir(newDir); end
                IOImat = fullfile(newDir,'IOI.mat');
            end
            [dir1 dummy] = fileparts(IOImat);
            filtNdownfname = fullfile(dir1,'filtNdown.mat');
            save(filtNdownfname,'filtNdownROI','filtNdownBrain');
            IOI.fcIOS.filtNdown.fname = filtNdownfname;
            save(IOImat,'IOI');
            toc
        end
        out.IOImat{SubjIdx} = IOImat;
        disp(['Subject ' int2str(SubjIdx) ' (' IOI.subj_name ')' ' complete']);
    catch exception
        disp(exception.identifier)
        disp(exception.stack(1))
        out.IOImat{SubjIdx} = job.IOImat{SubjIdx};
    end % End of try
end % End of main for
end % End of function
