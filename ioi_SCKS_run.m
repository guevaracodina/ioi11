function out = ioi_SCKS_run(job)
[O HPF LPF] = get_common_O_HDM_SCKS(job);
O.only_display = job.only_display;
show_figures = job.generate_figures;
%Big loop over subjects
for SubjIdx=1:length(job.IOImat)
    try
        tic
        clear IOI ROI SCKS
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
            disp(['No ROI available for subject ' int2str(SubjIdx) ' ... skipping SCKS']);
        else
            if ~isfield(IOI.res,'seriesOK')
                disp(['No extracted series on ROIs available for subject ' int2str(SubjIdx) ' ... skipping SCKS']);
            else
                if ~isfield(IOI.res,'SCKSOK') || job.force_redo
                    try
                        load(IOI.ROI.ROIfname);
                    catch
                        load(fullfile(dir_ioimat,'ROI.mat'));
                    end
                    [dir1 dummy] = fileparts(IOImat);
                    SCKSfname = fullfile(dir1,'SCKS.mat');
                    O = ioi_get_PS(IOI,O);
                    if O.all_sessions
                        O.selected_sessions = 1:length(IOI.sess_res);
                    end
                    if O.all_ROIs
                        O.selected_ROIs = 1:length(IOI.res.ROI);
                    end
                    clear SCKS
                    %loop over sessions
                    for s1=1:length(IOI.sess_res)
                        if O.all_sessions || sum(s1==O.selected_sessions)
                            %loop over ROIs
                            for r1=1:length(IOI.res.ROI)
                                if O.all_ROIs || sum(r1==O.selected_ROIs)
                                    if ~O.only_display
                                        %SCKS on ROI
                                        SCKS0 = [];
                                        SCKS0.O = O;
                                        SCKS0.SCKSparams = job.SCKSparams;
                                        %Data repetition rate (inverse of sampling frequency)
                                        SCKS0.dt = IOI.dev.TR;
                                        SCKS0.HPF = HPF; %High pass filter on data
                                        SCKS0.LPF = LPF;
                                        SCKS0 = ioi_get_data(ROI,SCKS0,r1,s1);
                                        SCKS0 = ioi_set_physiomodel(SCKS0);
                                        SCKS0 = ioi_set_priors(SCKS0);
                                        try %onsets may not be available, are optional
                                            U = IOI.Sess(s1).U;
                                        catch
                                            U = [];
                                        end
                                        SCKS0 = ioi_SCKS_set_SCKS(SCKS0,U);
                                        if job.SCKSparams.SCKSnoise
                                            SCKS0 = ioi_SCKS(SCKS0,0); %noise and parameters promoted to states
                                        else
                                            SCKS0 = ioi_SCKS2(SCKS0,0); %reduced version, noise and parameters are time independent
                                        end
                                        %Store and save results obtained so far
                                        SCKS{r1,s1} = SCKS0;
                                        save(SCKSfname,'SCKS');
                                    else
                                        try 
                                            if ~exist('SCKS','var')
                                                load(SCKSfname);
                                            end
                                            SCKS0 = SCKS{r1,s1};
                                        catch
                                            disp(['SCKS.mat not found for Session ' int2str(s1) ' and ROI ' int2str(r1)]);
                                        end
                                    end 
                                    %show results
                                    if show_figures
                                    end
                                end
                            end
                            disp(['SCKS for session ' int2str(s1) ' completed']);
                        end
                    end
                    IOI.res.SCKSOK = 1;
                    IOI.SCKS.fname = SCKSfname;
                    save(IOImat,'IOI');
                end
            end
        end
        toc
        disp(['Subject ' int2str(SubjIdx) ' complete']);
        out.IOImat{SubjIdx} = IOImat;
    catch exception
        disp(exception.identifier)
        disp(exception.stack(1))
        out.IOImat{SubjIdx} = job.IOImat{SubjIdx};
    end
end