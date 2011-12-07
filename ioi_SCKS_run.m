function out = ioi_SCKS_run(job)
%select a subset of sessions
if isfield(job.session_choice,'select_sessions')
    all_sessions = 0;
    selected_sessions = job.session_choice.select_sessions.selected_sessions;
else
    all_sessions = 1;
end
%select a subset of ROIs
if isfield(job.ROI_choice,'select_ROIs')
    all_ROIs = 0;
    selected_ROIs = job.ROI_choice.select_ROIs.selected_ROIs;
else
    all_ROIs = 1;
end
%HPF
if isfield(job.hpf_butter,'hpf_butter_On')
    HPF.hpf_butter_On = 1;
    HPF.hpf_butter_freq = job.hpf_butter.hpf_butter_On.hpf_butter_freq;
    HPF.hpf_butter_order = job.hpf_butter.hpf_butter_On.hpf_butter_order;
else
    HPF.hpf_butter_On = 0;
end
%modalities to include
includeHbR = job.includeHbR;
includeHbT = job.includeHbT;
includeFlow = job.includeFlow;
%Big loop over subjects
for SubjIdx=1:length(job.IOImat)
    try
        tic
        clear IOI ROI SCKS
        %Load IOI.mat information
        IOImat = job.IOImat{SubjIdx};
        load(IOImat);
        
        if ~isfield(IOI.res,'ROIOK')
            disp(['No ROI available for subject ' int2str(SubjIdx) ' ... skipping SCKS']);
        else
            if ~isfield(IOI.res,'seriesOK')
                disp(['No extracted series on ROIs available for subject ' int2str(SubjIdx) ' ... skipping SCKS']);
            else
                if ~isfield(IOI.res,'SCKSOK') || job.force_redo
                    [dir_ioimat dummy] = fileparts(job.IOImat{SubjIdx});
                    %load ROI.mat
%                     if ~isempty(job.ROImat)
%                         load(job.ROImat{SubjIdx});
%                     else
%                         load(IOI.ROI.ROIfname);
%                     end
                    try
                        load(IOI.ROI.ROIfname);
                    catch
                        load(fullfile(dir_ioimat,'ROI.mat'));
                    end
                    if isfield(job.IOImatCopyChoice,'IOImatCopy')
                        newDir = job.IOImatCopyChoice.IOImatCopy.NewIOIdir;
                        newDir = fullfile(dir_ioimat,newDir);
                        if ~exist(newDir,'dir'),mkdir(newDir); end
                        IOImat = fullfile(newDir,'IOI.mat');
                    end
                    [dir1 dummy] = fileparts(IOImat);
                    SCKSfname = fullfile(dir1,'SCKS.mat');
                    PS0 = ioi_get_PS(IOI,includeHbR,includeHbT,includeFlow,job.PhysioModel_Choice);
                                         
                    %loop over sessions
                    for s1=1:length(IOI.sess_res)
                        if all_sessions || sum(s1==selected_sessions)
                            %loop over ROIs
                            for r1=1:length(IOI.res.ROI)
                                if all_ROIs || sum(r1==selected_ROIs)
                                    %SCKS on ROI
                                    SCKS0 = [];
                                    SCKS0.SCKSparams = job.SCKSparams;
                                    %Data repetition rate (inverse of
                                    %sampling frequency)
                                    SCKS0.TR = IOI.dev.TR; 
                                    %data specification - which modalities to include: 
                                    SCKS0.PS = PS0; 
                                    SCKS0.HPF = HPF; %High pass filter on data                                   
                                    SCKS0 = ioi_get_data(ROI,SCKS0,r1,s1);
                                    SCKS0 = ioi_set_physiomodel(SCKS0);
                                    SCKS0 = ioi_set_priors(SCKS0);
                                    try %onsets may not be available, are optional 
                                        U = IOI.Sess(s1).U;
                                    catch
                                        U = [];
                                    end
                                    %set priors
                                    SCKS0.pE = SCKS0.PS.pE;
                                    SCKS0.pC = SCKS0.PS.pC;
                                    SCKS0 = ioi_SCKS_set_SCKS(SCKS0,U);
                                    if job.SCKSparams.SCKSnoise
                                        SCKS0 = ioi_SCKS(SCKS0,0); %noise and parameters promoted to states 
                                    else
                                        SCKS0 = ioi_SCKS2(SCKS0,0); %reduced version, noise and parameters are time independent  
                                    end
                                    %Store and save results obtained so far
                                    SCKS{r1,s1} = SCKS0;
                                    save(SCKSfname,'SCKS');                 
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