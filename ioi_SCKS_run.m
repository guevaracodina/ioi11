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
        clear IOI ROI
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
                    if ~isempty(job.ROImat)
                        load(job.ROImat{SubjIdx});
                    else
                        load(IOI.ROI.ROIfname);
                    end
                    %Loop to find HbO, HbR and Flow positions in ROI
                    for c1=1:length(IOI.color.eng)
                        if IOI.color.eng(c1)==IOI.color.HbO, cHbO = c1; end
                        if IOI.color.eng(c1)==IOI.color.HbR, cHbR = c1; end
                        if IOI.color.eng(c1)==IOI.color.flow, cFlow = c1; end
                    end
                        
                    %loop over sessions
                    for s1=1:length(IOI.sess_res)
                        if all_sessions || sum(s1==selected_sessions)
                            %loop over ROIs
                            for r1=1:length(IOI.res.ROI)
                                if all_ROIs || sum(r1==selected_ROIs)
                                    %SCKS on ROI
                                    SCKS = [];
                                    SCKS.SCKSparams = job.SCKSparams;
                                    %Data repetition rate (inverse of
                                    %sampling frequency)
                                    SCKS.TR = IOI.dev.TR; 
                                    %data specification - which modalities to include: 
                                    SCKS.PS.xY.cHbO = cHbO;
                                    SCKS.PS.xY.cHbR = cHbR;
                                    SCKS.PS.xY.cFlow = cFlow;
                                    SCKS.PS.xY.includeHbR = includeHbR;
                                    SCKS.PS.xY.includeHbT = includeHbT;
                                    SCKS.PS.xY.includeFlow = includeFlow;  
                                    SCKS.HPF = HPF; %High pass filter on data                                   
                                    SCKS = ioi_SCKS_get_data(ROI,SCKS,r1,s1);
                                    SCKS.PS.PhysioModel_Choice = job.PhysioModel_Choice;
                                    SCKS = ioi_SCKS_set_model(SCKS);
                                    SCKS = ioi_SCKS_set_priors(SCKS);
                                    try 
                                        U = IOI.Sess(s1).U;
                                    catch
                                        U = [];
                                    end
                                    SCKS = ioi_SCKS_set_SCKS(SCKS,U);
                                    if job.SCKSparams.SCKSnoise
                                        SCKS = ioi_SCKS(SCKS,0); 
                                    else
                                        SCKS = ioi_SCKS2(SCKS,0);  
                                    end
                                end
                            end
                            
                            disp(['SCKS for session ' int2str(s1) ' completed']);                            
                        end
                    end
                    
                    IOI.res.SCKSOK = 1;
                    if isfield(job.IOImatCopyChoice,'IOImatCopy')
                        newDir = job.IOImatCopyChoice.IOImatCopy.NewIOIdir;
                        newDir = fullfile(dir_ioimat,newDir);
                        if ~exist(newDir,'dir'),mkdir(newDir); end
                        IOImat = fullfile(newDir,'IOI.mat');
                    end
                    [dir1 dummy] = fileparts(IOImat);
                    ROIfname = fullfile(dir1,'ROI.mat');
                    save(ROIfname,'ROI');
                    IOI.ROI.ROIfname = ROIfname;
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
    end
end