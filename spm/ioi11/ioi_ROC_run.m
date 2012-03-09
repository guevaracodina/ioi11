function out = ioi_ROC_run(job)
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
%Big loop over subjects
for SubjIdx=1:length(job.IOImat)
    try
        tic
        clear IOI ROI SCKS
        %Load IOI.mat information
        IOImat = job.IOImat{SubjIdx};
        load(IOImat);
        if ~isfield(IOI.res,'SCKSOK')
            disp(['No SCKS available for subject ' int2str(SubjIdx) ' ... skipping ROC']);
        else
            if ~isfield(IOI.res,'ROC_OK') || job.force_redo
                [dir_ioimat dummy] = fileparts(job.IOImat{SubjIdx});
                if isfield(job.IOImatCopyChoice,'IOImatCopy')
                    newDir = job.IOImatCopyChoice.IOImatCopy.NewIOIdir;
                    newDir = fullfile(dir_ioimat,newDir);
                    if ~exist(newDir,'dir'),mkdir(newDir); end
                    IOImat = fullfile(newDir,'IOI.mat');
                end
                [dir1 dummy] = fileparts(IOImat);
                ROCfname = fullfile(dir1,'ROC.mat');
                load(IOI.SCKS.fname);
                %loop over sessions
                for s1=1:length(IOI.sess_res)
                    if all_sessions || sum(s1==selected_sessions)
                        %loop over ROIs
                        for r1=1:length(IOI.res.ROI)
                            if all_ROIs || sum(r1==selected_ROIs)
                                if ~isempty(SCKS{r1,s1})
                                    Options = [];
                                    ROC0 = ioi_ROC(SCKS{r1,s1},Options);
                                    %Store and save results obtained so far
                                    ROC{r1,s1} = ROC0;
                                    save(ROCfname,'ROC');
                                end
                            end
                        end
                        disp(['ROC for session ' int2str(s1) ' completed']);
                    end
                end
                IOI.res.ROC_OK = 1;
                IOI.SCKS.ROCfname = ROCfname;
                save(IOImat,'IOI');
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