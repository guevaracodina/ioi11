function out = ioi_HDM_run(job) 
%Based on SPM HDM 

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
%select onsets
if isfield(job.stim_choice,'electro_stims')
    electro_stims = 1;    
else
    %default stims
    electro_stims = 0;
end

%Big loop over subjects
for SubjIdx=1:length(job.IOImat)
    try
        tic
        clear IOI ROI SPM
        %Load IOI.mat information
        IOImat = job.IOImat{SubjIdx};
        load(IOImat);
        %Loop over sessions and ROIs

        if ~isfield(IOI.res,'ROIOK')
            disp(['No ROI available for subject ' int2str(SubjIdx) ' ... skipping series extraction']);
        else
            if ~isfield(IOI.res,'seriesOK')
                disp(['Extracted series not available for subject ' int2str(SubjIdx) ' ... skipping GLM']);
            else
                if ~isfield(IOI.res,'HDMOK') || job.force_redo
                    %check that onsets are available
                    if electro_stims
                        
                    
                    [dir_ioimat dummy] = fileparts(job.IOImat{SubjIdx});
                    if isfield(job.IOImatCopyChoice,'IOImatCopy')
                        newDir = job.IOImatCopyChoice.IOImatCopy.NewIOIdir;
                        newDir = fullfile(dir_ioimat,newDir);
                        if ~exist(newDir,'dir'),mkdir(newDir); end
                        IOImat = fullfile(newDir,'IOI.mat');
                        cdir = newDir; %current directory
                    else
                        cdir = dir_ioimat;
                    end
                    %loop over sessions
                    for s1=1:length(IOI.sess_res)
                        if all_sessions || sum(s1==selected_sessions)
                            %Choose onsets: stimulations or electrophysiology

                            %Electrophysiology, for each subject and session
                            if electro_stims
                                [pkh ons] = ioi_get_onsets(IOI,s1,E,cdir); %pk in seconds; pkh in arbitrary units
                                dur = 1;
                                name = '';
                                %                                 separate_slow_fast = 1;
                                %                                 if separate_slow_fast
                                %                                     name{1}
                                %                                 end
                            else
                                ons = IOI.sess_res{s1}.onsets{1}; %already in seconds *IOI.dev.TR;
                                dur = IOI.sess_res{s1}.durations{1}; %*IOI.dev.TR;
                                name = IOI.sess_res{s1}.names{1};
                            end


out=[];
M = []; %may need to fill this with subject-dependent information, etc.
M.Model_Choice = job.Model_Choice;
M.curveToFit = 3; %For now: 'HbR','HbT','Flow' 
%M.GroupPriors
%choose model
M=ioi_HDM_set_model(M); 
%choose priors
M=ioi_HDM_set_priors(M);


                  
                    IOI.res.GLMOK = 1;
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
end
