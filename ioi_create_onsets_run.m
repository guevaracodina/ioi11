function out = ioi_create_onsets_run(job)
%select onsets, default is stimulation based
[E stim_choice] = ioi_get_E_for_electrophysiology(job);

if isfield(job.stim_choice,'manual_stims')
    stim_choice = 2;
    if ~isempty(job.stim_choice.manual_stims.onset_type_info)
        oti = job.stim_choice.manual_stims.onset_type_info;
    else
        oti = [];
    end
end

%select a subset of sessions
[all_sessions selected_sessions] = ioi_get_sessions(job);

elDir = job.elDir;
%Big loop over subjects
for SubjIdx=1:length(job.IOImat)
    try
        tic
        clear IOI
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
        
        if ~isfield(IOI.res,'OnsetsOK') || job.force_redo
            %loop over sessions
            for s1=1:length(IOI.sess_res)              
                if all_sessions || sum(s1==selected_sessions)
                    if all_sessions
                        selected_sessions = 1:length(IOI.sess_res);
                    end
                    switch stim_choice
                        case 0
                            %Just check that default stimulations are present
                            if ~isfield(IOI.sess_res{s1},'onsets')
                                disp(['No onsets for subject ' int2str(SubjIdx) ', session ' int2str(s1)]);
                            end
                        case 1
                            %Electrophysiology, for each subject and session
                            %try to specify a valid folder for electrophysiology
                            elDir0 = [];
                            if ~isempty(elDir)
                                if length(job.IOImat) == length(elDir)
                                    elDir0 = elDir{SubjIdx};
                                end
                            end
                            [pkh ons dur] = ioi_get_onsets_from_electrophysiology(IOI,s1,E,newDir,elDir0); %pk in seconds; pkh in arbitrary units
                            ot = 1;
                            IOI.sess_res{s1}.E = E; %Electrophysiology structure used for detection
                            IOI.sess_res{s1}.names{ot} = E.electrophysiology_onset_name;
                            IOI.sess_res{s1}.onsets{ot} = ons;
                            if E.spkOn
                                IOI.sess_res{s1}.durations{ot} = IOI.dev.TR;
                            else
                                IOI.sess_res{s1}.durations{ot} = dur;
                            end
                            IOI.sess_res{s1}.parameters{ot} = pkh;
                        case 2
                            if isempty(oti)
                                spm_input(['Subject ' int2str(SubjIdx) ', Session ' int2str(s1)],'-1','d');
                                p = 1;
                                ot = 0;
                                while p
                                    p = spm_input('Add an onset type?',2,'y/n');
                                    if p == 'y', p = 1; else p = 0; end
                                    if p
                                        ot = ot + 1;
                                        IOI.sess_res{s1}.names{ot} = spm_input('Enter onset name',3,'s');
                                        IOI.sess_res{s1}.onsets{ot} = spm_input('Enter onsets in seconds',4,'e','',NaN);
                                        IOI.sess_res{s1}.durations{ot} = spm_input('Enter durations in seconds',5,'e','',NaN);
                                        IOI.sess_res{s1}.parameters{ot} = spm_input('Enter amplitude in seconds',6,'e','',NaN);
                                        try close(h2); end
                                    end
                                end
                            else
                                %find whether specified for each subject or to apply to all subjects
                                if length(oti) == length(job.IOImat) || length(oti) == 1
                                    %OK
                                    if length(oti) == 1
                                        subj = 1;
                                    else
                                        subj = SubjIdx;
                                    end
                                    if length(oti{subj}) == length(selected_sessions) || length(oti{subj}) == 1
                                        if length(oti{subj}) == 1
                                            sess = 1;
                                        else
                                            sess = find(s1==selected_sessions);
                                        end
                                        %loop over onset types
                                        for ot=1:length(oti{subj}{sess})
                                            IOI.sess_res{s1}.names{ot} = oti{subj}{sess}(ot).onset_name;
                                            IOI.sess_res{s1}.onsets{ot} = oti{subj}{sess}(ot).onset_times;
                                            IOI.sess_res{s1}.durations{ot} = oti{subj}{sess}(ot).onset_durations;
                                            IOI.sess_res{s1}.parameters{ot} = oti{subj}{sess}(ot).onset_amplitude;
                                        end
                                    else
                                        disp('Problem with number of sessions for which onsets were specified');
                                    end
                                else
                                    disp('Problem with number of subjects for which onsets were specified');
                                end
                            end
                    end
                end
            end
            IOI.res.OnsetsOK = 1;
            save(IOImat,'IOI');
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