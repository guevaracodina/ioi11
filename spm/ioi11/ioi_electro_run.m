function out = ioi_electro_run(job)
try
    for SubjIdx=1:length(job.IOImat)
        tic
        clear IOI
        %Load IOI.mat information
        IOImat = job.IOImat{SubjIdx};
        load(IOImat);
        %select a subset of sessions
        if isfield(job.session_choice,'select_sessions')
            all_sessions = 0;
            selected_sessions = job.session_choice.select_sessions.selected_sessions;
        else
            all_sessions = 1;
        end
        
        if ~isfield(IOI.res,'electroOK') || job.force_redo
            [dir_ioimat dummy] = fileparts(job.IOImat{SubjIdx});
            %Loop over sessions
            if isfield(IOI,'sess_raw')
                for s1=1:length(IOI.sess_raw)
                    if all_sessions || sum(s1==selected_sessions)
                        if isfield(IOI.sess_raw{s1},'electroOK')
                            %extract electro data and save
                            IOI = extract_electro(IOI,s1,IOI.sess_raw{s1}.electro{1});
                        end
                    end
                    disp(['Electrophysiology extraction for session ' int2str(s1) ' complete']);
                end
            end
            
            IOI.res.electroOK = 1;
            if isfield(job.IOImatCopyChoice,'IOImatCopy')
                newDir = job.IOImatCopyChoice.IOImatCopy.NewIOIdir;
                newDir = fullfile(dir_ioimat,newDir);
                if ~exist(newDir,'dir'),mkdir(newDir); end
                IOImat = fullfile(newDir,'IOI.mat');
            end
            save(IOImat,'IOI');
        end
        out.IOImat{SubjIdx} = IOImat;
        toc
        disp(['Subject ' int2str(SubjIdx) ' complete']);
    end
catch exception
    disp(exception.identifier)
    disp(exception.stack(1))
end


