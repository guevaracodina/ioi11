function out = ioi_get_elinfo_data_run(job)
%_______________________________________________________________________________
% Copyright (C) 2012 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________

for SubjIdx=1:length(job.IOImat)
    try
        tic
        %Load IOI.mat information
        [IOI IOImat dir_ioimat]= ioi_get_IOI(job,SubjIdx);     
        if ~isfield(IOI.fcIOS,'update_elinfo') || job.force_redo
            IOI = ioi_get_elinfo(IOI);
            % Updating elinfo data succesful!
            IOI.res.elinfoOK = true;
            save(IOImat,'IOI');
            disp(IOI.res.elinfo)
        elseif ~isfield(IOI.fcIOS.update_elinfo,'update_elinfoOK') || ~IOI.fcIOS.update_elinfo.update_elinfoOK
            IOI = ioi_get_elinfo(IOI);
            % Updating elinfo data succesful!
            IOI.res.elinfoOK = true;
            save(IOImat,'IOI');
            disp(IOI.res.elinfo)
        else
            % Do nothing, elinfo filenames are already saved
            disp(['elinfo already available for subject ' int2str(SubjIdx) ' (' IOI.subj_name ')']);
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
