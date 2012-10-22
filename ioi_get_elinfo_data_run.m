function out = ioi_get_elinfo_data_run(job)
%_______________________________________________________________________________
% Copyright (C) 2012 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________

for SubjIdx=1:length(job.IOImat)
    try
        % Load IOI.mat information
        [IOI IOImat dir_ioimat]= ioi_get_IOI(job,SubjIdx);
        if ~isfield(IOI,'elinfo') || job.force_redo
            % Updating elinfo data
            IOI = ioi_get_elinfo(IOI);
            % Turning off warning of ECG computation
            warning('off', 'MATLAB:polyfit:RepeatedPointsOrRescale')
            % Get ECG stats
            ecgStats{SubjIdx, 1} = ioi_get_ecg_stats(job, IOI);
            fprintf('%.2f ± %.2f bpm for subject %s\n', ...
                ecgStats{SubjIdx}.bpm_mean, ecgStats{SubjIdx}.bpm_std, IOI.subj_name)
            % Get temperature stats
            tempStats{SubjIdx, 1} = ioi_get_temp_stats(job, IOI);
            % Alt+0177 -> ±, Alt+0176 -> °.
            fprintf('%.2f ± %.2f °C for subject %s\n', ...
                tempStats{SubjIdx}.tMouse_mean, tempStats{SubjIdx}.tMouse_std, IOI.subj_name)
            subjectID{SubjIdx, 1} = IOI.subj_name;
        else
            % Do nothing, elinfo filenames are already saved
        end
        out.IOImat{SubjIdx} = IOImat;
    catch exception
        out.IOImat{SubjIdx} = job.IOImat{SubjIdx};
        disp(exception.identifier)
        disp(exception.stack(1))
    end % End of try
end % End of subjects loop
save (fullfile(IOI.dir.dir_group_res,'ecg_temp_stats_group.mat'),...
    'ecgStats', 'tempStats', 'subjectID');
end % End of function

% EOF
