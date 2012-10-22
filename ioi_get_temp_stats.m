function tempStats = ioi_get_temp_stats(job, IOI)
% Gets temperature stats per subject (mean and standard deviation)
% SYNTAX
% ecgStats = ioi_get_temp_stats(job, IOI)
% INPUT 
% job
% IOI
% OUTPUT 
% tempStats     Structure with information about temperature of both the mouse
%               and the warming pad.
%_______________________________________________________________________________
% Copyright (C) 2012 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________

% Get sessions info
[all_sessions selected_sessions] = ioi_get_sessions(job);

% Loop over sessions
for s1=1:length(IOI.sess_res)
    if all_sessions || sum(s1==selected_sessions)
        load(IOI.res.elinfo{s1})
        % Pad: ChanNames{1}{10} & Rat: ChanNames{1}{11}
        dt = ConvertedData.Data.MeasuredData(1,9).Property(1,3).Value;
        tMouse(s1) = nanmean(ConvertedData.Data.MeasuredData(1,9).Data);
        tPad(s1) = nanmean(ConvertedData.Data.MeasuredData(1,10).Data);
    end
end % sessions loop

% Cardiac rhythm statistics for subject
tempStats.tMouse_mean   = nanmean(tMouse);
tempStats.tMouse_std    = nanstd(tMouse);
tempStats.tPad_mean     = nanmean(tPad);
tempStats.tPad_std      = nanstd(tPad);

% EOF
