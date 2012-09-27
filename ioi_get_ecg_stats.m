function ecgStats = ioi_get_ecg_stats(job, IOI)
% Gets ECG stats per subject, such as bpm (mean and standard deviation)
% SYNTAX
% oct_fix_struct_filename(subjectFolder)
% INPUT 
% job
% IOI
% OUTPUT 
% ecgStats  Structure with information about cardiac rythm of the mouse,
%           specifically BPM (beats per minute) mean and standard deviation.
% 
%_______________________________________________________________________________
% Copyright (C) 2012 LIOM Laboratoire d'Imagerie Optique et Mol�culaire
%                    �cole Polytechnique de Montr�al
%_______________________________________________________________________________

% Get sessions info
% [all_sessions selected_sessions] = ioi_get_sessions(job);

all_sessions = 1;
selected_sessions = [];

% Add path to Find_ECG_peaks.m file 
% addpath(genpath('D:\spm8\toolbox\oct12'));

% Loop over sessions
for s1=1:length(IOI.sess_res)
    if all_sessions || sum(s1==selected_sessions)
        load(IOI.res.elinfo{s1})
        % ChanNames{1}{1}
        % ConvertedData.Data.MeasuredData(1,3).Property(1,3).Name
        dt = ConvertedData.Data.MeasuredData(1,3).Property(1,3).Value;
        % nSamples = ConvertedData.Data.MeasuredData(1,3).Total_Samples;
        ecgVector = ConvertedData.Data.MeasuredData(1,3).Data;
        % tVector = 0:dt:dt*(nSamples-1);
        bpm(s1) = ioi_find_ECG_peaks(ecgVector, dt);
    end
end % sessions loop

% Cardiac rhythm statistics for subject
ecgStats.bpm_mean = nanmean(bpm);
ecgStats.bpm_std = nanstd(bpm);

function bpm = ioi_find_ECG_peaks(ecgVector, dt)
% Wrapper function to use OCT tool Find_ECG_peaks
acqui_info.ecg_signal = ecgVector;

% ECG sampling period in seconds
acqui_info.line_period_us = dt*1e6;

% Number of A-lines per ramp
acqui_info.ramp_length = 945;

% Do not display plot
display_plots = false;

% Call oct12 function
[~, bpm] = Find_ECG_peaks(acqui_info,display_plots);

% EOF
