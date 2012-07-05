function Y = ButterBPF(fs, cutoff, FilterOrder, Y)
% Butterworth type Band-pass filter (1-D)
% SYNTAX
% function Y = ButterBPF(fs, cutoff, FilterOrder, Y)
% INPUTS
% fs            Sampling frequency(in seconds)
% cutoff        A two-element vector fn (in Hz) that must be 0.0 < fn < fs/2
% FilterOrder   An integer (usually N=4)
% Y             The 1-D signal to be filtered
% OUTPUTS
% Y             The filtered 1-D signal
%_______________________________________________________________________________
% Copyright (C) 2012 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________

Wn = cutoff*2 / fs;                 % Normalised cut-off frequency
[fb,fa] = butter(FilterOrder, Wn);	% Butterworth filter
Y = filtfilt(fb, fa, Y); 
end
