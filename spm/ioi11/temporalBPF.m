function Yfilt = temporalBPF(type, fs, cutoff, FilterOrder, Y)
% Butterworth type Band-pass filter (1-D)
% SYNTAX
% function Yfilt = ButterBPF(fs, cutoff, FilterOrder, Y)
% INPUTS
% type          String specifying the type of filter to use:
%               'butter'
%               'cheby1'
%               'cheby2'
%               'ellip'
% fs            Sampling frequency(in seconds)
% cutoff        A two-element vector fn (in Hz) that must be 0.0 < fn < fs/2
% FilterOrder   An integer (usually N=4)
% Y             The 1-D signal to be filtered
% OUTPUTS
% Yfilt         The filtered 1-D signal
%_______________________________________________________________________________
% Copyright (C) 2012 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________
try
    Wn = cutoff*2 / fs;                 % Normalised cut-off frequency
    
    switch type
        case 'butter'                   % Butterworth filter
            [z, p, k] = butter(FilterOrder, Wn, 'bandpass');
        case 'cheby1'                   % Chebyshev I filter
            % R dB of peak-to-peak ripple in the passband
            R = 0.5;
            [z, p, k] = cheby1(FilterOrder, R, Wn, 'bandpass');
        case 'cheby2'                   % Chebyshev II filter
            % stopband ripple R dB down from the peak passband value
            R = 80;
            [z, p, k] = cheby2(FilterOrder, R, Wn, 'bandpass');
        case 'ellip'                    % Elliptic filter
            % Rp dB of ripple in the passband, and a stopband Rs dB down from
            % the peak value in the passband
            Rp = .1; Rs = 80;
            % [b,a] = ellip(FilterOrder, Rp, Rs, Wn, 'bandpass');
            % Due to edge artifacts, one should use the [z,p,k] syntax to design
            % IIR filters.
            [z, p, k] = ellip(FilterOrder, Rp, Rs, Wn, 'bandpass');
        case 'yulewalk'                 % Recursive digital filter design
            f = [0 Wn 1];
            m = [0 1 1 0];
            [b,a] = yulewalk(FilterOrder, f, m);
        otherwise
            fprintf('Filter %s not available \n',type);
            Yfilt = Y;
            return
    end
    
    % Yfilt = filtfilt(b, a, Y);
    % Zero-phase forward and reverse digital filtering filtfilt disabled because
    % the method is not implemented in the dfilt object, only filter is.
    
    % Convert zero-pole-gain filter parameters to second-order sections form
    [sos, g] = zp2sos(z,p,k);
    % Create a dicrete-time filter object
    %       hFilt = dfilt.df2sos(sos,g);
    %     % Perform the filtering
    %     Yfilt = filter(hFilt, Y);
    
    Yfilt = ioi_filtfilt(sos, g, Y);
    
catch exception
    Yfilt = Y;
    disp(exception.identifier)
    disp(exception.stack(1))
end

% EOF
