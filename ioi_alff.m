function alff = ioi_alff (y, varargin)
% ioi_alff Computes Amplitude of Low Frequency Fluctuations
% SYNTAX
% alff = ioi_alff (y, varargin)
% INPUTS
% y         Time series vector sampled at fs
% yGlobal   [OPTIONAL] Time series of all pixels marked as brain, if given,
%           then ALFF is normalized by this value
% OUTPUT 
% alff  Amplitude of Low Frequency Fluctuations between 0.008 & 0.1 Hz
fs = 5;         % Hardcoded at 5 Hz
[X, freq] = ioi_positiveFFT(y, fs);
X = abs(X);
idxLow = find(freq<=0.008,1,'last');
idxHigh = find(freq>=0.1,1,'first');
switch nargin
    case 1
        alff = mean(sqrt(X(idxLow:idxHigh)));
    case 2
        yGlobal = varargin{1};
        alff = mean(sqrt(X(idxLow:idxHigh))) ./ yGlobal;
    otherwise
end
end

% EOF