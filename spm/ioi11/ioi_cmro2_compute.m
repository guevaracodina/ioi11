function cmro2 = ioi_cmro2_compute(dataFlow, dataHbO, dataHbR, varargin)
% Computes cerebral metabolic rate of oxygen (CMRO2) from blood flow, HbO and HbR data.
% Reference:
% M. Jones, J. Berwick, D. Johnston, and J. Mayhew, “Concurrent optical imaging
% spectroscopy and laser-Doppler flowmetry: the relationship between blood flow,
% oxygenation, and volume in rodent barrel cortex,” Neuroimage, vol. 13, no. 6
% Pt 1, pp. 1002–1015, Jun. 2001.
%_______________________________________________________________________________
% Copyright (C) 2012 LIOM Laboratoire d'Imagerie Optique et Moleculaire
%                    Ecole Polytechnique de Montreal
%_______________________________________________________________________________

% ------------------------------------------------------------------------------
% Optional inputs handling
% ------------------------------------------------------------------------------
% only want 1 optional input at most
numvarargs              = length(varargin);
if numvarargs > 2
    error('ioi_cmro2_compute:TooManyInputs', ...
        'Requires at most 2 optional inputs');
end
% set defaults for optional inputs
optargs                 = {ioi_get_defaults('cmro2.gammaT') ioi_get_defaults('cmro2.gammaR')};
% now put these defaults into the optargs cell array, and overwrite the ones
% specified in varargin.
optargs(1:numvarargs)	= varargin;
% Place optional args in memorable variable names
[gammaT gammaR]         = optargs{:};
% ------------------------------------------------------------------------------

try
    % Total hemoglobin
    dataHbT = dataHbO + dataHbR;
    % Baseline HbT
    HbT0 = ioi_get_defaults('conc1.baseline_hbt');
    % Baseline HbR
    HbR0 = ioi_get_defaults('conc1.baseline_hbr');
    % Baseline flow
    F0 = 1;
    % Estimation of cerebral metabolic rate of oxygen.
    cmro2 = (1 + dataFlow/F0).*(1 + gammaR*dataHbR/HbR0)./(1 + gammaT*dataHbT/HbT0) - 1;
catch exception
    disp(exception.identifier)
    disp(exception.stack(1))
    cmro2 = [];
end

% EOF
