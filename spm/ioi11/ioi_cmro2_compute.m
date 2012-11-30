function cmro2 = ioi_cmro2_compute(Flow, HbO, HbR, varargin)
% Computes cerebral metabolic rate of oxygen (CMRO2) from blood flow, HbO and HbR data.
% SYNTAX
% cmro2 = ioi_cmro2_compute(Flow, HbO, HbR, gammaT, gammaR)
% INPUTS
% Flow      Percent change of blood flow.
% HbO       Oxygenated hemoglobin concentrations
% HbR       Deoxygenated hemoglobin concentrations
% gammaT    [OPTIONAL] Vascular weighting constant gammaT (default: 1)
% gammaR    [OPTIONAL] Vascular weighting constant gammaR (default: 1)
% HbT0      [OPTIONAL] Total hemoglobin concentration (default: 100mM)
% HbR0      [OPTIONAL] Deoxy-hemoglobin baseline concentration (default: 40mM)
% OUTPUT
% cmro2     Yields \DeltaCMRO2/CMRO2_0
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
numvarargs                  = length(varargin);
if numvarargs > 4
    error('ioi_cmro2_compute:TooManyInputs', ...
        'Requires at most 4 optional inputs');
end
% set defaults for optional inputs
optargs                     = { ioi_get_defaults('cmro2.gammaT') ioi_get_defaults('cmro2.gammaR')...
                            ioi_get_defaults('conc1.baseline_hbt') ioi_get_defaults('conc1.baseline_hbr')};
% now put these defaults into the optargs cell array, and overwrite the ones
% specified in varargin.
optargs(1:numvarargs)       = varargin;
% Place optional args in memorable variable names
[gammaT gammaR HbT0 HbR0]   = optargs{:};
% ------------------------------------------------------------------------------

% Wavelength of the IR laser diode
lambda = 785e-9;
try
    % Using lambda = 785nm the decorrelation velocity is found from flow images,
    % technicallly we do not have \DeltaF/F_0, but since we high-pass filtered
    % flow time course it is equivalent to \DeltaF
    Flow = Flow*lambda/2*pi;
    % Total hemoglobin
    HbT = HbO + HbR;
    % Baseline flow?
    % F0 = 1;
    % Estimation of cerebral metabolic rate of oxygen.
    cmro2 = (1 + Flow).*(1 + gammaR*HbR/HbR0)./(1 + gammaT*HbT/HbT0) - 1;
catch exception
    disp(exception.identifier)
    disp(exception.stack(1))
    cmro2 = [];
end

% EOF
