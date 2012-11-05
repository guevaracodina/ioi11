function cmro2 = ioi_cmro2_compute()
% Computes CMRO2 from HbO, HbR and blood flow data.
% M. Jones, J. Berwick, D. Johnston, and J. Mayhew, “Concurrent optical imaging
% spectroscopy and laser-Doppler flowmetry: the relationship between blood flow,
% oxygenation, and volume in rodent barrel cortex,” Neuroimage, vol. 13, no. 6
% Pt 1, pp. 1002–1015, Jun. 2001.
%_______________________________________________________________________________
% Copyright (C) 2012 LIOM Laboratoire d'Imagerie Optique et Moleculaire
%                    Ecole Polytechnique de Montreal
%_______________________________________________________________________________
try
    cmro2 = true;
catch exception
    disp(exception.identifier)
    disp(exception.stack(1))
    cmro2 = [];
end

% EOF
