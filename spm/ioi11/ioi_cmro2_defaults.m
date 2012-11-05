function ioi_cmro2_defaults
% Sets the default values for the CMRO2 constants. 
% gammaR & gammaT are vascular weighting constant accounting for the fact that
% changes in hemoglobin are recorded in arterial and venous compartment while
% our computation is defined only in the venous compartment. The more
% physiologically plausible range is around 1 (0.75–1.25).
%_______________________________________________________________________________
% Copyright (C) 2012 LIOM Laboratoire d'Imagerie Optique et Moleculaire
%                    Ecole Polytechnique de Montreal
%_______________________________________________________________________________

global ioi11
ioi11.cmro2.gammaR = 1;
ioi11.cmro2.gammaT = 1;

% EOF
