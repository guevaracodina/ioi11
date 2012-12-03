function ioi_flow_defaults
% Sets the defaults for computing blood flow 
% (decorrelation velocity of the speckle imaging)
%_______________________________________________________________________________
% Copyright (C) 2011 LIOM Laboratoire d'Imagerie Optique et Moleculaire
%                    Ecole Polytechnique de Montreal
%_______________________________________________________________________________

global ioi11

ioi11.flow1.T           = 0.2;  % Integration time of the CCD camera (Units???)
ioi11.flow1.window_size = 5;    % Speckle contrast window size (In pixels)
ioi11.flow1.lambda      = 785e-9;% Wavelength of the IR laser diode (in m)

% EOF


