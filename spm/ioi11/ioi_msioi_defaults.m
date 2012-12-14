function ioi_msioi_defaults
%_______________________________________________________________________
% Copyright (C) 2011 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%______________________________________________________________________
% Sets the defaults for reading binary files

global ioi11

% main path
ioi11.msioi1.output_path_select.output_path = fullfile('D:/Users/');
ioi11.msioi1.force_redo = 0;
% Shrink factors to reduce impact on memory
ioi11.msioi1.shrink_x = 2;
ioi11.msioi1.shrink_y = 2;

% Contrast window size for speckle
ioi11.msioi1.window_size = 5;
