function mri12 = tbx_cfg_mri12
%_______________________________________________________________________
% Copyright (C) 2012 LIOM Laboratoire d'Imagerie Optique et Mol�culaire
%                    �cole Polytechnique de Montr�al
%______________________________________________________________________

addpath(fileparts(which(mfilename)));

%-----------------------------------------------------------------------
mri12        = cfg_choice;
mri12.name   = 'mri12';
mri12.tag    = 'mri12'; 
mri12.values = {}; 