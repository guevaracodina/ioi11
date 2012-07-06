function fcIOS = ioi_fcIOS_cfg
% Graphical interface configuration function for functional connectivity mapping
% with intrinsic optical signal (fcIOS) module
%_______________________________________________________________________________
% Copyright (C) 2012 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________

fcIOS        = cfg_choice;
fcIOS.name   = 'Functional connectivity mapping (fcIOS)';
fcIOS.tag    = 'fcIOS';
fcIOS.values = {ioi_brainmask_cfg ioi_spatial_LPF_cfg ioi_create_roi_cfg...
    ioi_extract_roi_time_series_cfg ioi_filtdown_cfg ioi_fc_GLM_on_ROI_cfg ...
    ioi_correlation_map_cfg};
fcIOS.help   = {'These modules perform resting-state functional connectivity mapping with intrinsic optical signal (fcIOS).'};
