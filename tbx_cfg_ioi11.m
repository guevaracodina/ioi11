function ioi11 = tbx_cfg_ioi11
%_______________________________________________________________________________
% Copyright (C) 2011 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________

addpath(fileparts(which(mfilename)));

%-----------------------------------------------------------------------
ioi11        = cfg_choice;
ioi11.name   = 'IOI11';
ioi11.tag    = 'ioi11'; % Careful, this tag ioi11 must be the same as
% the name of the toolbox and when called by spm_jobman in ioi11.m
ioi11.values = {ioi_shrink_cfg ioi_msioi_cfg ioi_concentrations_cfg ...
     ioi_flow_cfg ioi_cmro2_cfg ioi_fcIOS_cfg ioi_create_onsets_cfg ioi_cine2D_cfg ioi_cine_display_cfg ...
     ioi_create_roi_cfg ioi_extract_roi_time_series_cfg ioi_stim_mean_cfg ioi_stim_mean_image_cfg ...
     ioi_seizure_mean_roi_cfg ioi_GLM_on_ROI_cfg ioi_GLM_cfg ioi_HDM_cfg ioi_HDM_all_cfg ...
     ioi_SCKS_cfg ioi_ROC_cfg}; 

 % EOF
 
