function ioi11
% Based on spm_vbm8 from vbm8
% Toolbox wrapper to call ioi11 functions
%_______________________________________________________________________
% Copyright (C) 2011 LIOM Laboratoire d'Imagerie Optique et Mol�culaire
%                    �cole Polytechnique de Montr�al
%______________________________________________________________________
rev = '$Rev: 1 $';
SPMid = spm('FnBanner',mfilename,rev);
[~, Fgraph, ~] = spm('FnUIsetup','IOI11');
spm_help('!ContextHelp',mfilename);
spm_help('!Disp','ioi11.man','',Fgraph,'IOI11');
fig = spm_figure('GetWin','Interactive');
h0  = uimenu(fig,...
    'Label',	'IOI11',...
    'Separator',	'on',...
    'Tag',		'IOI11',...
    'HandleVisibility','on');

h01 = uimenu(h0,...
	'Label',	'Shrink specific colors -- old format only',...
	'Separator',	'off',...
	'Tag',		'Shrink specific colors -- old format only',...
	'CallBack','spm_jobman(''interactive'','''',''spm.tools.ioi11.shrink1'');',...
	'HandleVisibility','on');

h1  = uimenu(h0,...
	'Label',	'Read multispectral IOI data',...
	'Separator',	'off',...
	'Tag',		'Read multispectral IOI data',...
	'CallBack','spm_jobman(''interactive'','''',''spm.tools.ioi11.msioi1'');',...
	'HandleVisibility','on');

h2  = uimenu(h0,...
	'Label',	'Compute Concentrations',...
	'Separator',	'off',...
	'Tag',		'Compute Concentrations',...
    'CallBack','spm_jobman(''interactive'','''',''spm.tools.ioi11.conc1'');',...
	'HandleVisibility','on');

h3  = uimenu(h0,...
	'Label',	'Compute Flow',...
	'Separator',	'off',...
	'Tag',		'Compute Flow',...
    'CallBack','spm_jobman(''interactive'','''',''spm.tools.ioi11.flow1'');',...
	'HandleVisibility','on');
% For all the h4? items see ioi_fcIOS_cfg.m //EGC
h4  = uimenu(h0,...
    'Label',	'Functional Connectivity Mapping (fcIOS) Utilities',...
    'Separator',	'on',...
    'Tag',		'Functional Connectivity Mapping (fcIOS) Utilities',...
    'HandleVisibility','on');

h4a  = uimenu(h4,...
	'Label',	'Create Brain Mask',...
	'Separator',	'off',...
	'Tag',		'Create Brain Mask',...
    'CallBack','spm_jobman(''interactive'','''',''spm.tools.ioi11.fcIOS.mask1'');',...
	'HandleVisibility','on');

h4b  = uimenu(h4,...
	'Label',	'Create ROI (seeds) for fcIOS',...
	'Separator',	'off',...
	'Tag',		'Create Seeds for fcIOS',...
    'CallBack','spm_jobman(''interactive'','''',''spm.tools.ioi11.create_roi1'');',...
	'HandleVisibility','on');

h4c  = uimenu(h4,...
	'Label',	'Spatial 2-D low-pass filtering',...
	'Separator',	'off',...
	'Tag',		'Spatial 2-D low-pass filtering',...
    'CallBack','spm_jobman(''interactive'','''',''spm.tools.ioi11.fcIOS.spatial_LPF1'');',...
	'HandleVisibility','on');

h4d  = uimenu(h4,...
	'Label',	'Extract ROI (seeds) for fcIOS',...
	'Separator',	'off',...
	'Tag',		'Extract Seeds for fcIOS',...
    'CallBack','spm_jobman(''interactive'','''',''spm.tools.ioi11.extract_roi1'');',...
	'HandleVisibility','on');

h4e  = uimenu(h4,...
	'Label',	'Temporal filtering and downsampling of seeds & whole image series',...
	'Separator',	'off',...
	'Tag',		'Temporal filtering and downsampling seeds & whole image series',...
    'CallBack','spm_jobman(''interactive'','''',''spm.tools.ioi11.fcIOS.filtdown1'');',...
	'HandleVisibility','on');

h4f  = uimenu(h4,...
	'Label',	'GLM regression of global brain signal',...
	'Separator',	'off',...
	'Tag',		'GLM regression of global brain signal',...
    'CallBack','spm_jobman(''interactive'','''',''spm.tools.ioi11.fcIOS.fc_GLM_on_ROI1'');',...
	'HandleVisibility','on');

h4g  = uimenu(h4,...
	'Label',	'Functional connectivity (fcIOS) mapping',...
	'Separator',	'off',...
	'Tag',		'Functional connectivity (fcIOS) mapping',...
    'CallBack','spm_jobman(''interactive'','''',''spm.tools.ioi11.fcIOS.correlation_map1'');',...
	'HandleVisibility','on');

h4h  = uimenu(h4,...
	'Label',	'Bilateral correlation group comparison',...
	'Separator',	'on',...
	'Tag',		'Bilateral correlation group comparison',...
    'CallBack','spm_jobman(''interactive'','''',''spm.tools.ioi11.fcIOS.group_corr1'');',...
	'HandleVisibility','on');

h4i  = uimenu(h4,...
	'Label',	'Update elinfo data',...
	'Separator',	'on',...
	'Tag',		'Update elinfo data',...
    'CallBack','spm_jobman(''interactive'','''',''spm.tools.ioi11.fcIOS.elinfo1'');',...
	'HandleVisibility','on');

h5  = uimenu(h0,...
	'Label',	'Create onsets',...
	'Separator',	'on',...
	'Tag',		'Create onsets',...
    'CallBack','spm_jobman(''interactive'','''',''spm.tools.ioi11.create_onsets1'');',...
	'HandleVisibility','on');

h6  = uimenu(h0,...
	'Label',	'Cine 2D',...
	'Separator',	'off',...
	'Tag',		'Cine 2D',...
    'CallBack','spm_jobman(''interactive'','''',''spm.tools.ioi11.cine2D1'');',...
	'HandleVisibility','on');

h6b  = uimenu(h0,...
	'Label',	'Display Cine',...
	'Separator',	'off',...
	'Tag',		'Display Cine',...
    'CallBack','spm_jobman(''interactive'','''',''spm.tools.ioi11.cine_display1'');',...
	'HandleVisibility','on');

h7  = uimenu(h0,...
	'Label',	'Create ROI',...
	'Separator',	'off',...
	'Tag',		'Create ROI',...
    'CallBack','spm_jobman(''interactive'','''',''spm.tools.ioi11.create_roi1'');',...
	'HandleVisibility','on');

h8  = uimenu(h0,...
	'Label',	'Extract ROI',...
	'Separator',	'off',...
	'Tag',		'Extract ROI',...
    'CallBack','spm_jobman(''interactive'','''',''spm.tools.ioi11.extract_roi1'');',...
	'HandleVisibility','on');

h9  = uimenu(h0,...
	'Label',	'Average stimulations (on ROIs)',...
	'Separator',	'off',...
	'Tag',		'Average stimulations (on ROIs)',...
    'CallBack','spm_jobman(''interactive'','''',''spm.tools.ioi11.stim_mean1'');',...
	'HandleVisibility','on');

h9b  = uimenu(h0,...
	'Label',	'Average stimulations (on images)',...
	'Separator',	'off',...
	'Tag',		'Average stimulations (on images)',...
    'CallBack','spm_jobman(''interactive'','''',''spm.tools.ioi11.stim_mean_image1'');',...
	'HandleVisibility','on');

h10 = uimenu(h0,...
	'Label',	'GLM (on ROIs)',...
	'Separator',	'off',...
	'Tag',		'GLM (on ROIs)',...
    'CallBack','spm_jobman(''interactive'','''',''spm.tools.ioi11.glm_roi1'');',...
	'HandleVisibility','on');

h10b = uimenu(h0,...
	'Label',	'Full GLM (images or ROIs)',...
	'Separator',	'off',...
	'Tag',		'Full GLM (images or ROIs)',...
    'CallBack','spm_jobman(''interactive'','''',''spm.tools.ioi11.glm1'');',...
	'HandleVisibility','on');

h11  = uimenu(h0,...
	'Label',	'Hemodynamic Modelling',...
	'Separator',	'off',...
	'Tag',		'Hemodynamic Modelling',...
    'CallBack','spm_jobman(''interactive'','''',''spm.tools.ioi11.hdm1'');',...
	'HandleVisibility','on');

h11b  = uimenu(h0,...
	'Label',	'Hemodynamic Modelling (images or ROIs)',...
	'Separator',	'off',...
	'Tag',		'Hemodynamic Modelling (images or ROIs)',...
    'CallBack','spm_jobman(''interactive'','''',''spm.tools.ioi11.hdm_all1'');',...
	'HandleVisibility','on');

h12  = uimenu(h0,...
	'Label',	'SCKS Deconvolution',...
	'Separator',	'off',...
	'Tag',		'SCKS Deconvolution',...
    'CallBack','spm_jobman(''interactive'','''',''spm.tools.ioi11.SCKS1'');',...
	'HandleVisibility','on');

h13  = uimenu(h0,...
	'Label',	'ROC curves for SCKS',...
	'Separator',	'off',...
	'Tag',		'ROC curves for SCKS',...
    'CallBack','spm_jobman(''interactive'','''',''spm.tools.ioi11.ROC1'');',...
	'HandleVisibility','on');

h14  = uimenu(h0,...
    'Label',	'Full IOI study (stims)',...
    'Separator',	'on',...
    'Tag',		'Full IOI study (stims)',...
    'CallBack',  'spm_jobman(''interactive'',''ioi_fullstudy_stims_cfg.m'');',...
    'HandleVisibility','on');

h15  = uimenu(h0,...
    'Label',	'Full IOI study (spikes)',...
    'Separator',	'off',...
    'Tag',		'Full IOI study (spikes)',...
    'CallBack',  'spm_jobman(''interactive'',''ioi_fullstudy_spikes_cfg.m'');',...
    'HandleVisibility','on');

h16  = uimenu(h0,...
    'Label',	'IOI preprocessing (first 3 modules)',...
    'Separator',	'off',...
    'Tag',		'IOI preprocessing (first 3 modules)',...
    'CallBack',  'spm_jobman(''interactive'',''ioi_preprocessing_cfg.m'');',...
    'HandleVisibility','on');

h17  = uimenu(h0,...
    'Label',	'IOI study of stims (after preprocessing)',...
    'Separator',	'off',...
    'Tag',		'IOI study of stims (after preprocessing)',...
    'CallBack',  'spm_jobman(''interactive'',''ioi_partial_stims_cfg.m'');',...
    'HandleVisibility','on');

h18  = uimenu(h0,...
    'Label',	'IOI study of spikes (after preprocessing)',...
    'Separator',	'off',...
    'Tag',		'IOI study of spikes (after preprocessing)',...
    'CallBack',  'spm_jobman(''interactive'',''ioi_partial_spikes_cfg.m'');',...
    'HandleVisibility','on');
