function ioi11
rev = '$Rev: 1 $';
SPMid = spm('FnBanner',mfilename,rev);
[dummy,Fgraph,dummy2] = spm('FnUIsetup','IOI11');
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

h4  = uimenu(h0,...
	'Label',	'Create onsets',...
	'Separator',	'off',...
	'Tag',		'Create onsets',...
    'CallBack','spm_jobman(''interactive'','''',''spm.tools.ioi11.create_onsets1'');',...
	'HandleVisibility','on');

h5  = uimenu(h0,...
	'Label',	'Cine 2D',...
	'Separator',	'off',...
	'Tag',		'Cine 2D',...
    'CallBack','spm_jobman(''interactive'','''',''spm.tools.ioi11.cine2D1'');',...
	'HandleVisibility','on');

h5b  = uimenu(h0,...
	'Label',	'Display Cine',...
	'Separator',	'off',...
	'Tag',		'Display Cine',...
    'CallBack','spm_jobman(''interactive'','''',''spm.tools.ioi11.cine_display1'');',...
	'HandleVisibility','on');

h6  = uimenu(h0,...
	'Label',	'Select ROI',...
	'Separator',	'off',...
	'Tag',		'Select ROI',...
    'CallBack','spm_jobman(''interactive'','''',''spm.tools.ioi11.create_roi1'');',...
	'HandleVisibility','on');

h7  = uimenu(h0,...
	'Label',	'Extract ROI',...
	'Separator',	'off',...
	'Tag',		'Extract ROI',...
    'CallBack','spm_jobman(''interactive'','''',''spm.tools.ioi11.extract_roi1'');',...
	'HandleVisibility','on');

h8  = uimenu(h0,...
	'Label',	'Average stimulations',...
	'Separator',	'off',...
	'Tag',		'Average stimulations',...
    'CallBack','spm_jobman(''interactive'','''',''spm.tools.ioi11.stim_mean1'');',...
	'HandleVisibility','on');

h9  = uimenu(h0,...
	'Label',	'GLM on ROIs',...
	'Separator',	'off',...
	'Tag',		'GLM on ROIs',...
    'CallBack','spm_jobman(''interactive'','''',''spm.tools.ioi11.glm_roi1'');',...
	'HandleVisibility','on');

h9b  = uimenu(h0,...
	'Label',	'Full GLM (images or ROIs)',...
	'Separator',	'off',...
	'Tag',		'Full GLM (images or ROIs)',...
    'CallBack','spm_jobman(''interactive'','''',''spm.tools.ioi11.glm1'');',...
	'HandleVisibility','on');

h10  = uimenu(h0,...
	'Label',	'Hemodynamic Modelling',...
	'Separator',	'off',...
	'Tag',		'Hemodynamic Modelling',...
    'CallBack','spm_jobman(''interactive'','''',''spm.tools.ioi11.hdm1'');',...
	'HandleVisibility','on');

h10b  = uimenu(h0,...
	'Label',	'Hemodynamic Modelling (images or ROIs)',...
	'Separator',	'off',...
	'Tag',		'Hemodynamic Modelling (images or ROIs)',...
    'CallBack','spm_jobman(''interactive'','''',''spm.tools.ioi11.hdm_all1'');',...
	'HandleVisibility','on');

h11  = uimenu(h0,...
	'Label',	'SCKS Deconvolution',...
	'Separator',	'off',...
	'Tag',		'SCKS Deconvolution',...
    'CallBack','spm_jobman(''interactive'','''',''spm.tools.ioi11.SCKS1'');',...
	'HandleVisibility','on');

h12  = uimenu(h0,...
	'Label',	'ROC curves for SCKS',...
	'Separator',	'off',...
	'Tag',		'ROC curves for SCKS',...
    'CallBack','spm_jobman(''interactive'','''',''spm.tools.ioi11.ROC1'');',...
	'HandleVisibility','on');

h13  = uimenu(h0,...
    'Label',	'Full IOI study (stims)',...
    'Separator',	'on',...
    'Tag',		'Full IOI study (stims)',...
    'CallBack',  'spm_jobman(''interactive'',''ioi_fullstudy_stims_cfg.m'');',...
    'HandleVisibility','on');

h14  = uimenu(h0,...
    'Label',	'Full IOI study (spikes)',...
    'Separator',	'off',...
    'Tag',		'Full IOI study (spikes)',...
    'CallBack',  'spm_jobman(''interactive'',''ioi_fullstudy_spikes_cfg.m'');',...
    'HandleVisibility','on');

h15  = uimenu(h0,...
    'Label',	'IOI preprocessing (first 3 modules)',...
    'Separator',	'off',...
    'Tag',		'IOI preprocessing (first 3 modules)',...
    'CallBack',  'spm_jobman(''interactive'',''ioi_preprocessing_cfg.m'');',...
    'HandleVisibility','on');

h16  = uimenu(h0,...
    'Label',	'IOI study of stims (after preprocessing)',...
    'Separator',	'off',...
    'Tag',		'IOI study of stims (after preprocessing)',...
    'CallBack',  'spm_jobman(''interactive'',''ioi_partial_stims_cfg.m'');',...
    'HandleVisibility','on');

h17  = uimenu(h0,...
    'Label',	'IOI study of spikes (after preprocessing)',...
    'Separator',	'off',...
    'Tag',		'IOI study of spikes (after preprocessing)',...
    'CallBack',  'spm_jobman(''interactive'',''ioi_partial_spikes_cfg.m'');',...
    'HandleVisibility','on');