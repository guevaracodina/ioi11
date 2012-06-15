function stim_mean1 = ioi_stim_mean_cfg
%_______________________________________________________________________
% Copyright (C) 2011 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal

IOImat = ioi_dfg_IOImat(1);
redo1 = ioi_dfg_redo(0);
ROImat = ioi_dfg_ROImat(1);
IOImatCopyChoice = ioi_dfg_IOImatCopyChoice('Mean');
ROI_choice = ioi_cfg_ROI_choice;
session_choice = ioi_cfg_session_choice;
[window_before window_after window_offset] = ioi_dfg_window(3,20,0);

normalize_choice      = cfg_menu;
normalize_choice.tag  = 'normalize_choice';
normalize_choice.name = 'Normalization choice';
normalize_choice.labels = {'Median over window before','Time zero','Mean'};
normalize_choice.values = {1,2,3};
normalize_choice.val  = {1};
normalize_choice.help = {'Normalization choice. In one test,'
    'The mean standard deviation was higher by 10% or more'
    'when using time zero as the baseline, compared to taking '
    'an average (median) over the window before.'}';

which_onset_type = ioi_cfg_which_onset_type;
remove_stims = ioi_cfg_remove_stims;
use_stims = ioi_cfg_use_stims;
[generate_figures save_figures] = ioi_cfg_generate_figures;

extract_HRF      = cfg_menu;
extract_HRF.tag  = 'extract_HRF';
extract_HRF.name = 'Extract HRF';
extract_HRF.labels = {'Yes','No'};
extract_HRF.values = {1,0};
extract_HRF.val  = {1};
extract_HRF.help = {'Note: this should only be used for impulse-like responses,'
    'meaning short stimulations of 1 second or less.'
    'Extract 6 coefficients of hemodynamic response function,' 
    'by fitting the average curves to a difference of two gamma functions'}';

include_nlinfit      = cfg_menu;
include_nlinfit.tag  = 'include_nlinfit';
include_nlinfit.name = 'Include nlinfit';
include_nlinfit.labels = {'Yes','No'};
include_nlinfit.values = {1,0};
include_nlinfit.val  = {0};
include_nlinfit.help = {'Include nlinfit (option eventually to plot biophysical model response).'}';

fit_3_gamma      = cfg_menu;
fit_3_gamma.tag  = 'fit_3_gamma';
fit_3_gamma.name = 'Fit 3 gamma functions';
fit_3_gamma.labels = {'Yes','No'};
fit_3_gamma.values = {1,0};
fit_3_gamma.val  = {0};
fit_3_gamma.help = {'For Expectation-Maximization Fit, use 3 gamma functions, instead of 2'}';

generate_global      = cfg_menu;
generate_global.tag  = 'generate_global';
generate_global.name = 'Generate global data';
generate_global.labels = {'Yes','No'};
generate_global.values = {1,0};
generate_global.val  = {0};
generate_global.help = {'Generate data averaged over all sessions.'}';

add_error_bars      = cfg_menu;
add_error_bars.tag  = 'add_error_bars';
add_error_bars.name = 'Add error bars';
add_error_bars.labels = {'Yes','No'};
add_error_bars.values = {1,0};
add_error_bars.val  = {0};
add_error_bars.help = {'Add error bars.'}';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hpf_butter = ioi_dfg_hpf_butter(1,0.01,3);

% ---------------------------------------------------------------------
% lpf Low-pass filter
% ---------------------------------------------------------------------
fwhm1      = cfg_entry;
fwhm1.tag  = 'fwhm1';
fwhm1.name = 'FWHM in seconds';
fwhm1.val = {0.67};
fwhm1.strtype = 'r';  
fwhm1.num     = [1 1]; 
fwhm1.help    = {'FWHM in seconds.'}; 

lpf_gauss_On         = cfg_branch;
lpf_gauss_On.tag     = 'lpf_gauss_On';
lpf_gauss_On.name    = 'Gaussian LP filter';
lpf_gauss_On.val     = {fwhm1}; 
lpf_gauss_On.help    = {'Gaussian low-pass filter '
    '(applied forward then backward so that it does not create a time shift)'}';

lpf_Off         = cfg_branch;
lpf_Off.tag     = 'lpf_Off';
lpf_Off.name    = 'LP filter off';
lpf_Off.val     = {}; 
lpf_Off.help    = {'Low pass filter turned off.'};

lpf_choice      = cfg_choice;
lpf_choice.tag  = 'lpf_choice';
lpf_choice.name = 'Choose Low Pass Filter';
lpf_choice.values = {lpf_gauss_On lpf_Off};
lpf_choice.val = {lpf_Off};
lpf_choice.help = {'Choose whether to include a Low Pass Filter.'
        'Parameters'}';
    
remove_segment_drift      = cfg_menu;
remove_segment_drift.tag  = 'remove_segment_drift';
remove_segment_drift.name = 'Remove segment drift';
remove_segment_drift.labels = {'Yes','No'};
remove_segment_drift.values = {1,0};
remove_segment_drift.val  = {0};
remove_segment_drift.help = {'Remove linear drift separately on each segment.'}';

remove_stims_SD      = cfg_menu;
remove_stims_SD.tag  = 'remove_stims_SD';
remove_stims_SD.name = 'Remove stims on SD criterion';
remove_stims_SD.labels = {'Yes','No'};
remove_stims_SD.values = {1,0};
remove_stims_SD.val  = {0};
remove_stims_SD.help = {'Remove stims on SD criterion.'
    'This works in addition to and after other stims having been removed.'}';

IC = ioi_dfg_include_colors(0,1,1,1,0);

% Executable Branch
stim_mean1      = cfg_exbranch;       % This is the branch that has information about how to run this module
stim_mean1.name = 'Average stimulations';             % The display name
stim_mean1.tag  = 'stim_mean1'; %Very important: tag is used when calling for execution
stim_mean1.val  = {IOImat ROImat redo1 IOImatCopyChoice session_choice ...
    ROI_choice window_after window_before window_offset normalize_choice hpf_butter ...
    lpf_choice IC which_onset_type ...
    remove_stims use_stims remove_stims_SD extract_HRF fit_3_gamma include_nlinfit ...
    generate_global generate_figures save_figures add_error_bars ...
    remove_segment_drift };    % The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs
stim_mean1.prog = @ioi_stim_mean_run;  % A function handle that will be called with the harvested job to run the computation
stim_mean1.vout = @ioi_cfg_vout_stim_mean; % A function handle that will be called with the harvested job to determine virtual outputs
stim_mean1.help = {'Calculate average over stimulations.'};

return

%make IOI.mat available as a dependency
function vout = ioi_cfg_vout_stim_mean(job)
vout = cfg_dep;                     % The dependency object
vout.sname      = 'IOI.mat';       % Displayed dependency name
vout.src_output = substruct('.','IOImat'); %{1}); %,'IOImat');
%substruct('()',{1}); % The output subscript reference. This could be any reference into the output variable created during computation
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
