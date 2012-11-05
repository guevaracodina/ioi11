function stim_mean1 = ioi_stim_mean_cfg
%_______________________________________________________________________
% Copyright (C) 2011 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal

IOImat = ioi_dfg_IOImat(1);
redo1 = ioi_dfg_redo(0);
ROImat = ioi_dfg_ROImat(1);
IOImatCopyChoice = ioi_dfg_IOImatCopyChoice('Mean');
ROI_choice = ioi_dfg_ROI_choice;
session_choice = ioi_dfg_session_choice;
[window_before window_after window_offset] = ioi_dfg_window(3,20,0);

normalize_choice      = cfg_menu;
normalize_choice.tag  = 'normalize_choice';
normalize_choice.name = 'Additive normalization choice';
normalize_choice.labels = {'Median over window before','Time zero','Mean'};
normalize_choice.values = {1,2,3};
normalize_choice.val  = {1};
normalize_choice.help = {'Normalization choice. In one test,'
    'The mean standard deviation was higher by 10% or more'
    'when using time zero as the baseline, compared to taking '
    'an average (median) over the window before.'}';

Weigh_amplitude_of_threshold      = cfg_menu;
Weigh_amplitude_of_threshold.tag  = 'Weigh_amplitude_of_threshold';
Weigh_amplitude_of_threshold.name = 'Weigh_threshold';
Weigh_amplitude_of_threshold.labels = {'Yes','No'};
Weigh_amplitude_of_threshold.values = {1,0};
Weigh_amplitude_of_threshold.val  = {0};
Weigh_amplitude_of_threshold.help = {'Note: this is done when the amplitude of thrshold was changed during one session.'
    'When the amplitude of the threshold was changed slect Yes'
    'When the amplitude of the threshold was not changed select No'}';

mult_normalize_choice      = cfg_menu;
mult_normalize_choice.tag  = 'mult_normalize_choice';
mult_normalize_choice.name = 'Multiplicative normalization choice';
mult_normalize_choice.labels = {'Divide by constant','Divide by baseline of individual stim'};
mult_normalize_choice.values = {1,0};
mult_normalize_choice.val  = {0};
mult_normalize_choice.help = {'For oxy, deoxy and total hemoglobin only, '
    'choose to divide by a constant (e.g. 40, 60 or 100 microMolar'
    'or to divide by a different value for each stim, '
    'namely its baseline value in microMolar,'
    'however, this value may be too volatile and give nonsensical results.'}';

which_onset_type = ioi_dfg_which_onset_type;
remove_stims = ioi_dfg_remove_stims;
use_stims = ioi_dfg_use_stims;
[generate_figures save_figures] = ioi_dfg_generate_figures;

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

hpf_butter = ioi_dfg_hpf_butter(1,0.01,3);
lpf_choice = ioi_dfg_lpf_choice(0,0.67);
    
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

std_choice         = cfg_entry; 
std_choice.name    = 'Choice of standard deviation threshold';
std_choice.tag     = 'std_choice';       
std_choice.strtype = 'r';
std_choice.num     = [1 1];     
std_choice.val     = {1};
std_choice.help    = {'Choice of standard deviation threshold.'}';

IC = ioi_dfg_include_colors(0,1,1,1,0);

make_stim_figures      = cfg_menu;
make_stim_figures.tag  = 'make_stim_figures';
make_stim_figures.name = 'Make figures with all stims on the same figure';
make_stim_figures.labels = {'Yes','No'};
make_stim_figures.values = {1,0};
make_stim_figures.val  = {1};
make_stim_figures.help = {'Make figures with all stims on the same figure.'}';

make_timeCourse_figures      = cfg_menu;
make_timeCourse_figures.tag  = 'make_timeCourse_figures';
make_timeCourse_figures.name = 'Make figures of time courses';
make_timeCourse_figures.labels = {'Yes','No'};
make_timeCourse_figures.values = {1,0};
make_timeCourse_figures.val  = {1};
make_timeCourse_figures.help = {'Make figures of time courses.'}';

make_eachStim_figures      = cfg_menu;
make_eachStim_figures.tag  = 'make_eachStim_figures';
make_eachStim_figures.name = 'Make a figure for each stim';
make_eachStim_figures.labels = {'Yes','No'};
make_eachStim_figures.values = {1,0};
make_eachStim_figures.val  = {0};
make_eachStim_figures.help = {'Make a figure for each stim.'}';

write_excel_text      = cfg_menu;
write_excel_text.tag  = 'write_excel_text';
write_excel_text.name = 'Write data to Excel or to .txt files';
write_excel_text.labels = {'None','Excel','.txt','Both Excel and .txt'};
write_excel_text.values = {0,1,2,3};
write_excel_text.val  = {0};
write_excel_text.help = {'Make a figure for each stim.'}';

write_points_before      = cfg_menu;
write_points_before.tag  = 'write_points_before';
write_points_before.name = 'Write to Excel the time points before the onsets';
write_points_before.labels = {'Yes','No'};
write_points_before.values = {1,0};
write_points_before.val  = {0};
write_points_before.help = {'Write to Excel the time points before the onsets,'
    'used for the baseline.'}';

% Executable Branch
stim_mean1      = cfg_exbranch;       % This is the branch that has information about how to run this module
stim_mean1.name = 'Average stimulations';             % The display name
stim_mean1.tag  = 'stim_mean1'; %Very important: tag is used when calling for execution
stim_mean1.val  = {IOImat ROImat redo1 IOImatCopyChoice session_choice ...
    ROI_choice window_after window_before window_offset normalize_choice Weigh_amplitude_of_threshold ...
    mult_normalize_choice hpf_butter ...
    lpf_choice IC which_onset_type ...
    remove_stims use_stims remove_stims_SD std_choice extract_HRF fit_3_gamma include_nlinfit ...
    generate_global generate_figures save_figures make_timeCourse_figures ...
    make_stim_figures make_eachStim_figures add_error_bars ...
    remove_segment_drift write_excel_text write_points_before};    % The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs
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
