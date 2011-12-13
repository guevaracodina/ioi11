function stim_mean1 = ioi_stim_mean_cfg
%_______________________________________________________________________
% Copyright (C) 2011 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal

IOImat         = cfg_files; %Select NIRS.mat for this subject 
IOImat.name    = 'Select IOI.mat'; % The displayed name
IOImat.tag     = 'IOImat';       %file names
IOImat.filter = 'mat';
IOImat.ufilter = '^IOI.mat$';    
IOImat.num     = [1 Inf];     % Number of inputs required 
IOImat.help    = {'Select IOImat dependency if available. '
    'Otherwise, for each subject, select IOI.mat.'}'; % help text displayed

ROImat         = cfg_files; 
ROImat.name    = 'Select ROI.mat'; 
ROImat.tag     = 'ROImat';      
ROImat.filter = 'mat';
ROImat.ufilter = '^ROI.mat$';
ROImat.val     = {''};
ROImat.num     = [1 Inf];    
ROImat.help    = {'Optional: Select ROImat. This allows working on ROI data even'
    'if the paths are not correct in IOI.mat. If not specified, the ROI.mat '
    'specified in IOI.mat will be used.'
    'If several subjects are run, and ROImat is explicitly specified, then'
    'it should be specified for all subjects, in the same order as the list of IOI.mat'}'; 

redo1      = cfg_menu;
redo1.tag  = 'force_redo';
redo1.name = 'Force processing';
redo1.labels = {'False','True'};
redo1.values = {0,1};
redo1.val  = {0};
redo1.help = {'Force redoing this processing even when it has been done already.'
    'Use option below for treatment of previous ROIs.'}';

IOImatOverwrite         = cfg_branch;
IOImatOverwrite.tag     = 'IOImatOverwrite';
IOImatOverwrite.name    = 'Overwrite IOI.mat structure'; 
IOImatOverwrite.help    = {'Will not copy IOI structure.'
            'This will write over the previous NIRS.mat'}';

NewIOIdir         = cfg_entry;
NewIOIdir.name    = 'New directory for IOI.mat';
NewIOIdir.tag     = 'NewIOIdir';       
NewIOIdir.strtype = 's';
NewIOIdir.val{1}    = 'Mean';
NewIOIdir.num     = [1 Inf];     
NewIOIdir.help    = {'Directory for IOI.mat.'}'; 

IOImatCopy         = cfg_branch;
IOImatCopy.tag     = 'IOImatCopy';
IOImatCopy.name    = 'Create new directory and copy IOI structure there'; 
IOImatCopy.val     = {NewIOIdir};
IOImatCopy.help    = {'Create new directory and copy IOI structure there.'}';
        
%Common to most modules: for creating a new directory and copying IOI.mat
IOImatCopyChoice           = cfg_choice;
IOImatCopyChoice.name      = 'Choose IOI copy method';
IOImatCopyChoice.tag       = 'IOImatCopyChoice';
IOImatCopyChoice.values    = {IOImatOverwrite IOImatCopy}; 
IOImatCopyChoice.val       = {IOImatOverwrite}; 
IOImatCopyChoice.help      = {'Choose whether to overwrite the IOI.mat structure'
            'or to create a new directory'
            'and copy the IOI.mat structure there'}'; 

%%%%%%%%%%%

all_ROIs         = cfg_branch;
all_ROIs.tag     = 'all_ROIs';
all_ROIs.name    = 'All ROIs';
all_ROIs.val     = {};
all_ROIs.help    = {'All ROIs will be processed'};

selected_ROIs      = cfg_entry;
selected_ROIs.tag  = 'selected_ROIs';
selected_ROIs.name = 'Enter list of ROIs';
selected_ROIs.strtype  = 'r';
selected_ROIs.num = [1 Inf];
selected_ROIs.val{1} = 1;
selected_ROIs.help = {'Enter list of ROIs to process.'};

select_ROIs         = cfg_branch;
select_ROIs.tag     = 'select_ROIs';
select_ROIs.name    = 'Select ROIs';
select_ROIs.val     = {selected_ROIs};
select_ROIs.help    = {'Choose some ROIs to be processed'};

ROI_choice        = cfg_choice;
ROI_choice.name   = 'Choose ROI selection method';
ROI_choice.tag    = 'ROI_choice';
ROI_choice.values = {all_ROIs,select_ROIs};
ROI_choice.val    = {all_ROIs};
ROI_choice.help   = {'Choose ROI selection method'}';

%%%%%%%%%%%%%%%%%%%%%%%

all_sessions         = cfg_branch;
all_sessions.tag     = 'all_sessions';
all_sessions.name    = 'All sessions';
all_sessions.val     = {};
all_sessions.help    = {'All good enough sessions will be processed'};

selected_sessions      = cfg_entry;
selected_sessions.tag  = 'selected_sessions';
selected_sessions.name = 'Enter list of sessions';
selected_sessions.strtype  = 'r';
selected_sessions.num = [1 Inf];
selected_sessions.val{1} = 1;
selected_sessions.help = {'Enter list of sessions to process.'};

select_sessions         = cfg_branch;
select_sessions.tag     = 'select_sessions';
select_sessions.name    = 'Select sessions';
select_sessions.val     = {selected_sessions};
select_sessions.help    = {'Choose some sessions to be processed'};

session_choice        = cfg_choice;
session_choice.name   = 'Choose session selection method';
session_choice.tag    = 'session_choice';
session_choice.values = {all_sessions,select_sessions};
session_choice.val    = {all_sessions};
session_choice.help   = {'Choose session selection method'}';

window_before      = cfg_entry;
window_before.tag  = 'window_before';
window_before.name = 'Window before';
window_before.strtype  = 'r';
window_before.num = [1 1];
window_before.val  = {2};
window_before.help = {'Size of window to keep prior to each stimulation onset, in seconds.'};

window_after      = cfg_entry;
window_after.tag  = 'window_after';
window_after.name = 'Window after';
window_after.strtype  = 'r';
window_after.num = [1 1];
window_after.val  = {20};
window_after.help = {'Size of window to keep after each stimulation onset, in seconds.'};

normalize_choice      = cfg_menu;
normalize_choice.tag  = 'normalize_choice';
normalize_choice.name = 'Normalization choice';
normalize_choice.labels = {'Median over window before','Time zero'};
normalize_choice.values = {1,2};
normalize_choice.val  = {1};
normalize_choice.help = {'Normalization choice. In one test,'
    'The mean standard deviation was higher by 10% or more'
    'when using time zero as the baseline, compared to taking '
    'an average (median) over the window before.'}';

include_flow      = cfg_menu;
include_flow.tag  = 'include_flow';
include_flow.name = 'Include flow';
include_flow.labels = {'Yes','No'};
include_flow.values = {1,0};
include_flow.val  = {0};
include_flow.help = {'Include flow.'}';

include_HbT      = cfg_menu;
include_HbT.tag  = 'include_HbT';
include_HbT.name = 'Include HbT';
include_HbT.labels = {'Yes','No'};
include_HbT.values = {1,0};
include_HbT.val  = {1};
include_HbT.help = {'Include HbT.'}';

include_OD      = cfg_menu;
include_OD.tag  = 'include_OD';
include_OD.name = 'Include optical intensity';
include_OD.labels = {'Yes','No'};
include_OD.values = {1,0};
include_OD.val  = {0};
include_OD.help = {'If the optical intensity images (Green, Red, Yellow) have not been deleted'
    'previously, choose whether to generate movies for these colors.'}';

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

generate_figures      = cfg_menu;
generate_figures.tag  = 'generate_figures';
generate_figures.name = 'Show figures';
generate_figures.labels = {'Yes','No'};
generate_figures.values = {1,0};
generate_figures.val  = {0};
generate_figures.help = {'Show figures. When selecting this option, the figures will stay opened after the code has completed.'}';

save_figures      = cfg_menu;
save_figures.tag  = 'save_figures';
save_figures.name = 'Save figures';
save_figures.labels = {'Yes','No'};
save_figures.values = {1,0};
save_figures.val  = {1};
save_figures.help = {'Save figures.'}';

add_error_bars      = cfg_menu;
add_error_bars.tag  = 'add_error_bars';
add_error_bars.name = 'Add error bars';
add_error_bars.labels = {'Yes','No'};
add_error_bars.values = {1,0};
add_error_bars.val  = {0};
add_error_bars.help = {'Add error bars.'}';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hpf_butter_freq         = cfg_entry; 
hpf_butter_freq.name    = 'Cutoff frequency for HPF';
hpf_butter_freq.tag     = 'hpf_butter_freq';       
hpf_butter_freq.strtype = 'r';
hpf_butter_freq.num     = [1 1];     
hpf_butter_freq.val     = {0.01};
hpf_butter_freq.help    = {'Enter cutoff frequency in Hz for Butterworth HPF.'};

hpf_butter_order         = cfg_entry; 
hpf_butter_order.name    = 'Order of Butterworth HPF';
hpf_butter_order.tag     = 'hpf_butter_order';       
hpf_butter_order.strtype = 'r';
hpf_butter_order.num     = [1 1];     
hpf_butter_order.val     = {3};
hpf_butter_order.help    = {'Enter order of Butterworth HPF (preferred value = 3).'};

hpf_butter_On         = cfg_branch;
hpf_butter_On.tag     = 'hpf_butter_On';
hpf_butter_On.name    = 'Butterworth HP filter';
hpf_butter_On.val     = {hpf_butter_freq hpf_butter_order}; 
hpf_butter_On.help    = {'Butterworth high-pass filter.'};

hpf_butter_Off         = cfg_branch;
hpf_butter_Off.tag     = 'hpf_butter_Off';
hpf_butter_Off.name    = 'HP filter off';
hpf_butter_Off.val     = {}; 
hpf_butter_Off.help    = {'High pass filter turned off.'};

hpf_butter      = cfg_choice;
hpf_butter.tag  = 'hpf_butter';
hpf_butter.name = 'Butterworth High Pass Filter';
hpf_butter.values = {hpf_butter_On hpf_butter_Off};
hpf_butter.val = {hpf_butter_On};
hpf_butter.help = {'Choose whether to include a Butterworth High Pass Filter.'
        'Parameters are: order (e.g. 3) and frequency (e.g. 0.01 Hz)'}';
    
% ---------------------------------------------------------------------
% lpf Low-pass filter
% ---------------------------------------------------------------------
fwhm1      = cfg_entry;
fwhm1.tag  = 'fwhm1';
fwhm1.name = 'FWHM in seconds';
fwhm1.val = {1};
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

% Executable Branch
stim_mean1      = cfg_exbranch;       % This is the branch that has information about how to run this module
stim_mean1.name = 'Average stimulations';             % The display name
stim_mean1.tag  = 'stim_mean1'; %Very important: tag is used when calling for execution
stim_mean1.val  = {IOImat ROImat redo1 IOImatCopyChoice session_choice ...
    ROI_choice window_after window_before normalize_choice hpf_butter ...
    lpf_choice include_flow include_HbT include_OD extract_HRF fit_3_gamma include_nlinfit ...
    generate_global generate_figures save_figures add_error_bars ...
    remove_segment_drift};    % The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs
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
