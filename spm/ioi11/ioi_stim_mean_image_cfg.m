function stim_mean_image1 = ioi_stim_mean_image_cfg
%_______________________________________________________________________
% Copyright (C) 2011 LIOM Laboratoire d'Imagerie Optique et Mol�culaire
%                    �cole Polytechnique de Montr�al

IOImat         = cfg_files; %Select NIRS.mat for this subject 
IOImat.name    = 'Select IOI.mat'; % The displayed name
IOImat.tag     = 'IOImat';       %file names
IOImat.filter = 'mat';
IOImat.ufilter = '^IOI.mat$';    
IOImat.num     = [1 Inf];     % Number of inputs required 
IOImat.help    = {'Select IOImat dependency if available. '
    'Otherwise, for each subject, select IOI.mat.'}'; % help text displayed

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

window_offset      = cfg_entry;
window_offset.tag  = 'window_offset';
window_offset.name = 'Window offset';
window_offset.strtype  = 'r';
window_offset.num = [1 1];
window_offset.val  = {0};
window_offset.help = {'To look back in time, include an offset in seconds. '
    'A positive number corresponds to a shift back in time. '
    'This works in combination with variables window_after and window_before:'
    'Time 0 will be negative window_offset, window_after starts at time 0, and window_before ends at time 0.'};

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

which_onset_type         = cfg_entry; 
which_onset_type.name    = 'Enter onset type(s) to use';
which_onset_type.tag     = 'which_onset_type';       
which_onset_type.strtype = 'r';
which_onset_type.num     = [1 Inf];     
which_onset_type.val     = {1};
which_onset_type.help    = {'Enter which onset type(s) (relevant if there are'
    'several onset types.'
    'Enter (a list of) ordinal number(s) indicating the desired onset type(s) in the onset type list.'}';

remove_stims      = cfg_entry;
remove_stims.tag  = 'remove_stims';
remove_stims.name = 'Enter an array of time points (in seconds) from which to exclude onsets';
remove_stims.strtype  = 'r';
remove_stims.num = [0 Inf];
remove_stims.val{1} = '';
remove_stims.help = {'Onsets occuring within 1 second of any specified' 
    'time point will be removed from the averaging.'
    'The list of onsets kept will be written to IOI.mat'}';

use_stims      = cfg_entry;
use_stims.tag  = 'use_stims';
use_stims.name = 'Enter an array of stim onset numbers to use';
use_stims.strtype  = 'r';
use_stims.num = [0 Inf];
use_stims.val{1} = '';
use_stims.help = {'Use this option to specify, for example,'
    'that the first 10 onsets should be used, by entering the array 1:10.'
    'These onset numbers will be used for each onset type.'
    'Leave array empty in order to use all available onsets, except possibly '
    'some that are excluded by other mechanisms.'}';

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

%%%%%%%%%%%%%%%%%%
shrinkage_choice = ioi_cfg_shrinkage_choice;
%%%%%%%%%%%%%%%%%%%%%%%%%
spatial_LPF = ioi_cfg_spatial_LPF;

% Executable Branch
stim_mean_image1      = cfg_exbranch;       % This is the branch that has information about how to run this module
stim_mean_image1.name = 'Average stimulations (on images)';             % The display name
stim_mean_image1.tag  = 'stim_mean_image1'; %Very important: tag is used when calling for execution
stim_mean_image1.val  = {IOImat redo1 IOImatCopyChoice shrinkage_choice ...
    spatial_LPF session_choice ...
    window_after window_before window_offset normalize_choice hpf_butter ...
    lpf_choice include_flow include_HbT include_OD which_onset_type ...
    remove_stims use_stims remove_stims_SD ...
    generate_figures save_figures ...
    remove_segment_drift };    % The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs
stim_mean_image1.prog = @ioi_stim_mean_image_run;  % A function handle that will be called with the harvested job to run the computation
stim_mean_image1.vout = @ioi_cfg_vout_stim_mean_image; % A function handle that will be called with the harvested job to determine virtual outputs
stim_mean_image1.help = {'Calculate average over stimulations over full images.'};

return

%make IOI.mat available as a dependency
function vout = ioi_cfg_vout_stim_mean_image(job)
vout = cfg_dep;                     % The dependency object
vout.sname      = 'IOI.mat';       % Displayed dependency name
vout.src_output = substruct('.','IOImat'); %{1}); %,'IOImat');
%substruct('()',{1}); % The output subscript reference. This could be any reference into the output variable created during computation
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
