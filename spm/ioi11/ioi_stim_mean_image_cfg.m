function stim_mean_image1 = ioi_stim_mean_image_cfg
%_______________________________________________________________________
% Copyright (C) 2011 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal

IOImat = ioi_dfg_IOImat(1);
redo1 = ioi_dfg_redo(0);
IOImatCopyChoice = ioi_dfg_IOImatCopyChoice('IMean');

%%%%%%%%%%%%%%%%%%%%%%%
session_choice = ioi_dfg_session_choice;
[window_before window_after window_offset] = ioi_dfg_window(3,20,0);
window_start_delay = ioi_dfg_window_start_delay(0);
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

which_onset_type = ioi_dfg_which_onset_type;
remove_stims = ioi_dfg_remove_stims;
use_stims = ioi_dfg_use_stims;
[generate_figures save_figures] = ioi_dfg_generate_figures;

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

interactive_mode      = cfg_menu;
interactive_mode.tag  = 'interactive_mode';
interactive_mode.name = 'Interactive mode';
interactive_mode.labels = {'Yes','No'};
interactive_mode.values = {1,0};
interactive_mode.val  = {0};
interactive_mode.help = {'Interactive mode - let''s user choose stims.'
    'After results with one set of stims are displayed,'
    'The user can try a different set of stims, move ahead'
    'to the next image, or return to the automatic mode.'}';

%%%%%%%%%%%%%%%%%%
shrinkage_choice = ioi_dfg_shrinkage_choice;
%%%%%%%%%%%%%%%%%%%%%%%%%
spatial_LPF = ioi_dfg_spatial_LPF;
%%%%%%%%%%%%%%%%%%%%%%
IC = ioi_dfg_include_colors(0,1,1,1,1);

% Executable Branch
stim_mean_image1      = cfg_exbranch;       % This is the branch that has information about how to run this module
stim_mean_image1.name = 'Average stimulations (on images)';             % The display name
stim_mean_image1.tag  = 'stim_mean_image1'; %Very important: tag is used when calling for execution
stim_mean_image1.val  = {IOImat redo1 IOImatCopyChoice interactive_mode shrinkage_choice ...
    spatial_LPF session_choice window_start_delay ...
    window_after window_before window_offset normalize_choice hpf_butter ...
    lpf_choice IC ...
    which_onset_type ...
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
