function cine2D1 = ioi_cine2D_cfg
%_______________________________________________________________________
% Copyright (C) 2011 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal

IOImat = ioi_cfg_IOImat(1);
redo1 = ioi_dfg_redo(0);
IOImatCopyChoice = ioi_cfg_IOImatCopyChoice('Cine');

%%%%%%%%%%%%%%%%%%
shrinkage_choice = ioi_cfg_shrinkage_choice;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

onset_time      = cfg_entry;
onset_time.tag  = 'onset_time';
onset_time.name = 'Onset times in seconds';
onset_time.val = {1};
onset_time.strtype = 'r';  
onset_time.num     = [1 Inf]; 
onset_time.help    = {'Onset times in seconds.'
    'Note that the same values of window_after, _before and _offset will be used for each onset time specified. '}'; 

manual_onsets         = cfg_branch;
manual_onsets.tag     = 'manual_onsets';
manual_onsets.name    = 'Manual onsets';
manual_onsets.val     = {onset_time};
manual_onsets.help    = {'Define onset times here '
    'Note that the variables window_before, _after and _offset'
    'will apply. However, the variables group_onsets and which_onsets will no longer be used.'}';

available_onsets         = cfg_branch;
available_onsets.tag     = 'available_onsets';
available_onsets.name    = 'Available onsets';
available_onsets.val     = {};
available_onsets.help    = {'Use available onsets, as created earlier either in msioi or in create_onsets'};

stim_choice        = cfg_choice;
stim_choice.name   = 'Onset choice method';
stim_choice.tag    = 'stim_choice';
stim_choice.values = {available_onsets,manual_onsets};
stim_choice.val    = {available_onsets};
stim_choice.help   = {'Use available onsets or enter onsets.'}';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
session_choice = ioi_cfg_session_choice;
[window_before window_after window_offset] = ioi_cfg_window(3,20);

normalize_choice      = cfg_menu;
normalize_choice.tag  = 'normalize_choice';
normalize_choice.name = 'Normalization choice';
normalize_choice.labels = {'Median over window before','Time zero','Mean','No normalization'};
normalize_choice.values = {1,2,3,4};
normalize_choice.val  = {1};
normalize_choice.help = {'Normalization choice. In one test,'
    'The mean standard deviation was higher by 10% or more'
    'when using time zero as the baseline, compared to taking '
    'an average (median) over the window before.'}';

which_onset_type = ioi_cfg_which_onset_type;
remove_stims = ioi_cfg_remove_stims; %To add
use_stims = ioi_cfg_use_stims; %To add

low_limit         = cfg_entry; 
low_limit.name    = 'Enter low limit as percentage of min to max';
low_limit.tag     = 'low_limit';       
low_limit.strtype = 'r';
low_limit.num     = [1 1];     
low_limit.val     = {0};
low_limit.help    = {'Enter low limit as percentage of min to max.'}';

high_limit         = cfg_entry; 
high_limit.name    = 'Enter high limit as percentage of min to max';
high_limit.tag     = 'high_limit';       
high_limit.strtype = 'r';
high_limit.num     = [1 1];     
high_limit.val     = {10};
high_limit.help    = {'Enter high limit as percentage of min to max.'
    'For HbR, the code will invert min and max so the user does not have to worry about it.'}';

IC = ioi_cfg_include_colors(0,1,1,1,1);

skip_overlap      = cfg_menu;
skip_overlap.tag  = 'skip_overlap';
skip_overlap.name = 'Skip overlap';
skip_overlap.labels = {'Yes','No'};
skip_overlap.values = {1,0};
skip_overlap.val  = {1};
skip_overlap.help = {'For event-related onsets, such as spike. This option removes spikes'
    'that are followed by another spike at an interval shorter than the sum of '
    'window before + window after, as specified by the user'}';

group_onset_types      = cfg_menu;
group_onset_types.tag  = 'group_onset_types';
group_onset_types.name = 'Group onset types into the same type';
group_onset_types.labels = {'Yes','No'};
group_onset_types.values = {1,0};
group_onset_types.val  = {1};
group_onset_types.help = {'If there are several types of onsets, '
    'they will be grouped into the same type.'
    'This option overrides the which_onset_type to use.'}';

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
lpf_choice.help = {'Choose whether to include a Low Pass Filter.'}';
 
show_movie      = cfg_menu;
show_movie.tag  = 'show_movie';
show_movie.name = 'Show movie';
show_movie.labels = {'Yes','No'};
show_movie.values = {1,0};
show_movie.val  = {1};
show_movie.help = {'Show movie. In either case, the movie will be saved.'}';

% Executable Branch
cine2D1      = cfg_exbranch;       % This is the branch that has information about how to run this module
cine2D1.name = '2D Cine';             % The display name
cine2D1.tag  = 'cine2D1'; %Very important: tag is used when calling for execution
cine2D1.val  = {IOImat redo1 IOImatCopyChoice session_choice shrinkage_choice ...
    stim_choice window_after window_before window_offset skip_overlap ...
    normalize_choice group_onset_types which_onset_type ...
    high_limit low_limit IC ...
    lpf_choice show_movie };    % The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs
cine2D1.prog = @ioi_cine2D_run;  % A function handle that will be called with the harvested job to run the computation
cine2D1.vout = @ioi_cfg_vout_cine2D; % A function handle that will be called with the harvested job to determine virtual outputs
cine2D1.help = {'Generate a 2D movie'
    'Calculate average over stimulations or other previously specified onsets.'
    'This module also plays the role of a viewer for previously generated movies.'}';

return

%make IOI.mat available as a dependency
function vout = ioi_cfg_vout_cine2D(job)
vout = cfg_dep;                     % The dependency object
vout.sname      = 'IOI.mat';       % Displayed dependency name
vout.src_output = substruct('.','IOImat'); %{1}); %,'IOImat');
%substruct('()',{1}); % The output subscript reference. This could be any reference into the output variable created during computation
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
