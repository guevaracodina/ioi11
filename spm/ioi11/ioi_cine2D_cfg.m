function cine2D1 = ioi_cine2D_cfg
%_______________________________________________________________________
% Copyright (C) 2011 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal

IOImat = ioi_dfg_IOImat(1);
redo1 = ioi_dfg_redo(0);
IOImatCopyChoice = ioi_dfg_IOImatCopyChoice('Cine');

%%%%%%%%%%%%%%%%%%
shrinkage_choice = ioi_dfg_shrinkage_choice;

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
session_choice = ioi_dfg_session_choice;
[window_before window_after window_offset] = ioi_dfg_window(3,20,0);

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

which_onset_type = ioi_dfg_which_onset_type;
remove_stims = ioi_dfg_remove_stims; %To add
use_stims = ioi_dfg_use_stims; %To add

% low_limit         = cfg_entry; 
% low_limit.name    = 'Enter low limit as percentage of min to max';
% low_limit.tag     = 'low_limit';       
% low_limit.strtype = 'r';
% low_limit.num     = [1 1];     
% low_limit.val     = {0};
% low_limit.help    = {'Enter low limit as percentage of min to max.'}';

% high_limit         = cfg_entry; 
% high_limit.name    = 'Enter high limit as percentage of min to max';
% high_limit.tag     = 'high_limit';       
% high_limit.strtype = 'r';
% high_limit.num     = [1 1];     
% high_limit.val     = {10};
% high_limit.help    = {'Enter high limit as percentage of min to max.'
%     'For HbR, the code will invert min and max so the user does not have to worry about it.'}';

IC = ioi_dfg_include_colors(0,1,1,1,1);
hpf_butter = ioi_dfg_hpf_butter(0,0.01,3);
lpf_choice = ioi_dfg_lpf_choice(1,0.67);

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
 
show_movie      = cfg_menu;
show_movie.tag  = 'show_movie';
show_movie.name = 'Show movie';
show_movie.labels = {'Yes','No'};
show_movie.values = {1,0};
show_movie.val  = {1};
show_movie.help = {'Show movie. In either case, the movie will be saved.'}';

%*********by Cong on 12/08/28
generate_images      = cfg_menu;
generate_images.tag  = 'generate_images';
generate_images.name = 'Show images';
generate_images.labels = {'Yes','No'};
generate_images.values = {1,0};
generate_images.val  = {0};
generate_images.help = {'Show images.'}';

save_images      = cfg_menu;
save_images.tag  = 'save_images';
save_images.name = 'Save images';
save_images.labels = {'Yes','No'};
save_images.values = {1,0};
save_images.val  = {1};
save_images.help = {'Save images. '
    'two choices were selected: Yes or No'}';


interactive_mode      = cfg_menu;
interactive_mode.tag  = 'interactive_mode';
interactive_mode.name = 'Interactive mode';
interactive_mode.labels = {'Yes','No'};
interactive_mode.values = {1,0};
interactive_mode.val  = {1};
interactive_mode.help = {'Set min/max value interactively for the colorbar for each color.'}';

% save_choice        = cfg_choice;
% save_choice.name   = 'Save images';
% save_choice.tag    = 'stim_choice';
% save_choice.values = {available_onsets,manual_onsets};
% save_choice.val    = {available_onsets};
% save_choice.help   = {'Use available onsets or enter onsets.'}';


% interactive_mode      = cfg_menu;
% interactive_mode.tag  = 'interactive_mode';
% interactive_mode.name = 'Interactive mode';
% interactive_mode.labels = {'Yes','No'};
% interactive_mode.values = {1,0};
% interactive_mode.val  = {1};
% interactive_mode.help = {'Set min/max value for the colorbar for each color.'
%     'two choices was selested: Yes or NO'}';

% save_choice        = cfg_choice;
% save_choice.name   = 'Save images';
% save_choice.tag    = 'stim_choice';
% save_choice.values = {available_onsets,manual_onsets};
% save_choice.val    = {available_onsets};
% save_choice.help   = {'Use available onsets or enter onsets.'}';


%**********************end

downFact      = cfg_entry;
downFact.tag  = 'downFact';
downFact.name = 'Downsampling factor';
downFact.val = {5};
downFact.strtype = 'r';  
downFact.num     = [1 1]; 
downFact.help    = {'Downsampling factor.'}; 

%display options -- superpose activations on anatomical image
superpose_anatomical = ioi_dfg_superpose_anatomical;
superpose_ROIs = ioi_dfg_superpose_ROIs;

% Executable Branch
cine2D1      = cfg_exbranch;       % This is the branch that has information about how to run this module
cine2D1.name = '2D Cine';             % The display name
cine2D1.tag  = 'cine2D1'; %Very important: tag is used when calling for execution
cine2D1.val  = {IOImat redo1 IOImatCopyChoice session_choice shrinkage_choice downFact ...
    stim_choice window_after window_before window_offset skip_overlap ...
    normalize_choice group_onset_types which_onset_type ...
    IC interactive_mode ...
    hpf_butter lpf_choice show_movie generate_images save_images ...
    superpose_anatomical superpose_ROIs};    % The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs

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
