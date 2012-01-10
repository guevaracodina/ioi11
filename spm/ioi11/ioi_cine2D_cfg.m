function cine2D1 = ioi_cine2D_cfg
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
NewIOIdir.val{1}    = 'Cine';
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%

shrink_x      = cfg_entry;
shrink_x.tag  = 'shrink_x';
shrink_x.name = 'Shrink factor for x dimension';
shrink_x.strtype  = 'i';
shrink_x.num = [1 1];
shrink_x.val  = {2};
shrink_x.help = {'Data reduction factor in x.'};

shrink_y      = cfg_entry;
shrink_y.tag  = 'shrink_y';
shrink_y.name = 'Shrink factor for y dimension';
shrink_y.strtype  = 'i';
shrink_y.num = [1 1];
shrink_y.val = {2};
shrink_y.help = {'Data reduction factor in y.'};

force_shrink_recompute      = cfg_menu;
force_shrink_recompute.tag  = 'force_shrink_recompute';
force_shrink_recompute.name = 'Force shrink recompute';
force_shrink_recompute.labels = {'Yes','No'};
force_shrink_recompute.values = {1,0};
force_shrink_recompute.val  = {0};
force_shrink_recompute.help = {'This option is used when images of one shrunk size are present'
    'But a different size is desired. This forces generating new shrunk images. The old ones'
    'will be kept too, but will not be available from the new IOI.mat. Therefore, one should'
    'usually select the option to place the new IOI.mat in a new folder, so that the old IOI.mat'
    'can still be used to access the first set of shrunk images.'}';


configuration_shrink         = cfg_branch;
configuration_shrink.tag     = 'configuration_shrink';
configuration_shrink.name    = 'Configuration shrinkage';
configuration_shrink.val     = {shrink_x shrink_y force_shrink_recompute};
configuration_shrink.help    = {'Select values.'};

no_shrinkage         = cfg_branch;
no_shrinkage.tag     = 'no_shrinkage';
no_shrinkage.name    = 'No shrinkage';
no_shrinkage.val     = {};
no_shrinkage.help    = {};

shrinkage_choice        = cfg_choice;
shrinkage_choice.name   = 'Choose shrinkage method';
shrinkage_choice.tag    = 'shrinkage_choice';
shrinkage_choice.values = {no_shrinkage,configuration_shrink};
shrinkage_choice.val    = {configuration_shrink};
shrinkage_choice.help   = {'Choose whether to shrink the data. Images will then be stored. And could be reused later.'}';

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
normalize_choice.labels = {'Median over window before','Time zero','Mean','No normalization'};
normalize_choice.values = {1,2,3,4};
normalize_choice.val  = {1};
normalize_choice.help = {'Normalization choice. In one test,'
    'The mean standard deviation was higher by 10% or more'
    'when using time zero as the baseline, compared to taking '
    'an average (median) over the window before.'}';

which_onset_type         = cfg_entry; 
which_onset_type.name    = 'Enter onset type(s) to use';
which_onset_type.tag     = 'which_onset_type';       
which_onset_type.strtype = 'r';
which_onset_type.num     = [1 Inf];     
which_onset_type.val     = {1};
which_onset_type.help    = {'Enter which onset type(s) (relevant if there are'
    'several onset types.'
    'Enter (a list of) ordinal number(s) indicating the desired onset type(s) in the onset type list.'}';

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

include_HbO      = cfg_menu;
include_HbO.tag  = 'include_HbO';
include_HbO.name = 'Include HbO';
include_HbO.labels = {'Yes','No'};
include_HbO.values = {1,0};
include_HbO.val  = {1};
include_HbO.help = {'Include HbO.'}';

include_HbR      = cfg_menu;
include_HbR.tag  = 'include_HbR';
include_HbR.name = 'Include HbR';
include_HbR.labels = {'Yes','No'};
include_HbR.values = {1,0};
include_HbR.val  = {1};
include_HbR.help = {'Include HbR.'}';

include_OD      = cfg_menu;
include_OD.tag  = 'include_OD';
include_OD.name = 'Include optical intensity';
include_OD.labels = {'Yes','No'};
include_OD.values = {1,0};
include_OD.val  = {0};
include_OD.help = {'If the optical intensity images (Green, Red, Yellow) have not been deleted'
    'previously, choose whether to generate movies for these colors.'}';

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
    stim_choice window_after window_before window_offset skip_overlap normalize_choice group_onset_types which_onset_type ...
    high_limit low_limit include_flow include_OD include_HbT include_HbR include_HbO ...
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
