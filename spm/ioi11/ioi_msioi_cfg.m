function msioi1 = ioi_msioi_cfg
%_______________________________________________________________________
% Copyright (C) 2010 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%______________________________________________________________________
% Read IOI Multispectral files
% ---------------------------------------------------------------------

%Select top level directory with .bin files
top_bin_dir         = cfg_files;
top_bin_dir.tag     = 'top_bin_dir';
top_bin_dir.name    = 'Subject top level bin directory';
top_bin_dir.filter = 'dir'; 
top_bin_dir.num     = [1 Inf];
top_bin_dir.help    = {'For each subject, select the top level directory'
    'containing folders of .bin image files and folders of recording information'}';

redo1      = cfg_menu;
redo1.tag  = 'force_redo';
redo1.name = 'Force processing';
redo1.labels = {'False','True','Partial'};
redo1.values = {0,1,2};
redo1.def  = @(val)ioi_get_defaults('msioi1.force_redo', val{:});
redo1.help = {'If false and an IOI.mat exists, module will be skipped; '
    'If true and an IOI.mat exists, it will be ignored (no previous information kept);'
    'If partial and an IOI.mat exists, it will be loaded, and some of its fields will be modified'}';

subj         = cfg_branch;
subj.tag     = 'subj';
subj.name    = 'Subject';
subj.val     = {top_bin_dir redo1}; 
subj.help    = {};

%path structure
output_path_default         = cfg_branch;
output_path_default.name     = 'Default output path';
output_path_default.tag    = 'output_path_default';
output_path_default.val     = {}; 
output_path_default.help    = {'Default output path: '
    'Directory will be created one level above top-level selected directory.'
    'For example, if path to subject raw data is ..\MyGroupExpt\Subj1'
    'then a folder \Nifti\Subj1 will be created for the output'}';

output_path         = cfg_entry; %path
output_path.name    = 'path for .ioi output files';
output_path.tag     = 'output_path';       
output_path.strtype = 's';
output_path.num     = [1 Inf];     
output_path.def     = @(val)ioi_get_defaults('msioi1.output_path_select.output_path', val{:}); 
output_path.help    = {'Choose path for .ioi output files'}; 

output_path_select         = cfg_branch;
output_path_select.tag     = 'output_path_select';
output_path_select.name    = 'output_path_select';
output_path_select.val     = {output_path}; %{input1 input2 input3 input4 input5 input6 redo1};
output_path_select.help    = {};

output_path_choice        = cfg_choice;
output_path_choice.name   = 'Choose output path method';
output_path_choice.tag    = 'output_path_choice';
output_path_choice.values = {output_path_default,output_path_select};
output_path_choice.val    = {output_path_default};
output_path_choice.help   = {'Choose output_path_choice'}';

shrink_x      = cfg_entry;
shrink_x.tag  = 'shrink_x';
shrink_x.name = 'Shrink factor for x dimension';
shrink_x.strtype  = 'i';
shrink_x.num = [1 1];
shrink_x.def  = @(val)ioi_get_defaults('msioi1.shrink_x', val{:});
shrink_x.help = {'Data reduction factor in x.'};

shrink_y      = cfg_entry;
shrink_y.tag  = 'shrink_y';
shrink_y.name = 'Shrink factor for y dimension';
shrink_y.strtype  = 'i';
shrink_y.num = [1 1];
shrink_y.def  = @(val)ioi_get_defaults('msioi1.shrink_y', val{:});
shrink_y.help = {'Data reduction factor in y.'};

configuration_shrink         = cfg_branch;
configuration_shrink.tag     = 'configuration_shrink';
configuration_shrink.name    = 'Configuration shrinkage';
configuration_shrink.val     = {shrink_x shrink_y};
configuration_shrink.help    = {'Select values.'};

no_shrinkage         = cfg_branch;
no_shrinkage.tag     = 'no_shrinkage';
no_shrinkage.name    = 'No shrinkage';
no_shrinkage.val     = {};
no_shrinkage.help    = {};

configuration_choice        = cfg_choice;
configuration_choice.name   = 'Choose configuration method';
configuration_choice.tag    = 'configuration_choice';
configuration_choice.values = {no_shrinkage,configuration_shrink};
configuration_choice.val    = {no_shrinkage};
configuration_choice.help   = {'Choose output_path_choice'}';

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

sess_min_image_files         = cfg_entry; 
sess_min_image_files.name    = 'Minimum length of each session in seconds';
sess_min_image_files.tag     = 'sess_min_image_files';       
sess_min_image_files.strtype = 'r';
sess_min_image_files.num     = [1 1];
sess_min_image_files.val     = {60};
sess_min_image_files.help    = {'Minimum length of each session in seconds'
    'to declare a session valid'
    'Otherwise, session will be ignored'
    'and will be absent from the list of processed sessions'
    'Note that for this, an assumption is made that data is acquired at 5 Hz'
    'If that is not the case, change temp_TR in ioi_msioi_run'
    'Another assumption is made that there are approximately 80 images per'
    'raw binary file. If that is not the case, change temp_ImNum in ioi_msioi_run'};   

% save_choice        = cfg_choice;
% save_choice.name   = 'Choose saving method';
% save_choice.tag    = 'save_choice';
% save_choice.values = {one_file_per_session,one_file_per_block,one_file_per_image};
% save_choice.val    = {one_file_per_session};
% save_choice.help   = {'Choose saving method'}';

save_choice        = cfg_menu;
save_choice.name   = 'Choose saving method';
save_choice.tag    = 'save_choice';
save_choice.labels = {'One file per session','One file per block','One file per image'};
save_choice.values = {1,2,3};
save_choice.val    = {2};
save_choice.help   = {'Choose saving method'}';

memmapfileOn        = cfg_menu;
memmapfileOn.name   = 'Choose memory management method';
memmapfileOn.tag    = 'memmapfileOn';
memmapfileOn.labels = {'Keep all in memory','Use disk space for large structures'};
memmapfileOn.values = {0,1};
memmapfileOn.val    = {1};
memmapfileOn.help   = {'Select memory management method. Keeping all in memory'
    'is faster, but may require too much memory.'}';

forceProcessingOn        = cfg_menu;
forceProcessingOn.name   = 'Force processing of bad sessions';
forceProcessingOn.tag    = 'forceProcessingOn';
forceProcessingOn.labels = {'Yes','No'};
forceProcessingOn.values = {1,0};
forceProcessingOn.val    = {0};
forceProcessingOn.help   = {'Force processing of bad sessions: attempt will be'
    'made to process sessions with inconsistent number of files; '}';

% Executable Branch
msioi1      = cfg_exbranch;       % This is the branch that has information about how to run this module
msioi1.name = 'Read Multi-Spectral IOI';             % The display name
msioi1.tag  = 'msioi1'; %Very important: tag is used when calling for execution
msioi1.val  = {subj configuration_choice output_path_choice ...
    session_choice save_choice memmapfileOn sess_min_image_files forceProcessingOn};    % The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs
msioi1.prog = @ioi_msioi_run;  % A function handle that will be called with the harvested job to run the computation
msioi1.vout = @ioi_cfg_vout_msioi; % A function handle that will be called with the harvested job to determine virtual outputs
msioi1.help = {'Module to create .nifti images from .bin images'
    'for intrinsic optical imaging'}';
return

%make IOI.mat available as a dependency
function vout = ioi_cfg_vout_msioi(job)
vout = cfg_dep;                     % The dependency object
vout.sname      = 'IOI.mat';       % Displayed dependency name
vout.src_output = substruct('.','IOImat'); %{1}); %,'IOImat');
%substruct('()',{1}); % The output subscript reference. This could be any reference into the output variable created during computation
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
