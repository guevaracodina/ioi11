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

expeditive_mode         = cfg_branch;
expeditive_mode.tag     = 'expeditive_mode';
expeditive_mode.name    = 'Expeditive mode';
expeditive_mode.val     = {};
expeditive_mode.help    = {'Expeditive mode of data processing: '
    'only the electrophysiology will be processed'}';

standard_mode         = cfg_branch;
standard_mode.tag     = 'standard_mode';
standard_mode.name    = 'Standard mode';
standard_mode.val     = {};
standard_mode.help    = {'Standard mode of data processing: all images will be processed'};

treatment_mode = cfg_choice;
treatment_mode.name   = 'Choose treatment method';
treatment_mode.tag    = 'treatment_mode';
treatment_mode.values = {standard_mode,expeditive_mode};
treatment_mode.val    = {standard_mode};
treatment_mode.help   = {'Choose treatment method: standard, or expeditive processing'}';

%Special redo
redo1      = cfg_entry;
redo1.tag  = 'force_redo';
redo1.name = 'Force processing';
redo1.strtype = 'r';
redo1.num     = [1 Inf];
redo1.val  = {0};
redo1.help = {'Enter an array of 0 (do not force redo), 1 (force redo) or 2 (partial redo)'
    'with one entry for each subject. If only one number is entered, it will '
    'be applied to all subjects'
    'If 0 is entered, and an IOI.mat exists, module will be skipped; '
    'If 1 is entered and an IOI.mat exists, it will be ignored (no previous information kept);'
    'If 2 is entered and an IOI.mat exists, it will be loaded, and some of its fields will be modified'}';

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
%output_path.val     = {fullfile('D:/Users/')};
output_path.def     = @(val)ioi_get_defaults('msioi1.output_path_select.output_path', val{:}); 
output_path.help    = {'Choose path for .ioi output files'}; 

output_path_select         = cfg_branch;
output_path_select.tag     = 'output_path_select';
output_path_select.name    = 'output_path_select';
output_path_select.val     = {output_path}; 
output_path_select.help    = {};

output_path_choice        = cfg_choice;
output_path_choice.name   = 'Choose output path method';
output_path_choice.tag    = 'output_path_choice';
output_path_choice.values = {output_path_default,output_path_select};
output_path_choice.val    = {output_path_default};
output_path_choice.help   = {'Choose output_path_choice'}';

shrinkage_choice = ioi_dfg_shrinkage_choice;
session_choice = ioi_dfg_session_choice;

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
    'raw binary file. If that is not the case, change temp_ImNum in ioi_msioi_run'
    'Applies to old data format only'}; 

stim_cutoff         = cfg_entry; 
stim_cutoff.name    = 'Cutoff to use on stim data';
stim_cutoff.tag     = 'stim_cutoff';       
stim_cutoff.strtype = 'r';
stim_cutoff.num     = [1 1];
stim_cutoff.val     = {1.5};
stim_cutoff.help    = {'Enter voltage cutoff. Sometimes it is 5 V, but 0.5 V '
    'is also possible'
    'Only applicable for new recording system.'}'; 

acq_freq         = cfg_entry; 
acq_freq.name    = 'Acquisition frequency (of all colors) in Hz';
acq_freq.tag     = 'acq_freq';       
acq_freq.strtype = 'r';
acq_freq.num     = [1 1];
acq_freq.val     = {20};
acq_freq.help    = {'Acquisition frequency'
    'Enter 20 Hz if this is the frequency of the camera. '
    'If there are 4 colors recorded, then the TR will be 0.2 s.'}'; 

color_number         = cfg_entry; 
color_number.name    = 'Number of colors acquired';
color_number.tag     = 'color_number';       
color_number.strtype = 'r';
color_number.num     = [1 1];
color_number.val     = {4};
color_number.help    = {'Number of colors acquired'}'; 

% order_colors         = cfg_entry; 
% order_colors.name    = 'Color order';
% order_colors.tag     = 'order_colors';       
% order_colors.strtype = 's';
% order_colors.num     = [1 1];
% order_colors.val     = {'RVJL'};
% order_colors.help    = {'Enter acquisition colors in French.'}'; 

save_choice        = cfg_menu;
save_choice.name   = 'Choose saving method';
save_choice.tag    = 'save_choice';
save_choice.labels = {'One file per session','One file per block','One file per image'};
save_choice.values = {1,2,3};
save_choice.val    = {2};
save_choice.help   = {'Choose saving method'
    'Applies to old data format only; new format: always one file per block'}';

memmapfileOn        = cfg_menu;
memmapfileOn.name   = 'Choose memory management method';
memmapfileOn.tag    = 'memmapfileOn';
memmapfileOn.labels = {'Keep all in memory','Use disk space for large structures'};
memmapfileOn.values = {0,1};
memmapfileOn.val    = {1};
memmapfileOn.help   = {'Select memory management method. Keeping all in memory'
    'is faster, but may require too much memory.'
    'Applies to old data format only; new format: always using disk space'}';

forceProcessingOn        = cfg_menu;
forceProcessingOn.name   = 'Force processing of bad sessions';
forceProcessingOn.tag    = 'forceProcessingOn';
forceProcessingOn.labels = {'Yes','No'};
forceProcessingOn.values = {1,0};
forceProcessingOn.val    = {0};
forceProcessingOn.help   = {'Force processing of bad sessions: attempt will be'
    'made to process sessions with inconsistent number of files; '
    'Applies to old data format only'}';

% Executable Branch
msioi1      = cfg_exbranch;       % This is the branch that has information about how to run this module
msioi1.name = 'Read Multi-Spectral IOI';             % The display name
msioi1.tag  = 'msioi1'; %Very important: tag is used when calling for execution
msioi1.val  = {top_bin_dir treatment_mode redo1 shrinkage_choice output_path_choice ...
    session_choice save_choice memmapfileOn sess_min_image_files ...
    color_number stim_cutoff acq_freq forceProcessingOn};    % The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs
msioi1.prog = @ioi_msioi_run;  % A function handle that will be called with the harvested job to run the computation
msioi1.vout = @ioi_cfg_vout_msioi; % A function handle that will be called with the harvested job to determine virtual outputs
msioi1.help = {'Module to create .nifti images from .bin images'
    'for intrinsic optical imaging'
    'This module is used for both recording formats (new and old)'
    'Note that some options work only with the old format'
    '(Saving Method, Memory Management)'
    'some only with the new format'
    'It is recommended to remove unwanted sessions in raw Rep folder '
    'before beginning processing.'}';
return

%make IOI.mat available as a dependency
function vout = ioi_cfg_vout_msioi(job)
vout = cfg_dep;                     % The dependency object
vout.sname      = 'IOI.mat';       % Displayed dependency name
vout.src_output = substruct('.','IOImat'); %{1}); %,'IOImat');
%substruct('()',{1}); % The output subscript reference. This could be any reference into the output variable created during computation
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
