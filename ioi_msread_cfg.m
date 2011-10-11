function msread1 = ioi_msread_cfg
%_______________________________________________________________________
% Copyright (C) 2011 LIOM Laboratoire d'Imagerie Optique et Moléculaire
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

redo1      = cfg_entry;
redo1.tag  = 'force_redo';
redo1.name = 'Force processing';
redo1.strtype = 'r';
redo1.num     = [1 Inf];
redo1.val  = {0};
redo1.help = {'Enter an array of 0 (do not force redo) and 1 (force redo) '
    'with one entry for each subject. If only one number is entered, it will '
    'be applied to all subjects'}';

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
output_path.val     = {fullfile('D:/Users/')};
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

% Executable Branch
msread1      = cfg_exbranch;       % This is the branch that has information about how to run this module
msread1.name = 'Read Multi-Spectral IOI, new, interlaced format';             % The display name
msread1.tag  = 'msread1'; %Very important: tag is used when calling for execution
msread1.val  = {top_bin_dir redo1 output_path_choice};    % The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs
msread1.prog = @ioi_msread_run;  % A function handle that will be called with the harvested job to run the computation
msread1.vout = @ioi_cfg_vout_msread; % A function handle that will be called with the harvested job to determine virtual outputs
msread1.help = {'This is for the new format - interlaced colors and TDMS'
    'Module to create .nifti images from .bin images'
    'for intrinsic optical imaging'}';
return

%make IOI.mat available as a dependency
function vout = ioi_cfg_vout_msread(job)
vout = cfg_dep;                     % The dependency object
vout.sname      = 'IOI.mat';       % Displayed dependency name
vout.src_output = substruct('.','IOImat'); %{1}); %,'IOImat');
%substruct('()',{1}); % The output subscript reference. This could be any reference into the output variable created during computation
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
