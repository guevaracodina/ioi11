function flow1 = ioi_flow_cfg
%_______________________________________________________________________
% Copyright (C) 2011 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%______________________________________________________________________
% Computes blood flow from speckle contrast
% ---------------------------------------------------------------------

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
redo1.help = {'Force redoing this processing even when it has been done already'};

inttime1      = cfg_entry;
inttime1.tag  = 'integ_time';
inttime1.name = 'Integration Time';
inttime1.strtype  = 'r';
inttime1.num = [1 1];
inttime1.def  = @(val)ioi_get_defaults('flow1.T', val{:});
inttime1.help = {'Camera integration time.'};

windowsize1      = cfg_entry;
windowsize1.tag  = 'window_size';
windowsize1.name = 'Speckle Window Size';
windowsize1.strtype  = 'i';
windowsize1.num = [1 1];
windowsize1.def  = @(val)ioi_get_defaults('flow1.window_size', val{:});
windowsize1.help = {'Window size for speckle contrast computation.'};

configuration         = cfg_branch;
configuration.tag     = 'configuration';
configuration.name    = 'Configuration options';
configuration.val     = {inttime1 windowsize1};
configuration.help    = {'Select values.'};

IOImatOverwrite         = cfg_branch;
IOImatOverwrite.tag     = 'IOImatOverwrite';
IOImatOverwrite.name    = 'Overwrite IOI.mat structure'; 
IOImatOverwrite.help    = {'Will not copy IOI structure.'
            'This will write over the previous NIRS.mat'}';

NewIOIdir         = cfg_entry;
NewIOIdir.name    = 'New directory for IOI.mat';
NewIOIdir.tag     = 'NewIOIdir';       
NewIOIdir.strtype = 's';
NewIOIdir.val{1}    = 'Flow';
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

RemoveLC      = cfg_menu;
RemoveLC.tag  = 'RemoveLC';
RemoveLC.name = 'Remove Laser/Contrast nifti images';
RemoveLC.labels = {'Yes','No'};
RemoveLC.values = {1,0};
RemoveLC.val  = {1};
RemoveLC.help = {'After flow images are obtained'
    'Laser and contrast nifti images (LC) are usually no longer required.'}';

% Executable Branch
flow1      = cfg_exbranch;       % This is the branch that has information about how to run this module
flow1.name = 'Compute Flow';             % The display name
flow1.tag  = 'flow1'; %Very important: tag is used when calling for execution
flow1.val  = {IOImat redo1 IOImatCopyChoice configuration ...
    session_choice RemoveLC};    % The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs
flow1.prog = @ioi_flow_run;  % A function handle that will be called with the harvested job to run the computation
flow1.vout = @ioi_cfg_vout_flow; % A function handle that will be called with the harvested job to determine virtual outputs
flow1.help = {'Flow computations.'};

return

%make IOI.mat available as a dependency
function vout = ioi_cfg_vout_flow(job)
vout = cfg_dep;                     % The dependency object
vout.sname      = 'IOI.mat';       % Displayed dependency name
vout.src_output = substruct('.','IOImat'); %{1}); %,'IOImat');
%substruct('()',{1}); % The output subscript reference. This could be any reference into the output variable created during computation
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
