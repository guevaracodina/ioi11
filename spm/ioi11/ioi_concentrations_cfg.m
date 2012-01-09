function conc1 = ioi_concentrations_cfg
%_______________________________________________________________________
% Copyright (C) 2010 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%______________________________________________________________________
% Read IOI Multispectral files
% ---------------------------------------------------------------------

IOImat         = cfg_files; %Select NIRS.mat for this subject 
IOImat.name    = 'Select IOI.mat'; % The displayed name
IOImat.tag     = 'IOImat';       %file names
IOImat.filter = 'mat';
IOImat.ufilter = '^IOI.mat$';    
IOImat.num     = [1 Inf];     % Number of inputs required 
IOImat.help    = {'Select IOI.mat for this subject.'}; % help text displayed

redo1      = cfg_menu;
redo1.tag  = 'force_redo';
redo1.name = 'Force processing';
redo1.labels = {'False','True'};
redo1.values = {0,1};
redo1.val  = {0};
redo1.help = {'Force redoing this processing even when it has been done already'};

pathlength1      = cfg_menu;
pathlength1.tag  = 'pathlength';
pathlength1.name = 'Pathlength Factor';
pathlength1.labels = {'Kohl','Dunn'};
pathlength1.values = {'Kohl','Dunn'};
pathlength1.def  = @(val)ioi_get_defaults('conc1.pathlength', val{:});
pathlength1.help = {'Which pathlength factor to use from literature'};

camera1      = cfg_menu;
camera1.tag  = 'camera_corr';
camera1.name = 'Camera correction';
camera1.labels = {'False','True'};
camera1.values = {0,1};
camera1.def  = @(val)ioi_get_defaults('conc1.camera_correction', val{:});
camera1.help = {'Correct for camera wavelength sensitivity'};

led1      = cfg_menu;
led1.tag  = 'leds_conv';
led1.name = 'Leds Convolution';
led1.labels = {'False','True'};
led1.values = {0,1};
led1.def  = @(val)ioi_get_defaults('conc1.leds_spectra', val{:});
led1.help = {'Do we convolve by LEDs spectra curves to get effective epsilons?'};

basehbt1      = cfg_entry;
basehbt1.tag  = 'HbT0';
basehbt1.name = 'HbT0';
basehbt1.strtype  = 'r';
basehbt1.num = [1 1];
basehbt1.def  = @(val)ioi_get_defaults('conc1.baseline_hbt', val{:});
basehbt1.help = {'Baseline total hemoglobin, in microM, default 100 microM.'};

configuration         = cfg_branch;
configuration.tag     = 'configuration';
configuration.name    = 'Configuration options';
configuration.val     = {pathlength1 camera1 led1 basehbt1};
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
NewIOIdir.val{1}  = 'Conc';
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

MemoryManagementMenu      = cfg_menu;
MemoryManagementMenu.tag  = 'MemoryManagementMenu';
MemoryManagementMenu.name = 'Memory Management';
MemoryManagementMenu.labels = {'Load all at once','Load one at a time'};
MemoryManagementMenu.values = {1,0};
MemoryManagementMenu.val  = {1};
MemoryManagementMenu.help = {'Load all images at once (faster but requires'
    'more memory, or load one image at a time, to compute concentrations.'}';

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

%%%%%%%%%%%%%%%%%%%%
% Normalization choice
%%%%%%%%%%%%%%%%%%%%%%%%%%

current_sessions         = cfg_branch;
current_sessions.tag     = 'current_sessions';
current_sessions.name    = 'Current session';
current_sessions.val     = {};
current_sessions.help    = {'Use normalization calculated previously' 
    'for this session'}';

selected_norm_session      = cfg_entry;
selected_norm_session.tag  = 'selected_norm_session';
selected_norm_session.name = 'Enter session to use for normalization';
selected_norm_session.strtype  = 'r';
selected_norm_session.num = [1 1];
selected_norm_session.val{1} = 1;
selected_norm_session.help = {'Enter session to use for normalization.'
    'The median for this session will be used for all the other sessions.'}';

select_norm_session         = cfg_branch;
select_norm_session.tag     = 'select_norm_session';
select_norm_session.name    = 'Select session median to use for normalization';
select_norm_session.val     = {selected_norm_session};
select_norm_session.help    = {'Select session median to use for normalization'};

normalization_choice        = cfg_choice;
normalization_choice.name   = 'Choose session normalization method';
normalization_choice.tag    = 'normalization_choice';
normalization_choice.values = {current_sessions,select_norm_session};
normalization_choice.val    = {current_sessions};
normalization_choice.help   = {'Choose session selection method'}';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RemoveRGY      = cfg_menu;
RemoveRGY.tag  = 'RemoveRGY';
RemoveRGY.name = 'Remove RGY nifti images';
RemoveRGY.labels = {'Yes','No'};
RemoveRGY.values = {1,0};
RemoveRGY.val  = {1};
RemoveRGY.help = {'After concentration concentrations are obtained'
    'Optical intensity images (RGY) are usually no longer required.'}';

% Executable Branch
conc1      = cfg_exbranch;       % This is the branch that has information about how to run this module
conc1.name = 'Compute Concentrations';             % The display name
conc1.tag  = 'conc1'; %Very important: tag is used when calling for execution
conc1.val  = {IOImat redo1 IOImatCopyChoice normalization_choice configuration ...
    MemoryManagementMenu session_choice RemoveRGY};    % The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs
conc1.prog = @ioi_concentrations_run;  % A function handle that will be called with the harvested job to run the computation
conc1.vout = @ioi_cfg_vout_concentrations; % A function handle that will be called with the harvested job to determine virtual outputs
conc1.help = {'Concentration computations.'};

return

%make IOI.mat available as a dependency
function vout = ioi_cfg_vout_concentrations(job)
vout = cfg_dep;                     % The dependency object
vout.sname      = 'IOI.mat';       % Displayed dependency name
vout.src_output = substruct('.','IOImat'); %{1}); %,'IOImat');
%substruct('()',{1}); % The output subscript reference. This could be any reference into the output variable created during computation
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
