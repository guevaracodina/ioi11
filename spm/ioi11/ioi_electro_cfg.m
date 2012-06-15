function electro1 = ioi_electro_cfg
%_______________________________________________________________________
% Copyright (C) 2011 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%______________________________________________________________________

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

IOImatCopyChoice = ioi_cfg_IOImatCopyChoice('Flow');

session_choice = ioi_cfg_session_choice;

% Executable Branch
electro1      = cfg_exbranch;       % This is the branch that has information about how to run this module
electro1.name = 'Electrophysiology extraction';             % The display name
electro1.tag  = 'electro1'; %Very important: tag is used when calling for execution
electro1.val  = {IOImat redo1 IOImatCopyChoice ...
    session_choice};    % The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs
electro1.prog = @ioi_electro_run;  % A function handle that will be called with the harvested job to run the computation
electro1.vout = @ioi_cfg_vout_electro; % A function handle that will be called with the harvested job to determine virtual outputs
electro1.help = {'Electrophysiology extraction.'};

return

%make IOI.mat available as a dependency
function vout = ioi_cfg_vout_electro(job)
vout = cfg_dep;                     % The dependency object
vout.sname      = 'IOI.mat';       % Displayed dependency name
vout.src_output = substruct('.','IOImat'); %{1}); %,'IOImat');
%substruct('()',{1}); % The output subscript reference. This could be any reference into the output variable created during computation
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
