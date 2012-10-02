function elinfo1    = ioi_get_elinfo_data_cfg
% Graphical interface configuration function for ioi_update_elinfo_run
%_______________________________________________________________________________
% Copyright (C) 2010 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________

% Select IOI.mat
IOImat              = ioi_dfg_IOImat(1);
% Force processing
redo1               = ioi_dfg_redo(0);
% IOI copy/overwrite method
IOImatCopyChoice    = ioi_dfg_IOImatCopyChoice('elinfoUpdate');
% Choose session selection method (all/selected)
session_choice      = ioi_dfg_session_choice;
% Generate / save figures
[generate_figures ...
    save_figures]   = ioi_dfg_generate_figures;

% Executable Branch
elinfo1             = cfg_exbranch; % This is the branch that has information about how to run this module
elinfo1.name        = 'Update elinfo';             % The display name
elinfo1.tag         = 'elinfo1'; %Very important: tag is used when calling for execution
elinfo1.val         = {IOImat redo1 IOImatCopyChoice session_choice generate_figures save_figures};    % The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs
elinfo1.prog        = @ioi_get_elinfo_data_run;  % A function handle that will be called with the harvested job to run the computation
elinfo1.vout        = @ioi_cfg_vout_elinfo; % A function handle that will be called with the harvested job to determine virtual outputs
elinfo1.help        = {'Update elinfo data in IOI matrix'};

return

% Make IOI.mat available as a dependency
function vout = ioi_cfg_vout_elinfo(job)
vout                = cfg_dep; % The dependency object
vout.sname          = 'IOI.mat'; % Displayed dependency name
vout.src_output     = substruct('.','IOImat'); %{1}); %,'IOImat');
%substruct('()',{1}); % The output subscript reference. This could be any reference into the output variable created during computation
vout.tgt_spec       = cfg_findspec({{'filter','mat','strtype','e'}});

% EOF

