function elinfostats1    = ioi_get_elinfo_data_cfg
% Graphical interface configuration function for ioi_get_elinfo_data_run
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
elinfostats1        = cfg_exbranch; % This is the branch that has information about how to run this module
elinfostats1.name   = 'Get ECG/Temp. stats';             % The display name
elinfostats1.tag    = 'elinfostats1'; %Very important: tag is used when calling for execution
elinfostats1.val    = {IOImat redo1 IOImatCopyChoice session_choice generate_figures save_figures};    % The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs
elinfostats1.prog   = @ioi_get_elinfo_data_run;  % A function handle that will be called with the harvested job to run the computation
elinfostats1.vout   = @ioi_cfg_vout_elinfostats; % A function handle that will be called with the harvested job to determine virtual outputs
elinfostats1.help   = {'Get stats on ecg (bpm) and temperature.'};

return

% Make IOI.mat available as a dependency
function vout = ioi_cfg_vout_elinfostats(job)
vout                = cfg_dep; % The dependency object
vout.sname          = 'IOI.mat'; % Displayed dependency name
vout.src_output     = substruct('.','IOImat'); %{1}); %,'IOImat');
%substruct('()',{1}); % The output subscript reference. This could be any reference into the output variable created during computation
vout.tgt_spec       = cfg_findspec({{'filter','mat','strtype','e'}});

% EOF

