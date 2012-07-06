function correlation_map1 = ioi_correlation_map_cfg
% Graphical interface configuration function for ioi_correlation_map_run
% A functional connectivity (fcIOS) map is made by correlating the seed/ROI with
% all other brain (non-masked) pixels
%_______________________________________________________________________________
% Copyright (C) 2012 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________

% Select IOI.mat
IOImat = ioi_dfg_IOImat(1);
% Force processing
redo1 = ioi_dfg_redo(0);
% IOI copy/overwrite method
IOImatCopyChoice = ioi_dfg_IOImatCopyChoice('corrMap');
% Choose ROI selection method (all/selected)
ROI_choice = ioi_dfg_ROI_choice;
% Choose session selection method (all/selected)
session_choice = ioi_dfg_session_choice;
% Colors to include (OD,HbO,HbR,HbT,Flow)
IC = ioi_dfg_include_colors(0,1,1,1,1);
% Generate / save figures
[generate_figures save_figures] = ioi_dfg_generate_figures;


% Executable Branch
correlation_map1        = cfg_exbranch;       % This is the branch that has information about how to run this module
correlation_map1.name   = 'Functional connectivity (fcIOS) map';             % The display name
correlation_map1.tag    = 'correlation_map1'; %Very important: tag is used when calling for execution
correlation_map1.val    = {IOImat redo1 IOImatCopyChoice ROI_choice session_choice IC generate_figures save_figures};    % The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs
correlation_map1.prog   = @ioi_correlation_map_run;  % A function handle that will be called with the harvested job to run the computation
correlation_map1.vout   = @ioi_cfg_vout_correlation_map; % A function handle that will be called with the harvested job to determine virtual outputs
correlation_map1.help   = {'A functional connectivity (fcIOS) map is made by correlating the seed/ROI with all other brain (non-masked) pixels'}';

return

% Make IOI.mat available as a dependency
function vout = ioi_cfg_vout_correlation_map(job)
vout = cfg_dep;                     % The dependency object
vout.sname      = 'IOI.mat';       % Displayed dependency name
vout.src_output = substruct('.','IOImat'); %{1}); %,'IOImat');
%substruct('()',{1}); % The output subscript reference. This could be any reference into the output variable created during computation
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});

% EOF
