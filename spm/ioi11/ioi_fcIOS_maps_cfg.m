function fcIOS_maps         = ioi_fcIOS_maps_cfg
% Graphical interface configuration function for ioi_fcIOS_maps_run
%_______________________________________________________________________________
% Copyright (C) 2012 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________

% Select IOI.mat
IOImat                      = ioi_dfg_IOImat(2);
% IOI copy/overwrite method
IOImatCopyChoice            = ioi_dfg_IOImatCopyChoice('fcIOS_maps');
% Colors to include (OD,HbO,HbR,HbT,Flow,CMRO2)
IC                          = ioi_dfg_include_colors(0,1,1,1,1);
% Choose session selection method (all/selected)
session_choice              = ioi_dfg_session_choice;

% Select directory to save global results
parent_results_dir          = cfg_files;
parent_results_dir.tag      = 'parent_results_dir';
parent_results_dir.name     = 'Top directory to save group results';
parent_results_dir.filter   = 'dir';
parent_results_dir.num      = [1 1];
parent_results_dir.help     = {'Select the directory where figures will be saved.'}';

% Executable Branch
fcIOS_maps                  = cfg_exbranch; % This is the branch that has information about how to run this module
fcIOS_maps.name             = 'fcIOS maps'; % The display name
fcIOS_maps.tag              = 'fcIOS_maps'; %Very important: tag is used when calling for execution
fcIOS_maps.val              = {IOImat IOImatCopyChoice IC session_choice parent_results_dir };    % The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs
fcIOS_maps.prog             = @ioi_fcIOS_maps_run; % A function handle that will be called with the harvested job to run the computation
fcIOS_maps.vout             = @ioi_cfg_vout_fcIOS_maps; % A function handle that will be called with the harvested job to determine virtual outputs
fcIOS_maps.help             = {'Plots multiple correlation maps at the same scale'};

return

% Make IOI.mat available as a dependency
function vout = ioi_cfg_vout_fcIOS_maps(job)
vout                        = cfg_dep;                  % The dependency object
vout.sname                  = 'IOI.mat';                % Displayed dependency name
vout.src_output             = substruct('.','IOImat');  %{1}); %,'IOImat');
vout.tgt_spec               = cfg_findspec({{'filter','mat','strtype','e'}});

% EOF
