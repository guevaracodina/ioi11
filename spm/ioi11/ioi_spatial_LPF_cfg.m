function spatial_LPF1 = ioi_spatial_LPF_cfg
% Graphical interface configuration function for ioi_spatial_LPF_run
% Low-pass filtering of 2-D images with a rotationally symmetric gaussian kernel
%_______________________________________________________________________________
% Copyright (C) 2012 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________

% Select IOI.mat
IOImat = ioi_dfg_IOImat(1);
% Force processing
redo1 = ioi_dfg_redo(0);
% IOI copy/overwrite method
IOImatCopyChoice = ioi_dfg_IOImatCopyChoice('LPF');
% Choose session selection method (all/selected)
session_choice = ioi_dfg_session_choice;
% Colors to include (OD,HbO,HbR,HbT,Flow)
IC = ioi_dfg_include_colors(0,1,1,1,1);
% Spatial low-pass filter options
spatial_LPF_options = ioi_dfg_spatial_LPF;

% Executable Branch
spatial_LPF1      = cfg_exbranch;       % This is the branch that has information about how to run this module
spatial_LPF1.name = 'Spatial low-pass filtering';             % The display name
spatial_LPF1.tag  = 'spatial_LPF1'; %Very important: tag is used when calling for execution
spatial_LPF1.val  = {IOImat redo1 IOImatCopyChoice session_choice IC spatial_LPF_options};    % The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs
spatial_LPF1.prog = @ioi_spatial_LPF_run;  % A function handle that will be called with the harvested job to run the computation
spatial_LPF1.vout = @ioi_cfg_vout_spatial_LPF; % A function handle that will be called with the harvested job to determine virtual outputs
spatial_LPF1.help = {'Low-pass filtering of 2-D images using a rotationally symmetric gaussian kernel.'
    'Usually done after computing concentrations/flow.'}';

return

% Make IOI.mat available as a dependency
function vout = ioi_cfg_vout_spatial_LPF(job)
vout = cfg_dep;                     % The dependency object
vout.sname      = 'IOI.mat';       % Displayed dependency name
vout.src_output = substruct('.','IOImat'); %{1}); %,'IOImat');
%substruct('()',{1}); % The output subscript reference. This could be any reference into the output variable created during computation
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});

% EOF
