function filtdown1 = ioi_fc_GLM_on_ROI_cfg
% Graphical interface configuration function for ioi_filtdown_cfg
% GLM regression of global brain signal in resting-state from ROI/seeds time
% trace in order to remove global source of variance. 
%_______________________________________________________________________________
% Copyright (C) 2012 LIOM Laboratoire d'Imagerie Optique et Mol�culaire
%                    �cole Polytechnique de Montr�al
%_______________________________________________________________________________

% Select IOI.mat
IOImat = ioi_dfg_IOImat(1);
% Force processing
redo1 = ioi_dfg_redo(0);
% IOI copy/overwrite method
IOImatCopyChoice = ioi_dfg_IOImatCopyChoice('GLM');
% Choose ROI selection method (all/selected)
ROI_choice = ioi_dfg_ROI_choice;
% Choose session selection method (all/selected)
session_choice = ioi_dfg_session_choice;
% Colors to include (OD,HbO,HbR,HbT,Flow)
IC = ioi_dfg_include_colors(0,1,1,1,1);

% Global brain signal as a regressor
regressBrainSignal          = cfg_menu;
regressBrainSignal.tag      = 'regressBrainSignal';
regressBrainSignal.name     = 'Regress global brain signal';
regressBrainSignal.labels   = {'Yes','No'};
regressBrainSignal.values   = {1,0};
regressBrainSignal.val      = {1};
regressBrainSignal.help     = {'Include global brain signal as a regressor'};

% Executable Branch
filtdown1      = cfg_exbranch;       % This is the branch that has information about how to run this module
filtdown1.name = 'GLM regression';             % The display name
filtdown1.tag  = 'fc_GLM_on_ROI1'; %Very important: tag is used when calling for execution
filtdown1.val  = {IOImat redo1 IOImatCopyChoice ROI_choice session_choice IC regressBrainSignal};    % The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs
filtdown1.prog = @ioi_fc_GLM_on_ROI_run;  % A function handle that will be called with the harvested job to run the computation
filtdown1.vout = @ioi_cfg_vout_fc_GLM_on_ROI; % A function handle that will be called with the harvested job to determine virtual outputs
filtdown1.help = {'GLM regression of global brain signal from ROI/seeds time trace in order to remove global source of variance'};

return

% Make IOI.mat available as a dependency
function vout = ioi_cfg_vout_fc_GLM_on_ROI(job)
vout = cfg_dep;                     % The dependency object
vout.sname      = 'IOI.mat';       % Displayed dependency name
vout.src_output = substruct('.','IOImat'); %{1}); %,'IOImat');
%substruct('()',{1}); % The output subscript reference. This could be any reference into the output variable created during computation
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
