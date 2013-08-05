function mask1 = ioi_brainmask_cfg
% Batch job configuration system for ioi_brainmask_run
% Manual segmentation of the brain to provide a mask for fcOIS analysis
%_______________________________________________________________________________
% Copyright (C) 2012 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________

% Select IOI.mat
IOImat                      = ioi_dfg_IOImat(1);
% Force processing
redo1                       = ioi_dfg_redo(0);
% IOI copy/overwrite method
IOImatCopyChoice            = ioi_dfg_IOImatCopyChoice('brainMask');

% Choose whether interface should ask whether to use a previously saved list of ROI
SelectPreviousMask          = cfg_menu;
SelectPreviousMask.tag      = 'SelectPreviousMask';
SelectPreviousMask.name     = 'Ask for a possible previous brain mask';
SelectPreviousMask.labels   = {'No','Yes'};
SelectPreviousMask.values   = {false,  true};
SelectPreviousMask.val      = {false};
SelectPreviousMask.help     = {'If this option is selected, then before selection of brain mask,'
    'The interface will ask the user if they want to use a brain mask from elsewhere.'
    'This brain mask can then be selected by locating the IOI.mat structure that contains the information.'
    'Useful when re-running the pipeline'}';

% Executable Branch
mask1                       = cfg_exbranch; % This is the branch that has information about how to run this module
mask1.name                  = 'Create Brain Mask'; % The display name
mask1.tag                   = 'mask1'; %Very important: tag is used when calling for execution
mask1.val                   = {IOImat redo1 IOImatCopyChoice SelectPreviousMask};    % The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs
mask1.prog                  = @ioi_brainmask_run;  % A function handle that will be called with the harvested job to run the computation
mask1.vout                  = @ioi_cfg_vout_brainmask; % A function handle that will be called with the harvested job to determine virtual outputs
mask1.help                  = {'Manual segmentation of the brain to provide a mask for fcOIS analysis'};

return

% make IOI.mat available as a dependency
function vout               = ioi_cfg_vout_brainmask(job)
vout                        = cfg_dep; % The dependency object
vout.sname                  = 'IOI.mat'; % Displayed dependency name
vout.src_output             = substruct('.','IOImat'); %{1}); %,'IOImat');
vout.tgt_spec               = cfg_findspec({{'filter','mat','strtype','e'}});

% EOF
