function mask1 = ioi_fcmap1_cfg
% Batch job configuration system for ioi_fcmap1_run
%_______________________________________________________________________
% Copyright (C) 2010 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%______________________________________________________________________

IOImat = ioi_dfg_IOImat(1);
redo1 = ioi_dfg_redo(0);

% Executable Branch
mask1      = cfg_exbranch;       % This is the branch that has information about how to run this module
mask1.name = 'Generate Network Mask';             % The display name
mask1.tag  = 'mask1'; %Very important: tag is used when calling for execution
mask1.val  = {IOImat redo1};    % The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs
mask1.prog = @ioi_brainmask_run;  % A function handle that will be called with the harvested job to run the computation
mask1.vout = @ioi_cfg_vout_brainmask; % A function handle that will be called with the harvested job to determine virtual outputs
mask1.help = {'Manual segmentation of the brain to provide a mask for fcOIS analysis'};

return

% make IOI.mat available as a dependency
function vout = ioi_cfg_vout_brainmask(job)
vout = cfg_dep;                     % The dependency object
vout.sname      = 'IOI.mat';       % Displayed dependency name
vout.src_output = substruct('.','IOImat'); %{1}); %,'IOImat');
%substruct('()',{1}); % The output subscript reference. This could be any reference into the output variable created during computation
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
