function ROC1 = ioi_ROC_cfg
%_______________________________________________________________________
% Copyright (C) 2011 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal

IOImat = ioi_cfg_IOImat(1);
redo1 = ioi_dfg_redo(0);
ROImat = ioi_cfg_ROImat(1);
IOImatCopyChoice = ioi_cfg_IOImatCopyChoice('ROC');
ROI_choice = ioi_cfg_ROI_choice;
session_choice = ioi_cfg_session_choice;
[generate_figures save_figures] = ioi_cfg_generate_figures;

% Executable Branch
ROC1      = cfg_exbranch;       % This is the branch that has information about how to run this module
ROC1.name = 'ROC curves';             % The display name
ROC1.tag  = 'ROC1'; %Very important: tag is used when calling for execution
ROC1.val  = {IOImat redo1 IOImatCopyChoice session_choice ROI_choice ...
     generate_figures save_figures};    % The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs
ROC1.prog = @ioi_ROC_run;  % A function handle that will be called with the harvested job to run the computation
ROC1.vout = @ioi_cfg_vout_ROC; % A function handle that will be called with the harvested job to determine virtual outputs
ROC1.help = {'Receiver operating characteristic curves.'};

return

%make IOI.mat available as a dependency
function vout = ioi_cfg_vout_ROC(job)
vout = cfg_dep;                     % The dependency object
vout.sname      = 'IOI.mat';       % Displayed dependency name
vout.src_output = substruct('.','IOImat'); %{1}); %,'IOImat');
%substruct('()',{1}); % The output subscript reference. This could be any reference into the output variable created during computation
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
