function extract_roi1 = ioi_extract_roi_time_series_cfg
%_______________________________________________________________________
% Copyright (C) 2011 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
IOImat = ioi_cfg_IOImat(1);
redo1 = ioi_dfg_redo(0);
IOImatCopyChoice = ioi_cfg_IOImatCopyChoice('Series');
ROI_choice = ioi_cfg_ROI_choice;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
session_choice = ioi_cfg_session_choice;

% Colors to include
IC = ioi_cfg_include_colors(0,1,1,1,1);

% Executable Branch
extract_roi1      = cfg_exbranch;       % This is the branch that has information about how to run this module
extract_roi1.name = 'Extract ROI';             % The display name
extract_roi1.tag  = 'extract_roi1'; %Very important: tag is used when calling for execution
extract_roi1.val  = {IOImat redo1 IOImatCopyChoice session_choice ROI_choice IC};    % The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs
extract_roi1.prog = @ioi_extract_roi_time_series_run;  % A function handle that will be called with the harvested job to run the computation
extract_roi1.vout = @ioi_cfg_vout_extract_roi; % A function handle that will be called with the harvested job to determine virtual outputs
extract_roi1.help = {'Create regions of interest.'};

return

%make IOI.mat available as a dependency
function vout = ioi_cfg_vout_extract_roi(job)
vout = cfg_dep;                     % The dependency object
vout.sname      = 'IOI.mat';       % Displayed dependency name
vout.src_output = substruct('.','IOImat'); %{1}); %,'IOImat');
%substruct('()',{1}); % The output subscript reference. This could be any reference into the output variable created during computation
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
