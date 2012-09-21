function extract_roi1 = ioi_extract_roi_time_series_cfg
%_______________________________________________________________________________
% Copyright (C) 2011 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________

% Select IOI.mat
IOImat                  = ioi_dfg_IOImat(1);
% Force processing
redo1                   = ioi_dfg_redo(0);
% IOI copy/overwrite method
IOImatCopyChoice        = ioi_dfg_IOImatCopyChoice('Series');
% Choose ROI selection method (all/selected)
ROI_choice              = ioi_dfg_ROI_choice;
% Choose session selection method (all/selected)
session_choice          = ioi_dfg_session_choice;
% Colors to include
IC                      = ioi_dfg_include_colors(0,1,1,1,1);

% Extract mean signal (from the brain mask pixels)
extractBrainMask        = cfg_menu;
extractBrainMask.tag    = 'extractBrainMask';
extractBrainMask.name   = 'Extract brain mask signal';
extractBrainMask.labels = {'Yes','No'};
extractBrainMask.values = {1,0};
extractBrainMask.val    = {0};   % Default value = 0
extractBrainMask.help   = {'Extract mean signal from the non-masked brain pixels'};
% Extract mean signal (from the brain mask pixels) -- this time using an
% activation map
activMask_choice = ioi_dfg_activation_mask_choice(1);
[generate_figures save_figures] = ioi_dfg_generate_figures;

% Executable Branch
extract_roi1            = cfg_exbranch;       % This is the branch that has information about how to run this module
extract_roi1.name       = 'Extract ROI/seed';             % The display name
extract_roi1.tag        = 'extract_roi1'; %Very important: tag is used when calling for execution
extract_roi1.val        = {IOImat redo1 IOImatCopyChoice session_choice ...
    ROI_choice IC extractBrainMask activMask_choice ...
    generate_figures save_figures};    % The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs
extract_roi1.prog       = @ioi_extract_roi_time_series_run;  % A function handle that will be called with the harvested job to run the computation
extract_roi1.vout       = @ioi_cfg_vout_extract_roi; % A function handle that will be called with the harvested job to determine virtual outputs
extract_roi1.help       = {'Create regions of interest.'};

return

%make IOI.mat available as a dependency
function vout = ioi_cfg_vout_extract_roi(job)
vout                    = cfg_dep; % The dependency object
vout.sname              = 'IOI.mat'; % Displayed dependency name
vout.src_output         = substruct('.','IOImat'); %{1}); %,'IOImat');
%substruct('()',{1}); % The output subscript reference. This could be any reference into the output variable created during computation
vout.tgt_spec           = cfg_findspec({{'filter','mat','strtype','e'}});

% EOF
