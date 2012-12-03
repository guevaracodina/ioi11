function flow1 = ioi_flow_cfg
% Computes blood flow from speckle contrast
%_______________________________________________________________________________
% Copyright (C) 2011 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________


% Select IOI.mat
IOImat              = ioi_dfg_IOImat(1);
% Force processing
redo1               = ioi_dfg_redo(0);
% IOI copy/overwrite method
IOImatCopyChoice    = ioi_dfg_IOImatCopyChoice('Flow');
% Choose session selection method (all/selected)
session_choice      = ioi_dfg_session_choice;

% Integration time of the CCD camera
inttime1            = cfg_entry;
inttime1.tag        = 'integ_time';
inttime1.name       = 'Integration Time';
inttime1.strtype    = 'r';
inttime1.num        = [1 1];
inttime1.def        = @(val)ioi_get_defaults('flow1.T', val{:});
inttime1.help       = {'Camera integration time.'};

% Window size for speckle contrast
windowsize1         = cfg_entry;
windowsize1.tag     = 'window_size';
windowsize1.name    = 'Speckle Window Size';
windowsize1.strtype = 'i';
windowsize1.num     = [1 1];
windowsize1.def     = @(val)ioi_get_defaults('flow1.window_size', val{:});
windowsize1.help    = {'Window size for speckle contrast computation.'};

configuration       = cfg_branch;
configuration.tag   = 'configuration';
configuration.name  = 'Configuration options';
configuration.val   = {inttime1 windowsize1};
configuration.help  = {'Select values.'};

RemoveLC            = cfg_menu;
RemoveLC.tag        = 'RemoveLC';
RemoveLC.name       = 'Remove Laser/Contrast NIfTI images';
RemoveLC.labels     = {'Yes','No'};
RemoveLC.values     = {1,0};
RemoveLC.val        = {1};
RemoveLC.help       = {'After flow images are obtained'
    'Laser and contrast nifti images (LC) are usually no longer required.'}';

% Option to normalize decorrelation velocity Vc = lambda/(2*pi*Tc)
flow_lambda_norm    = ioi_flow_norm_cfg(0);

% Executable Branch
flow1               = cfg_exbranch;       % This is the branch that has information about how to run this module
flow1.name          = 'Compute Flow';             % The display name
flow1.tag           = 'flow1'; %Very important: tag is used when calling for execution
flow1.val           = {IOImat redo1 IOImatCopyChoice session_choice configuration ...
     RemoveLC flow_lambda_norm};    % The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs
flow1.prog          = @ioi_flow_run;  % A function handle that will be called with the harvested job to run the computation
flow1.vout          = @ioi_cfg_vout_flow; % A function handle that will be called with the harvested job to determine virtual outputs
flow1.help          = {'Flow computations.'};

return

%make IOI.mat available as a dependency
function vout       = ioi_cfg_vout_flow(job)
vout = cfg_dep;                     % The dependency object
vout.sname          = 'IOI.mat';       % Displayed dependency name
vout.src_output     = substruct('.','IOImat'); %{1}); %,'IOImat');
%substruct('()',{1}); % The output subscript reference. This could be any reference into the output variable created during computation
vout.tgt_spec       = cfg_findspec({{'filter','mat','strtype','e'}});

% EOF
