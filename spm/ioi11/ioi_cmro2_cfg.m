function cmro2          = ioi_cmro2_cfg
% Graphical interface configuration function for ioi_cmro2_run.
% This code is part of a batch job configuration system for MATLAB. See help
% matlabbatch for a general overview.
%_______________________________________________________________________________
% Copyright (C) 2012 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________

% ------------------------------------------------------------------------------
% General options
% ------------------------------------------------------------------------------
% Select IOI.mat
IOImat              = ioi_dfg_IOImat(1);
% Force processing
redo1               = ioi_dfg_redo(0);
% IOI copy/overwrite method
IOImatCopyChoice    = ioi_dfg_IOImatCopyChoice('CMRO2');
% Choose session selection method (all/selected)
session_choice      = ioi_dfg_session_choice;
% Memomry management type
MemoryManagementMenu= ioi_dfg_MemoryManagement;
% Bandpass filtering
bpf                 = ioi_bpf_cfg(1, [0.009 0.5], 4, 'butter');
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
% Vascular weighting constants
% ------------------------------------------------------------------------------
gammaR              = cfg_entry;
gammaR.tag          = 'gammaR';
gammaR.name         = 'gammaR';
gammaR.strtype      = 'r';
gammaR.num          = [1 1];
gammaR.def          = @(val)ioi_get_defaults('cmro2.gammaR', val{:});
gammaR.help         = {'Vascular weighting constant \gammaR.'};

gammaT              = cfg_entry;
gammaT.tag          = 'gammaT';
gammaT.name         = 'gammaT';
gammaT.strtype      = 'r';
gammaT.num          = [1 1];
gammaT.def          = @(val)ioi_get_defaults('cmro2.gammaT', val{:});
gammaT.help         = {'Vascular weighting constant gammaT'};

constants           = cfg_branch;
constants.tag       = 'constants';
constants.name      = 'Constants';
constants.val       = {gammaT gammaR};
constants.help      = {'The more physiologically plausible range for gammaR & gammaT is around 1 (0.75-1.25)'};
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
% Check previously computed files
% ------------------------------------------------------------------------------
keepFiles           = cfg_menu;
keepFiles.tag       = 'keepFiles';
keepFiles.name      = 'Keep files';
keepFiles.labels    = {'False','True'};
keepFiles.values    = {0,1};
keepFiles.val       = {1};
keepFiles.help      = {'When batch is interrupted, keeps previous NIfTI files'};
% ------------------------------------------------------------------------------


% ------------------------------------------------------------------------------
% Executable Branch
% ------------------------------------------------------------------------------
cmro2               = cfg_exbranch;         % This is the branch that has information about how to run this module
cmro2.name          = 'Compute CMRO2';      % The display name
cmro2.tag           = 'cmro2';              % Very important: tag is used when calling for execution
cmro2.val           = {IOImat redo1 ...     % The items that belong to this branch. 
    keepFiles IOImatCopyChoice...     
    session_choice MemoryManagementMenu...  % All items must be filled before this 
    bpf constants};                         % branch can run or produce virtual outputs
cmro2.prog          = @ioi_cmro2_run;       % A function handle that will be called with the harvested job to run the computation
cmro2.vout          = @ioi_cfg_vout_cmro2;  % A function handle that will be called with the harvested job to determine virtual outputs
cmro2.help          = {'CMRO2 computation. gammaR and gammaT can be varied over a broad range (0.1-5), but the more physiologically plausible range is around 1 (0.75-1.25)'};
% ------------------------------------------------------------------------------
return

% Make IOI.mat available as a dependency
function vout       = ioi_cfg_vout_cmro2(job)
vout                = cfg_dep;              % The dependency object
vout.sname          = 'IOI.mat';            % Displayed dependency name
vout.src_output     = substruct('.','IOImat');
vout.tgt_spec       = cfg_findspec({{'filter','mat','strtype','e'}});
