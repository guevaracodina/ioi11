function network1 = ioi_network_analyses_cfg
% Graphical interface configuration function for ioi_network_analyses_run
%_______________________________________________________________________________
% Copyright (C) 2012 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________

% Select IOI.mat
IOImat                  = ioi_dfg_IOImat(2);

% String identifying Control (NaCl) group [NC]
controlString           = cfg_entry;
controlString.name      = 'Control group ID';
controlString.tag       = 'controlString';       
controlString.strtype   = 's';
controlString.val       = {'NC'}; 
controlString.num       = [2 2];     
controlString.help      = {'String to identify Control Group.'}'; 

% String identifying treatment (CaCl2) group [CC]
treatmentString         = cfg_entry;
treatmentString.name    = 'Treatment group ID';
treatmentString.tag     = 'treatmentString';       
treatmentString.strtype = 's';
treatmentString.val     = {'CC'}; 
treatmentString.num     = [2 2];     
treatmentString.help    = {'String to identify Treatment Group.'}'; 

% Force processing
redo1                   = ioi_dfg_redo(0);
% IOI copy/overwrite method
IOImatCopyChoice        = ioi_dfg_IOImatCopyChoice('network');
% Choose ROI selection method (all/selected)
ROI_choice              = ioi_dfg_ROI_choice;
% Choose session selection method (all/selected)
session_choice          = ioi_dfg_session_choice;
% Colors to include (OD,HbO,HbR,HbT,Flow)
IC                      = ioi_dfg_include_colors(0,1,1,1,1,1);

% Select directory to save global results
results_dir             = cfg_files;
results_dir.tag         = 'results_dir';
results_dir.name        = 'Top directory to save group results';
results_dir.filter      = 'dir';
results_dir.ufilter     = '.*';
results_dir.num         = [1 1];
results_dir.help        = {'Select the directory where figures will be saved.'}';

% Measures Type
measures                = cfg_menu;
measures.tag            = 'measures';
measures.name           = 'Measures';
measures.labels         = {'Efficiency/cost', 'PathDistance/Clustering'};
measures.values         = {1, 2};
measures.val            = {1};
measures.help           = {'Choose whether to perform'
                            '   1: Efficiency/cost measures.'
                            '   2: PathDistance/Clustering measures.'}';

                        % Normalization Type
normalType           	= cfg_menu;
normalType.tag          = 'normalType';
normalType.name         = 'Normalization';
normalType.labels       = {'Raw', 'Z-scores', 'Cost'};
normalType.values       = {0, 1, 2};
normalType.val          = {1};                  % Default value
normalType.help         = { 'Choose normalization type: '
                            '0: Uses raw (Fisher-transformed) correlation coefficients;'
                            '1: Normalizes across subjects (transform to z-scores)'
                            '2: Uses Cost'}';
% Threshold
threshold              	= cfg_entry;
threshold.tag           = 'threshold';         	% file names
threshold.name          = 'Threshold';          % The displayed name
threshold.strtype       = 'r';                  % Real numbers
threshold.num           = [0 1];                % Number of inputs required
threshold.val           = {0.5};                % Default value
threshold.help          = {'Threshold value to compute adjacency matrix'};

% Generate / save figures
[generate_figures ...
    save_figures]       = ioi_dfg_generate_figures;

% Executable Branch
network1                = cfg_exbranch; % This is the branch that has information about how to run this module
network1.name           = 'Network analyses'; % The display name
network1.tag            = 'network1'; %Very important: tag is used when calling for execution
network1.val            = {IOImat controlString treatmentString redo1 IOImatCopyChoice ROI_choice session_choice IC results_dir measures normalType threshold generate_figures save_figures};    % The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs
network1.prog           = @ioi_network_analyses_run; % A function handle that will be called with the harvested job to run the computation
network1.vout           = @ioi_cfg_vout_network_analyses; % A function handle that will be called with the harvested job to determine virtual outputs
network1.help           = {'Network analyses.'}';
return

% Make IOI.mat available as a dependency
function vout = ioi_cfg_vout_network_analyses(job)
vout                    = cfg_dep;                  % The dependency object
vout.sname              = 'IOI.mat';                % Displayed dependency name
vout.src_output         = substruct('.','IOImat');  %{1}); %,'IOImat');
vout.tgt_spec           = cfg_findspec({{'filter','mat','strtype','e'}});

% EOF
