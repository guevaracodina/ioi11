function network1 = ioi_network_analyses_cfg
% Graphical interface configuration function for ioi_network_analyses_run
%_______________________________________________________________________________
% Copyright (C) 2012 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________

% Select IOI.mat
IOImat                  = ioi_dfg_IOImat(2);
% Force processing
redo1                   = ioi_dfg_redo(0);
% IOI copy/overwrite method
IOImatCopyChoice        = ioi_dfg_IOImatCopyChoice('network');
% Choose ROI selection method (all/selected)
ROI_choice              = ioi_dfg_ROI_choice;
% Choose session selection method (all/selected)
session_choice          = ioi_dfg_session_choice;
% Colors to include (OD,HbO,HbR,HbT,Flow)
IC                      = ioi_dfg_include_colors(0,1,1,1,1);
% Generate / save figures
[generate_figures ...
    save_figures]       = ioi_dfg_generate_figures;

eff_cost                = cfg_branch;
eff_cost.tag            = 'eff_cost';
eff_cost.name           = 'Efficiency/cost';
eff_cost.val            = {};
eff_cost.help           = {'Efficiency/cost measures.'};

path_cluster            = cfg_branch;
path_cluster.tag        = 'path_cluster';
path_cluster.name       = 'PathDistance/Clustering';
path_cluster.val        = {};
path_cluster.help       = {'PathDistance/Clustering measures.'};

measures                = cfg_choice;
measures.tag            = 'measures';
measures.name           = 'Measures';
measures.values         = {eff_cost path_cluster};
measures.val            = {eff_cost};
measures.help           = {'Choose whether to perform'
                            '   1: Efficiency/cost measures.'
                            '   2: PathDistance/Clustering measures.'}';

% p-Values
pValue                  = cfg_entry;
pValue.tag              = 'pValue';                 % file names
pValue.name             = 'p-value';                % The displayed name
pValue.strtype          = 'r';                      % Real numbers
pValue.num              = [1 1];                    % Number of inputs required
pValue.val              = {0.05};                   % Default value
pValue.help             = {'p-value for testing the hypothesis of no correlation against the alternative that there is a nonzero correlation. If p-value is small, say less than 0.05, then the correlation r is significantly different from zero.'};

% Fisher's z transform
fisherZ                 = cfg_menu;
fisherZ.tag             = 'fisherZ';
fisherZ.name            = 'Fisher''s z transform';
fisherZ.labels          = {'No', 'Yes'};
fisherZ.values          = {0, 1};
fisherZ.val             = {1};                      % Default value
fisherZ.help            = {'Choose whether to perform a Fisher''s z transform of correlation coefficients. The correlation coefficient need to be transformed to the normal distribution by Fisher''s z transform before performing the random effect t-tests'}';

% Seeds correlation matrix
seed2seedCorrMat        = cfg_menu;
seed2seedCorrMat.tag    = 'seed2seedCorrMat';
seed2seedCorrMat.name   = 'seed2seedCorrMat';
seed2seedCorrMat.labels = {'No', 'Yes'};
seed2seedCorrMat.values = {0, 1};
seed2seedCorrMat.val    = {1};                      % Default value
seed2seedCorrMat.help   = {'Choose whether to compute a seed-to-seed correlation matrix'}';

% Correlation on 1st derivative
derivative              = cfg_menu;
derivative.tag          = 'derivative';
derivative.name         = '1st derivative';
derivative.labels       = {'No', 'Yes'};
derivative.values       = {0, 1};
derivative.val          = {1};                      % Default value
derivative.help         = {'Choose whether to perform correlation analysis on 1st derivative of seeds/pixels time-course'}';

% Executable Branch
network1                = cfg_exbranch; % This is the branch that has information about how to run this module
network1.name           = 'Network analyses'; % The display name
network1.tag            = 'network1'; %Very important: tag is used when calling for execution
network1.val            = {IOImat redo1 IOImatCopyChoice ROI_choice session_choice IC measures pValue fisherZ seed2seedCorrMat derivative generate_figures save_figures};    % The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs
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
