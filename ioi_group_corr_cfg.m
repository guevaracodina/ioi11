function group_corr1 = ioi_group_corr_cfg
% Graphical interface configuration function for ioi_group_corr_run
%_______________________________________________________________________________
% Copyright (C) 2012 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________

% Select IOI.mat
IOImat                      = ioi_dfg_IOImat(1);
% Force processing
redo1                       = ioi_dfg_redo(0);
% IOI copy/overwrite method
IOImatCopyChoice            = ioi_dfg_IOImatCopyChoice('groupCorr');
% Colors to include (OD,HbO,HbR,HbT,Flow)
IC                          = ioi_dfg_include_colors(0,1,1,1,1);
% Generate / save figures
[generate_figures ...
    save_figures]           = ioi_dfg_generate_figures;

% CONTROL / TREATMENT SESSIONS SHOULD BE ASKED INSIDE THE SUBJECT LOOP

% Control sessions
control_sessions            = cfg_entry;
control_sessions.name       = 'Control sessions';	% The displayed name
control_sessions.tag        = 'control_sessions';	% file names
control_sessions.strtype    = 'r';                  % Real numbers
control_sessions.num        = [1 Inf];              % Number of inputs required
control_sessions.val        = {1};                  % Default value
control_sessions.help       = {'Choose the number of control sessions (Before 4-AP injection). See resumeExp.txt'};

% Treatment sessions
treatment_sessions          = cfg_entry;
treatment_sessions.name     = 'Treatment sessions'; % The displayed name
treatment_sessions.tag      = 'treatment_sessions'; % file names
treatment_sessions.strtype  = 'r';                  % Real numbers
treatment_sessions.num      = [1 Inf];              % Number of inputs required
treatment_sessions.val      = {2};                  % Default value
treatment_sessions.help     = {'Choose the number of treatment sessions (After 4-AP injection). See resumeExp.txt'};

% Paired seeds
paired_seeds                = cfg_entry;
paired_seeds.name           = 'Paired seeds';       % The displayed name
paired_seeds.tag            = 'paired_seeds';       % file names
paired_seeds.strtype        = 'r';                  % Real numbers
paired_seeds.num            = [Inf 2];              % Number of inputs required
paired_seeds.val            = {[(1:2:12)', (2:2:12)']};                  % Default value
paired_seeds.help           = {'Choose the pairs of seeds to compare'};

% alpha significance level
alpha                       = cfg_entry;
alpha.name                  = 'alpha';            % The displayed name
alpha.tag                   = 'alpha';             % file names
alpha.strtype               = 'r';                  % Real numbers
alpha.num                   = [1 1];                % Number of inputs required
alpha.val                   = {0.05};               % Default value
alpha.help                  = {'Performs the test at the significance level (100*alpha)%.' 
    ' alpha must be a scalar'};

% Select directory to save global results
parent_results_dir          = cfg_files;
parent_results_dir.tag      = 'parent_results_dir';
parent_results_dir.name     = 'Top directory to save group results';
parent_results_dir.filter   = 'dir';
parent_results_dir.num      = [1 1];
parent_results_dir.help     = {'Select the directory where consolidated results will be saved.'}';

% Executable Branch
group_corr1                 = cfg_exbranch; % This is the branch that has information about how to run this module
group_corr1.name            = 'Bilateral correlation group comparison'; % The display name
group_corr1.tag             = 'group_corr1'; %Very important: tag is used when calling for execution
group_corr1.val             = {IOImat redo1 IOImatCopyChoice IC  paired_seeds alpha parent_results_dir generate_figures save_figures};    % The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs
group_corr1.prog            = @ioi_group_corr_run; % A function handle that will be called with the harvested job to run the computation
group_corr1.vout            = @ioi_cfg_vout_group_corr; % A function handle that will be called with the harvested job to determine virtual outputs
group_corr1.help            = {'Gets the correlation between each seed and its contralateral homologue. Then performs a paired t-test for each seed set, to have a group comparison.'}';

return

% Make IOI.mat available as a dependency
function vout = ioi_cfg_vout_group_corr(job)
vout                    = cfg_dep;                  % The dependency object
vout.sname              = 'IOI.mat';                % Displayed dependency name
vout.src_output         = substruct('.','IOImat');  %{1}); %,'IOImat');
vout.tgt_spec           = cfg_findspec({{'filter','mat','strtype','e'}});

% EOF
