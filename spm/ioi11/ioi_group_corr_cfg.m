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

% CONTROL / TREATMENT SESSIONS ASKED INSIDE THE SUBJECT LOOP FOR MANUAL MODE

% Control sessions
control_sessions            = cfg_entry;
control_sessions.name       = 'Control sessions';	% The displayed name
control_sessions.tag        = 'control_sessions';	% file names
control_sessions.strtype    = 'e';                  % Real numbers
control_sessions.val        = { {[1 2 3]; [1 2]; [1 2]; [1 2]; [1 2]; [1 2]} }; % Default value
control_sessions.help       = {'Choose the number of control sessions (Before 4-AP injection).' 
    'Enter the control session numbers as a cell of vectors.'  
    'For more info see [IOI.dir.dir_subj_raw ''resumeExp.txt'']'};

% Treatment sessions
treatment_sessions          = cfg_entry;
treatment_sessions.name     = 'Treatment sessions'; % The displayed name
treatment_sessions.tag      = 'treatment_sessions'; % file names
treatment_sessions.strtype  = 'e';                  % Real numbers
treatment_sessions.val      = { {[4 5]; [3 4]; [3]; [3 4]; [3 4]; [3 4]} }; % Default value
treatment_sessions.help     = {'Choose the treatment sessions (After 4-AP injection).'
    'Enter the 4-AP session numbers as a cell of vectors.' 
    'For more info see [IOI.dir.dir_subj_raw ''resumeExp.txt'']'};

AutoSessions                = cfg_branch;
AutoSessions.tag            = 'AutoSessions';
AutoSessions.name           = 'Automatic session selection'; 
AutoSessions.val            = {control_sessions treatment_sessions};
AutoSessions.help           = {'Automatic selection of control and 4-AP sessions.'}';
        
ManualSessions              = cfg_branch;
ManualSessions.tag          = 'ManualSessions';
ManualSessions.name         = 'Manual session selection: GUI'; 
ManualSessions.val          = {};
ManualSessions.help         = {'Manual selection of control and 4-AP sessions.'}';

AutoSessionChoice           = cfg_choice;
AutoSessionChoice.name      = 'Choose session selection mode';
AutoSessionChoice.tag       = 'AutoSessionChoice';
AutoSessionChoice.values    = {ManualSessions AutoSessions}; 
AutoSessionChoice.val       = {AutoSessions}; 
AutoSessionChoice.help      = {'Choose whether to select sessions manually or automatically'}'; 

% Paired seeds
paired_seeds                = cfg_entry;
paired_seeds.name           = 'Paired seeds';       % The displayed name
paired_seeds.tag            = 'paired_seeds';       % file names
paired_seeds.strtype        = 'r';                  % Real numbers
paired_seeds.num            = [Inf 2];              % Number of inputs required
paired_seeds.val            = {[(1:2:12)', (2:2:12)']}; % Default value
paired_seeds.help           = {'Choose the pairs of seeds to compare. Usually:' 
                                '  1,2: Frontal'
                                '  3,4: Motor'
                                '  5,6: Cingulate'
                                '  7,8: Somatosensory'
                                ' 9,10: Retrosplenial'
                                '11,12: Visual'};
                            
% Paired t-test
ttest1                      = cfg_menu;
ttest1.tag                  = 'ttest1';
ttest1.name                 = 't-test';
ttest1.labels               = {'No','Yes'};
ttest1.values               = {0,1};
ttest1.val                  = {1};
ttest1.help                 = {'Perform a paired-sample t-test'}';

% Wilcoxon rank sum test
wilcoxon1                   = cfg_menu;
wilcoxon1.tag               = 'wilcoxon1';
wilcoxon1.name              = 'Wilcoxon test';
wilcoxon1.labels            = {'No','Yes'};
wilcoxon1.values            = {0,1};
wilcoxon1.val               = {1};
wilcoxon1.help              = {'Perform a Wilcoxon rank sum test'}';

% alpha significance level
alpha                       = cfg_entry;
alpha.name                  = 'alpha';              % The displayed name
alpha.tag                   = 'alpha';              % file names
alpha.strtype               = 'r';                  % Real numbers
alpha.num                   = [1 1];                % Number of inputs required
alpha.val                   = {0.05};               % Default value
alpha.help                  = {'Performs the test at the significance level (100*alpha)%.' 
    'alpha must be a scalar'};

% Correlation on 1st derivative
derivative              = cfg_menu;
derivative.tag          = 'derivative';
derivative.name         = '1st derivative';
derivative.labels       = {'No', 'Yes'};
derivative.values       = {0, 1};
derivative.val          = {1};                      % Default value
derivative.help         = {'Choose whether to perform correlation analysis on 1st derivative of seeds/pixels time-course'}';

% Show standard error bar
stderror                   = cfg_menu;
stderror.tag               = 'stderror';
stderror.name              = 'Std. error bars';
stderror.labels            = {'No','Yes'};
stderror.values            = {0,1};
stderror.val               = {1};
stderror.help              = {'Show standard error bars: sigma/sqrt(N)'}';

% Select directory to save global results
parent_results_dir          = cfg_files;
parent_results_dir.tag      = 'parent_results_dir';
parent_results_dir.name     = 'Top directory to save group results';
parent_results_dir.filter   = 'dir';
parent_results_dir.num      = [1 1];
parent_results_dir.help     = {'Select the directory where consolidated results will be saved.'}';

% Generate / save figures
[generate_figures ...
    save_figures]           = ioi_dfg_generate_figures;

% Executable Branch
group_corr1                 = cfg_exbranch; % This is the branch that has information about how to run this module
group_corr1.name            = 'Bilateral correlation group comparison'; % The display name
group_corr1.tag             = 'group_corr1'; %Very important: tag is used when calling for execution
group_corr1.val             = {IOImat redo1 IOImatCopyChoice IC AutoSessionChoice paired_seeds ttest1 wilcoxon1 alpha derivative stderror parent_results_dir generate_figures save_figures};    % The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs
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
