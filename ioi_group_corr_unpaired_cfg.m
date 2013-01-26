function group_corr2        = ioi_group_corr_unpaired_cfg
% Graphical interface configuration function for ioi_group_corr_unpaired_run
%_______________________________________________________________________________
% Copyright (C) 2012 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________

% Select IOI.mat (2 files minimum)
IOImat                      = ioi_dfg_IOImat(2);
% Force processing
redo1                       = ioi_dfg_redo(1);
% IOI copy/overwrite method
IOImatCopyChoice            = ioi_dfg_IOImatCopyChoice('groupCorrUnpaired');
% Colors to include (OD,HbO,HbR,HbT,Flow)
IC                          = ioi_dfg_include_colors(0,1,1,1,1,1);

% String identifying Control (NaCl) group [NC]
controlString               = cfg_entry;
controlString.name          = 'Control group ID';
controlString.tag           = 'controlString';       
controlString.strtype       = 's';
controlString.val{1}        = 'NC'; 
controlString.num           = [2 2];     
controlString.help          = {'String to identify Control Group.'}'; 

% String identifying treatment (CaCl2) group [CC]
treatmentString             = cfg_entry;
treatmentString.name        = 'Treatment group ID';
treatmentString.tag         = 'treatmentString';       
treatmentString.strtype     = 's';
treatmentString.val{1}      = 'CC'; 
treatmentString.num         = [2 2];     
treatmentString.help        = {'String to identify Treatment Group.'}'; 

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
                            
% Multiple comparisons correc
bonferroni              	= cfg_menu;
bonferroni.tag              = 'bonferroni';
bonferroni.name             = 'Bonferroni correction';
bonferroni.labels           = {'No','Yes'};
bonferroni.values           = {0,1};
bonferroni.val              = {1};
bonferroni.help             = {'Perform Bonferroni correction for multiple comparisons.'}';
                            
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
derivative                  = cfg_menu;
derivative.tag              = 'derivative';
derivative.name             = '1st derivative';
derivative.labels           = {'No', 'Yes'};
derivative.values           = {0, 1};
derivative.val              = {1};                      % Default value
derivative.help             = {'Choose whether to perform correlation analysis on 1st derivative of seeds/pixels time-course'}';

% Correlation on raw data time course (before filtering, downsampling and GLM regression)
rawData                     = cfg_menu;
rawData.tag                 = 'rawData';
rawData.name                = 'raw time course';
rawData.labels              = {'No', 'Yes'};
rawData.values              = {0, 1};
rawData.val                 = {1};                      % Default value
rawData.help                = {'Choose whether to perform correlation analysis on seeds raw time course'}';

% Show standard error bar
stderror                    = cfg_menu;
stderror.tag                = 'stderror';
stderror.name               = 'Std. error bars';
stderror.labels             = {'No','Yes'};
stderror.values             = {0,1};
stderror.val                = {1};
stderror.help               = {'Show standard error bars: sigma/sqrt(N)'}';

% Select directory to save global results
parent_results_dir          = cfg_files;
parent_results_dir.tag      = 'parent_results_dir';
parent_results_dir.name     = 'Top directory to save group results';
parent_results_dir.filter   = 'dir';
parent_results_dir.num      = [1 1];
parent_results_dir.help     = {'Select the directory where consolidated results will be saved.'}';

% Figure size
figSize                     = cfg_entry;
figSize.tag                 = 'figSize';
figSize.name                = 'Figure size';
figSize.strtype             = 'r';
figSize.num                 = [1 2];
figSize.val{1}              = [3.25 3.25];
figSize.help                = {'Enter figure size in inches.'};

% Figure resolution
figRes                      = cfg_entry;
figRes.tag                  = 'figRes';
figRes.name                 = 'Figure resolution';
figRes.strtype              = 'r';
figRes.num                  = [1 1];
figRes.val{1}               = 300;
figRes.help                 = {'Enter figure resolution in dpi [300-1200]'};

% ------------------------------------------------------------------------------
% Choose axis limits
% ------------------------------------------------------------------------------
yLimValue                   = cfg_entry;
yLimValue.tag               = 'yLimValue';
yLimValue.name              = 'Y axis limits';
yLimValue.strtype           = 'r';
yLimValue.num               = [1 2];
yLimValue.val{1}            = [-0.9 1.3];
yLimValue.help              = {'Enter limits for Y axis'};

yLimManual                  = cfg_branch;
yLimManual.tag              = 'yLimManual';
yLimManual.name             = 'Manual Ylim';
yLimManual.val              = {yLimValue};
yLimManual.help             = {'Manual limits for Y axis.'};

yLimAuto                    = cfg_branch;
yLimAuto.tag                = 'yLimAuto';
yLimAuto.name               = 'Auto Ylim';
yLimAuto.val                = {};
yLimAuto.help               = {'Auto limits for Y axis.'};

yLimits                     = cfg_choice;
yLimits.tag                 = 'yLimits';
yLimits.name                = 'Y axis limits';
yLimits.values              = {yLimManual yLimAuto};
yLimits.val                 = {yLimAuto};
yLimits.help                = {'Choose whether to set manual limits to Y axis'};
% ------------------------------------------------------------------------------

% Generate / save figures
[generate_figures ...
    save_figures]           = ioi_dfg_generate_figures;

% Executable Branch
group_corr2                 = cfg_exbranch; % This is the branch that has information about how to run this module
group_corr2.name            = 'Bilateral correlation group comparison (unpaired)'; % The display name
group_corr2.tag             = 'group_corr2'; %Very important: tag is used when calling for execution
group_corr2.val             = {IOImat redo1 IOImatCopyChoice IC controlString ...
    treatmentString paired_seeds bonferroni ttest1 wilcoxon1 alpha derivative rawData ...
    stderror parent_results_dir figSize figRes yLimits generate_figures save_figures};    % The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs
group_corr2.prog            = @ioi_group_corr_unpaired_run; % A function handle that will be called with the harvested job to run the computation
group_corr2.vout            = @ioi_cfg_vout_group_corr_unpaired; % A function handle that will be called with the harvested job to determine virtual outputs
group_corr2.help            = {'Gets the correlation between each seed and its contralateral homologue. Then performs a non-paired t-test for each seed set, to have a group comparison.'}';

return

% Make IOI.mat available as a dependency
function vout = ioi_cfg_vout_group_corr_unpaired(job)
vout                        = cfg_dep;                  % The dependency object
vout.sname                  = 'IOI.mat';                % Displayed dependency name
vout.src_output             = substruct('.','IOImat');  %{1}); %,'IOImat');
vout.tgt_spec               = cfg_findspec({{'filter','mat','strtype','e'}});

% EOF
