function hdm1 = ioi_HDM_cfg
%_______________________________________________________________________
% Copyright (C) 2011 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal

IOImat         = cfg_files; %Select NIRS.mat for this subject 
IOImat.name    = 'Select IOI.mat'; % The displayed name
IOImat.tag     = 'IOImat';       %file names
IOImat.filter = 'mat';
IOImat.ufilter = '^IOI.mat$';    
IOImat.num     = [1 Inf];     % Number of inputs required 
IOImat.help    = {'Select IOImat dependency if available. '
    'Otherwise, for each subject, select IOI.mat.'}'; % help text displayed

ROImat         = cfg_files; 
ROImat.name    = 'Select ROI.mat'; 
ROImat.tag     = 'ROImat';      
ROImat.filter = 'mat';
ROImat.ufilter = '^ROI.mat$';
ROImat.val     = {''};
ROImat.num     = [1 Inf];    
ROImat.help    = {'Optional: Select ROImat. This allows working on ROI data even'
    'if the paths are not correct in IOI.mat. If not specified, the ROI.mat '
    'specified in IOI.mat will be used.'
    'If several subjects are run, and ROImat is explicitly specified, then'
    'it should be specified for all subjects, in the same order as the list of IOI.mat'}'; 

redo1      = cfg_menu;
redo1.tag  = 'force_redo';
redo1.name = 'Force processing';
redo1.labels = {'False','True'};
redo1.values = {0,1};
redo1.val  = {0};
redo1.help = {'Force redoing this processing even when it has been done already.'
    'Use option below for treatment of previous ROIs.'}';

only_display      = cfg_menu;
only_display.tag  = 'only_display';
only_display.name = 'Only display';
only_display.labels = {'Yes','No'};
only_display.values = {1,0};
only_display.val  = {0};
only_display.help = {'Only display:'
    'Use this option if this HDM has already been generated.'
    'With this option, it will not be regenerated, but display options can be changed.'
    'It will overwrite figures, so if you want to keep the old figures, you should rename'
    'the old figure directory, but not the location where the HDM.mat structure is located.'}';

IOImatOverwrite         = cfg_branch;
IOImatOverwrite.tag     = 'IOImatOverwrite';
IOImatOverwrite.name    = 'Overwrite IOI.mat structure'; 
IOImatOverwrite.help    = {'Will not copy IOI structure.'
            'This will write over the previous NIRS.mat'}';

NewIOIdir         = cfg_entry;
NewIOIdir.name    = 'New directory for IOI.mat';
NewIOIdir.tag     = 'NewIOIdir';       
NewIOIdir.strtype = 's';
NewIOIdir.val{1}  = 'HDM';
NewIOIdir.num     = [1 Inf];     
NewIOIdir.help    = {'Directory for IOI.mat.'}'; 

IOImatCopy         = cfg_branch;
IOImatCopy.tag     = 'IOImatCopy';
IOImatCopy.name    = 'Create new directory and copy IOI structure there'; 
IOImatCopy.val     = {NewIOIdir};
IOImatCopy.help    = {'Create new directory and copy IOI structure there.'}';
        
%Common to most modules: for creating a new directory and copying IOI.mat
IOImatCopyChoice           = cfg_choice;
IOImatCopyChoice.name      = 'Choose IOI copy method';
IOImatCopyChoice.tag       = 'IOImatCopyChoice';
IOImatCopyChoice.values    = {IOImatOverwrite IOImatCopy}; 
IOImatCopyChoice.val       = {IOImatOverwrite}; 
IOImatCopyChoice.help      = {'Choose whether to overwrite the IOI.mat structure'
            'or to create a new directory'
            'and copy the IOI.mat structure there'}'; 

includeHbR      = cfg_menu;
includeHbR.tag  = 'includeHbR';
includeHbR.name = 'Include HbR';
includeHbR.labels = {'Yes','No'};
includeHbR.values = {1,0};
includeHbR.val  = {1};
includeHbR.help = {'Include modality HbR.'}';

includeHbT      = cfg_menu;
includeHbT.tag  = 'includeHbT';
includeHbT.name = 'Include HbT';
includeHbT.labels = {'Yes','No'};
includeHbT.values = {1,0};
includeHbT.val  = {1};
includeHbT.help = {'Include modality HbT.'}';

includeFlow      = cfg_menu;
includeFlow.tag  = 'includeFlow';
includeFlow.name = 'Include Flow';
includeFlow.labels = {'Yes','No'};
includeFlow.values = {1,0};
includeFlow.val  = {1};
includeFlow.help = {'Include modality Flow.'}';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hemodynamic Model choice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% BF         = cfg_branch;
% BF.tag     = 'BF';
% BF.name    = 'Buxton-Friston';
% BF.val     = {};
% BF.help    = {'Buxton-Friston'};
% 
% ZM         = cfg_branch;
% ZM.tag     = 'ZM';
% ZM.name    = 'Zheng-Mayhew';
% ZM.val     = {};
% ZM.help    = {'Zheng-Mayhew'};
% 
% BH         = cfg_branch;
% BH.tag     = 'BH';
% BH.name    = 'Boas-Huppert';
% BH.val     = {};
% BH.help    = {'Boas-Huppert'};
% 
% Model_Choice        = cfg_choice;
% Model_Choice.name   = 'Choose hemodynamic model';
% Model_Choice.tag    = 'Model_Choice';
% Model_Choice.values = {BF,ZM,BH};
% Model_Choice.val    = {BF};
% Model_Choice.help   = {'Choose hemodynamic model'}';

PhysioModel_Choice        = cfg_menu;
PhysioModel_Choice.name   = 'Choose hemodynamic model';
PhysioModel_Choice.tag    = 'PhysioModel_Choice';
PhysioModel_Choice.labels = {'Buxton-Friston','Zheng-Mayhew','Boas-Huppert'};
PhysioModel_Choice.values = {0,1,2};
PhysioModel_Choice.val    = {0};
PhysioModel_Choice.help   = {'Choose hemodynamic model'}';

all_ROIs         = cfg_branch;
all_ROIs.tag     = 'all_ROIs';
all_ROIs.name    = 'All ROIs';
all_ROIs.val     = {};
all_ROIs.help    = {'All ROIs will be processed'};

selected_ROIs      = cfg_entry;
selected_ROIs.tag  = 'selected_ROIs';
selected_ROIs.name = 'Enter list of ROIs';
selected_ROIs.strtype  = 'r';
selected_ROIs.num = [1 Inf];
selected_ROIs.val{1} = 1;
selected_ROIs.help = {'Enter list of ROIs to process.'};

select_ROIs         = cfg_branch;
select_ROIs.tag     = 'select_ROIs';
select_ROIs.name    = 'Select ROIs';
select_ROIs.val     = {selected_ROIs};
select_ROIs.help    = {'Choose some ROIs to be processed'};

ROI_choice        = cfg_choice;
ROI_choice.name   = 'Choose ROI selection method';
ROI_choice.tag    = 'ROI_choice';
ROI_choice.values = {all_ROIs,select_ROIs};
ROI_choice.val    = {all_ROIs};
ROI_choice.help   = {'Choose ROI selection method'}';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

all_sessions         = cfg_branch;
all_sessions.tag     = 'all_sessions';
all_sessions.name    = 'All sessions';
all_sessions.val     = {};
all_sessions.help    = {'All good enough sessions will be processed'};

selected_sessions      = cfg_entry;
selected_sessions.tag  = 'selected_sessions';
selected_sessions.name = 'Enter list of sessions';
selected_sessions.strtype  = 'r';
selected_sessions.num = [1 Inf];
selected_sessions.val{1} = 1;
selected_sessions.help = {'Enter list of sessions to process.'};

select_sessions         = cfg_branch;
select_sessions.tag     = 'select_sessions';
select_sessions.name    = 'Select sessions';
select_sessions.val     = {selected_sessions};
select_sessions.help    = {'Choose some sessions to be processed'};

session_choice        = cfg_choice;
session_choice.name   = 'Choose session selection method';
session_choice.tag    = 'session_choice';
session_choice.values = {all_sessions,select_sessions};
session_choice.val    = {all_sessions};
session_choice.help   = {'Choose session selection method'}';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

use_onset_amplitudes      = cfg_menu;
use_onset_amplitudes.tag  = 'use_onset_amplitudes';
use_onset_amplitudes.name = 'Use onset amplitudes to weigh the hemodynamic response';
use_onset_amplitudes.labels = {'Yes','No'};
use_onset_amplitudes.values = {1,0};
use_onset_amplitudes.val  = {0};
use_onset_amplitudes.help = {'Use onset amplitudes as parameters to weigh the hemodynamic response.'}';


hpf_butter_freq         = cfg_entry; 
hpf_butter_freq.name    = 'Cutoff frequency for HPF';
hpf_butter_freq.tag     = 'hpf_butter_freq';       
hpf_butter_freq.strtype = 'r';
hpf_butter_freq.num     = [1 1];     
hpf_butter_freq.val     = {0.01};
hpf_butter_freq.help    = {'Enter cutoff frequency in Hz for Butterworth HPF.'};

hpf_butter_order         = cfg_entry; 
hpf_butter_order.name    = 'Order of Butterworth HPF';
hpf_butter_order.tag     = 'hpf_butter_order';       
hpf_butter_order.strtype = 'r';
hpf_butter_order.num     = [1 1];     
hpf_butter_order.val     = {3};
hpf_butter_order.help    = {'Enter order of Butterworth HPF (preferred value = 3).'};

hpf_butter_On         = cfg_branch;
hpf_butter_On.tag     = 'hpf_butter_On';
hpf_butter_On.name    = 'Butterworth HP filter';
hpf_butter_On.val     = {hpf_butter_freq hpf_butter_order}; 
hpf_butter_On.help    = {'Butterworth high-pass filter.'};

hpf_butter_Off         = cfg_branch;
hpf_butter_Off.tag     = 'hpf_butter_Off';
hpf_butter_Off.name    = 'HP filter off';
hpf_butter_Off.val     = {}; 
hpf_butter_Off.help    = {'High pass filter turned off.'};

hpf_butter      = cfg_choice;
hpf_butter.tag  = 'hpf_butter';
hpf_butter.name = 'Butterworth High Pass Filter';
hpf_butter.values = {hpf_butter_On hpf_butter_Off};
hpf_butter.val = {hpf_butter_On};
hpf_butter.help = {'Choose whether to include a Butterworth High Pass Filter.'
        'Parameters are: order (e.g. 3) and frequency (e.g. 0.01 Hz)'}';
  
% ---------------------------------------------------------------------
% lpf Low-pass filter
% ---------------------------------------------------------------------
fwhm1      = cfg_entry;
fwhm1.tag  = 'fwhm1';
fwhm1.name = 'FWHM in seconds';
fwhm1.val = {0.67};
fwhm1.strtype = 'r';  
fwhm1.num     = [1 1]; 
fwhm1.help    = {'FWHM in seconds.'}; 

lpf_gauss_On         = cfg_branch;
lpf_gauss_On.tag     = 'lpf_gauss_On';
lpf_gauss_On.name    = 'Gaussian LP filter';
lpf_gauss_On.val     = {fwhm1}; 
lpf_gauss_On.help    = {'Gaussian low-pass filter '
    '(applied forward then backward so that it does not create a time shift)'}';

lpf_Off         = cfg_branch;
lpf_Off.tag     = 'lpf_Off';
lpf_Off.name    = 'LP filter off';
lpf_Off.val     = {}; 
lpf_Off.help    = {'Low pass filter turned off.'};

lpf_choice      = cfg_choice;
lpf_choice.tag  = 'lpf_choice';
lpf_choice.name = 'Choose Low Pass Filter';
lpf_choice.values = {lpf_gauss_On lpf_Off};
lpf_choice.val = {lpf_gauss_On};
lpf_choice.help = {'Choose whether to include a Low Pass Filter.'
        'Parameters'}';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
no_baseline_correction         = cfg_branch;
no_baseline_correction.tag     = 'no_baseline_correction';
no_baseline_correction.name    = 'No baseline correction';
no_baseline_correction.val     = {};
no_baseline_correction.help    = {'No correction will be made to baseline'};

baseline_percentile_HbR      = cfg_entry;
baseline_percentile_HbR.tag  = 'baseline_percentile_HbR';
baseline_percentile_HbR.name = 'Choose data percentile to set baseline for HbR';
baseline_percentile_HbR.strtype  = 'r';
baseline_percentile_HbR.num = [Inf Inf];
baseline_percentile_HbR.val{1} = 90;
baseline_percentile_HbR.help = {'Enter percentile of data (after filtering) to set baseline at, for HbR.'
    'Enter either a single number,'
    'to be applied to each selected session and ROI, or '
    'enter an array of offsets, by selected session and by selected ROI.'}';

baseline_percentile_HbT      = cfg_entry;
baseline_percentile_HbT.tag  = 'baseline_percentile_HbT';
baseline_percentile_HbT.name = 'Choose data percentile to set baseline for HbT';
baseline_percentile_HbT.strtype  = 'r';
baseline_percentile_HbT.num = [Inf Inf];
baseline_percentile_HbT.val{1} = 10;
baseline_percentile_HbT.help = {'Enter percentile of data (after filtering) to set baseline at, for HbT.'
    'Enter either a single number,'
    'to be applied to each selected session and ROI, or '
    'enter an array of offsets, by selected session and by selected ROI.'}';

baseline_percentile_flow      = cfg_entry;
baseline_percentile_flow.tag  = 'baseline_percentile_flow';
baseline_percentile_flow.name = 'Choose data percentile to set baseline for flow';
baseline_percentile_flow.strtype  = 'r';
baseline_percentile_flow.num = [Inf Inf];
baseline_percentile_flow.val{1} = 10;
baseline_percentile_flow.help = {'Enter percentile of data (after filtering) to set baseline at, for flow.'
    'Enter either a single number,'
    'to be applied to each selected session and ROI, or '
    'enter an array of offsets, by selected session and by selected ROI.'}';

baseline_percentile_choice         = cfg_branch;
baseline_percentile_choice.tag     = 'baseline_percentile_choice';
baseline_percentile_choice.name    = 'Baseline percentile choice';
baseline_percentile_choice.val     = {baseline_percentile_HbR baseline_percentile_HbT baseline_percentile_flow};
baseline_percentile_choice.help    = {'Set baseline to a chosen percentile'}';

baseline_offset_HbR      = cfg_entry;
baseline_offset_HbR.tag  = 'baseline_offset_HbR';
baseline_offset_HbR.name = 'Baseline offset for HbR';
baseline_offset_HbR.strtype  = 'r';
baseline_offset_HbR.num = [Inf Inf];
baseline_offset_HbR.val{1} = 0;
baseline_offset_HbR.help = {'Enter baseline offset for HbR. Enter either a single number,'
    'to be applied to each selected session and ROI, or  '
    'enter an array of offsets, by selected session and by selected ROI.'}';

baseline_offset_HbT      = cfg_entry;
baseline_offset_HbT.tag  = 'baseline_offset_HbT';
baseline_offset_HbT.name = 'Baseline offset for HbT';
baseline_offset_HbT.strtype  = 'r';
baseline_offset_HbT.num = [Inf Inf];
baseline_offset_HbT.val{1} = 0;
baseline_offset_HbT.help = {'Enter baseline offset for HbT. Enter either a single number,'
    'to be applied to each selected session and ROI, or  '
    'enter an array of offsets, by selected session and by selected ROI.'}';

baseline_offset_flow      = cfg_entry;
baseline_offset_flow.tag  = 'baseline_offset_flow';
baseline_offset_flow.name = 'Baseline offset for Flow';
baseline_offset_flow.strtype  = 'r';
baseline_offset_flow.num = [Inf Inf];
baseline_offset_flow.val{1} = 0;
baseline_offset_flow.help = {'Enter baseline offset for flow. Enter either a single number,'
    'to be applied to each selected session and ROI, or  '
    'enter an array of offsets, by selected session and by selected ROI.'}';

baseline_offset_choice         = cfg_branch;
baseline_offset_choice.tag     = 'baseline_offset_choice';
baseline_offset_choice.name    = 'Baseline offset choice';
baseline_offset_choice.val     = {baseline_offset_HbR baseline_offset_HbT baseline_offset_flow};
baseline_offset_choice.help    = {'Offset baseline by a chosen amount'};

baseline_choice        = cfg_choice;
baseline_choice.name   = 'Choose baseline selection method';
baseline_choice.tag    = 'baseline_choice';
baseline_choice.values = {no_baseline_correction,baseline_percentile_choice,baseline_offset_choice};
baseline_choice.val    = {baseline_percentile_choice};
baseline_choice.help   = {'Choose baseline selection method'}';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EM parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Niterations         = cfg_entry; 
Niterations.name    = 'Maximum number of EM iterations';
Niterations.tag     = 'Niterations';       
Niterations.strtype = 'r';
Niterations.num     = [1 1];     
Niterations.val     = {128};
Niterations.help    = {'Maximum number of EM iterations. 128 is the basic number.'
    'Increase to 512 or more if necessary.'}';

dFcriterion         = cfg_entry; 
dFcriterion.name    = 'Convergence criterion';
dFcriterion.tag     = 'dFcriterion';       
dFcriterion.strtype = 'r';
dFcriterion.num     = [1 1];     
dFcriterion.val     = {1};
dFcriterion.help    = {'Convergence criterion on changes of free energy F.'
    'Changes in F less than this value are required for convergence.'
    '1e-2 is the SPM8 default value.'}';

LogAscentRate         = cfg_entry; 
LogAscentRate.name    = 'Initial log ascent rate';
LogAscentRate.tag     = 'LogAscentRate';       
LogAscentRate.strtype = 'r';
LogAscentRate.num     = [1 1];     
LogAscentRate.val     = {-2};
LogAscentRate.help    = {'Initial log ascent rate: control initial rate of movement in parameter space'}';

spm_integrator      = cfg_menu;
spm_integrator.tag  = 'spm_integrator';
spm_integrator.name = 'Choose ODE integrator';
spm_integrator.labels = {'spm_int','spm_int_ode','spm_int_J'};
spm_integrator.values = {'spm_int','spm_int_ode','spm_int_J'};
spm_integrator.val  = {'spm_int'};
spm_integrator.help = {'Choose integrator to use for ordinary differential equations.'
    'spm_int is fastest, spm_int_ode most precise, spm_int_J in between'}';

Mstep_iterations         = cfg_entry; 
Mstep_iterations.name    = 'Maximum number of M step iterations';
Mstep_iterations.tag     = 'Mstep_iterations';       
Mstep_iterations.strtype = 'r';
Mstep_iterations.num     = [1 1];     
Mstep_iterations.val     = {8};
Mstep_iterations.help    = {'Maximum number of M-step iterations. 8 is the standard choice.'}';

EM_parameters         = cfg_branch;
EM_parameters.tag     = 'EM_parameters';
EM_parameters.name    = 'Parameters of EM algorithm';
EM_parameters.val     = {Niterations spm_integrator dFcriterion LogAscentRate Mstep_iterations}; 
EM_parameters.help    = {'Parameters of Expectation Maximization algorithm.'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

show_normalized_parameters      = cfg_menu;
show_normalized_parameters.tag  = 'show_normalized_parameters';
show_normalized_parameters.name = 'show normalized parameters';
show_normalized_parameters.labels = {'Yes','No'};
show_normalized_parameters.values = {1,0};
show_normalized_parameters.val  = {0};
show_normalized_parameters.help = {'Show normalized values parameters.'}';

generate_figures      = cfg_menu;
generate_figures.tag  = 'generate_figures';
generate_figures.name = 'Show figures';
generate_figures.labels = {'Yes','No'};
generate_figures.values = {1,0};
generate_figures.val  = {0};
generate_figures.help = {'Show figures. When selecting this option, the figures will stay opened after the code has completed.'}';

save_figures      = cfg_menu;
save_figures.tag  = 'save_figures';
save_figures.name = 'Save figures';
save_figures.labels = {'Yes','No'};
save_figures.values = {1,0};
save_figures.val  = {0};
save_figures.help = {'Save figures.'}';

show_mse      = cfg_menu;
show_mse.tag  = 'show_mse';
show_mse.name = 'Show MSE';
show_mse.labels = {'Yes','No'};
show_mse.values = {1,0};
show_mse.val  = {1};
show_mse.help = {'Show mean square error on figures.'}';

plot_algebraic_CMRO2      = cfg_menu;
plot_algebraic_CMRO2.tag  = 'plot_algebraic_CMRO2';
plot_algebraic_CMRO2.name = 'Plot algebraic CMRO2';
plot_algebraic_CMRO2.labels = {'Yes','No'};
plot_algebraic_CMRO2.values = {1,0};
plot_algebraic_CMRO2.val  = {1};
plot_algebraic_CMRO2.help = {'Plot algebraic CMRO2.'}';

simuA         = cfg_entry;
simuA.tag     = 'simuA';
simuA.name    = 'Signal Amplitude';
simuA.help    = {'Enter signal amplitude, as a percentage of the BOLD signal (e.g. enter 1 for a 1% amplitude)'};
simuA.strtype = 'e';
simuA.num     = [1 1];
simuA.val     = {1};

simuS         = cfg_entry;
simuS.tag     = 'simuS';
simuS.name    = 'Stimuli to simulate';
simuS.help    = {'Enter array of stimuli types to simulated.'
    'Enter 0 to include all stimuli types.'}';
simuS.strtype = 'e';
simuS.num     = [1 Inf];
simuS.val     = {0};

simuP         = cfg_entry;
simuP.tag     = 'simuP';
simuP.name    = 'Parameters to randomize';
simuP.help    = {'Enter array of parameters to be sampled.'
    'Enter 0 to randomize all parameters.'}';
simuP.strtype = 'e';
simuP.num     = [1 Inf];
simuP.val     = {1};

simuR         = cfg_entry;
simuR.tag     = 'simuR';
simuR.name    = 'Parameter range';
simuR.help    = {'For each parameter specified, enter range to be sampled as a percentage of the prior value.'
    'Parameters will be sampled randomly uniformly between 0.75 and 1.25 times the prior value.'
    'If only one number is specified, it will be applied to all parameters to be sampled.'}';
simuR.strtype = 'e';
simuR.num     = [1 Inf];
simuR.val     = {25};

simuPrior         = cfg_entry;
simuPrior.tag     = 'simuPrior';
simuPrior.name    = 'Prior values of parameters';
simuPrior.help    = {'Enter array of prior values of the previously specified parameters to be sampled.'
    'If nothing is entered, the default prior values will be used.'}';
simuPrior.strtype = 'e';
simuPrior.num     = [0 Inf];
simuPrior.val     = {''};

simuIt         = cfg_entry;
simuIt.tag     = 'simuIt';
simuIt.name    = 'Number of simulations';
simuIt.help    = {'Enter number of simulations'};
simuIt.strtype = 'e';
simuIt.num     = [1 1];
simuIt.val     = {1};

simuNoise           = cfg_menu;
simuNoise.name      = 'Include baseline noise';
simuNoise.tag       = 'simuNoise';
simuNoise.labels    = {'Yes' 'No'};
simuNoise.values    = {1,0};
simuNoise.val       = {1};
simuNoise.help      = {'Include noise background scans; if No, the simulated data will be noiseless, i.e. on 0 background.'}';

simuUpsample         = cfg_entry;
simuUpsample.tag     = 'simuUpsample';
simuUpsample.name    = 'Data upsampling factor';
simuUpsample.help    = {'Enter an upsampling factor (max = 16)'};
simuUpsample.strtype = 'e';
simuUpsample.num     = [1 1];
simuUpsample.val     = {1};

simuYes         = cfg_branch;
simuYes.tag     = 'simuYes';
simuYes.name    = 'HDM on simulated data';
simuYes.val     = {simuIt simuA simuS simuP simuPrior simuR simuUpsample simuNoise };
simuYes.help    = {'Perform HDM on real data'}';

simuNo         = cfg_branch;
simuNo.tag     = 'simuNo';
simuNo.name    = 'HDM on real data';
simuNo.val     = {};
simuNo.help    = {'No simulations; perform HDM on real data'}';

simuOn         = cfg_choice;
simuOn.tag     = 'simuOn';
simuOn.name    = 'Perform simulations';
simuOn.values  = {simuNo simuYes};
simuOn.val     = {simuNo};
simuOn.help    = {'Perform simulations'}';

% Executable Branch
hdm1      = cfg_exbranch;       % This is the branch that has information about how to run this module
hdm1.name = 'HDM on ROI';             % The display name
hdm1.tag  = 'hdm1'; %Very important: tag is used when calling for execution
hdm1.val  = {IOImat ROImat redo1 only_display IOImatCopyChoice session_choice ROI_choice...
     PhysioModel_Choice use_onset_amplitudes includeHbR includeHbT includeFlow ... 
     hpf_butter lpf_choice baseline_choice EM_parameters show_normalized_parameters ...
     generate_figures save_figures show_mse plot_algebraic_CMRO2 simuOn};    % The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs
hdm1.prog = @ioi_HDM_run;  % A function handle that will be called with the harvested job to run the computation
hdm1.vout = @ioi_cfg_vout_HDM; % A function handle that will be called with the harvested job to determine virtual outputs
hdm1.help = {'Run hemodynamic modelling (HDM) on average ROI time courses.'};

return

%make IOI.mat available as a dependency
function vout = ioi_cfg_vout_HDM(job)
vout = cfg_dep;                     % The dependency object
vout.sname      = 'IOI.mat';       % Displayed dependency name
vout.src_output = substruct('.','IOImat'); %{1}); %,'IOImat');
%substruct('()',{1}); % The output subscript reference. This could be any reference into the output variable created during computation
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
