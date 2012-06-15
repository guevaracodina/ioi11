function hdm1 = ioi_HDM_cfg
[IOImat ROImat redo1 IOImatCopyChoice IC ...
    PhysioModel_Choice ROI_choice session_choice baseline_choice ...
    lpf_choice hpf_butter] = ioi_common_fields_SCKS_HDM('HDM');
%_______________________________________________________________________
% Copyright (C) 2011 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
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

use_onset_amplitudes      = cfg_menu;
use_onset_amplitudes.tag  = 'use_onset_amplitudes';
use_onset_amplitudes.name = 'Use onset amplitudes to weigh the hemodynamic response';
use_onset_amplitudes.labels = {'Yes','No'};
use_onset_amplitudes.values = {1,0};
use_onset_amplitudes.val  = {0};
use_onset_amplitudes.help = {'Use onset amplitudes as parameters to weigh the hemodynamic response.'}';

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

[generate_figures save_figures] = ioi_cfg_generate_figures;

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
     PhysioModel_Choice use_onset_amplitudes IC ... 
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
