function SCKS1 = ioi_SCKS_cfg
%_______________________________________________________________________
% Copyright (C) 2011 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
[IOImat ROImat redo1 IOImatCopyChoice includeHbR includeHbO ...
    includeHbT includeFlow PhysioModel_Choice ROI_choice session_choice baseline_choice ...
    lpf_choice hpf_butter] = ioi_common_fields_SCKS_HDM('SCKS');

only_display      = cfg_menu;
only_display.tag  = 'only_display';
only_display.name = 'Only display';
only_display.labels = {'Yes','No'};
only_display.values = {1,0};
only_display.val  = {0};
only_display.help = {'Only display:'
    'Use this option if this SCKS has already been generated.'
    'With this option, it will not be regenerated, but display options can be changed.'
    'It will overwrite figures, so if you want to keep the old figures, you should rename'
    'the old figure directory, but not the location where the SCKS.mat structure is located.'}';

SCKSnoise      = cfg_menu;
SCKSnoise.tag  = 'SCKSnoise';
SCKSnoise.name = 'Promote parameters and noise to states';
SCKSnoise.labels = {'Yes','No'};
SCKSnoise.values = {1,0};
SCKSnoise.val  = {0};
SCKSnoise.help = {'Promote parameters and noise to time-dependent states.'}';
    
State_annealing         = cfg_entry; 
State_annealing.name    = 'Annealing factor on states';
State_annealing.tag     = 'State_annealing';       
State_annealing.strtype = 'r';
State_annealing.num     = [1 1];     
State_annealing.val     = {0.9995};
State_annealing.help    = {'Enter annealing factor on states. This is a number between 0 and 1.'
    'A value closer to 1 corresponds to more annealing. A value of 0.5 corresponds to very little annealing.'}';

Parameter_annealing         = cfg_entry; 
Parameter_annealing.name    = 'Annealing factor on parameters';
Parameter_annealing.tag     = 'Parameter_annealing';       
Parameter_annealing.strtype = 'r';
Parameter_annealing.num     = [1 1];     
Parameter_annealing.val     = {0.9995};
Parameter_annealing.help    = {'Enter annealing factor on parameters.' 
    'This is a number between 0 and 1.'
    'A value closer to 1 corresponds to more annealing. A value of 0.5 corresponds to very little annealing.'}';

SCKSparams         = cfg_branch;
SCKSparams.tag     = 'SCKSparams';
SCKSparams.name    = 'SCKS parameters';
SCKSparams.val     = {SCKSnoise State_annealing Parameter_annealing}; 
SCKSparams.help    = {'User-controlled SCKS parameters.'};

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

% Executable Branch
SCKS1      = cfg_exbranch;       % This is the branch that has information about how to run this module
SCKS1.name = 'SCKS Deconvolution';             % The display name
SCKS1.tag  = 'SCKS1'; %Very important: tag is used when calling for execution
SCKS1.val  = {IOImat ROImat only_display redo1 IOImatCopyChoice session_choice ROI_choice ...
    PhysioModel_Choice includeHbR includeHbO includeHbT includeFlow hpf_butter lpf_choice baseline_choice ...
    SCKSparams generate_figures save_figures};    % The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs
SCKS1.prog = @ioi_SCKS_run;  % A function handle that will be called with the harvested job to run the computation
SCKS1.vout = @ioi_cfg_vout_SCKS; % A function handle that will be called with the harvested job to determine virtual outputs
SCKS1.help = {'Create regions of interest.'};

return

%make IOI.mat available as a dependency
function vout = ioi_cfg_vout_SCKS(job)
vout = cfg_dep;                     % The dependency object
vout.sname      = 'IOI.mat';       % Displayed dependency name
vout.src_output = substruct('.','IOImat'); %{1}); %,'IOImat');
%substruct('()',{1}); % The output subscript reference. This could be any reference into the output variable created during computation
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
