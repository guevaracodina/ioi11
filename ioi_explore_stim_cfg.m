function explore_stim1 = ioi_explore_stim_cfg
%_______________________________________________________________________
% Copyright (C) 2011 LIOM Laboratoire d'Imagerie Optique et Mol�culaire
%                    �cole Polytechnique de Montr�al

IOImat = ioi_dfg_IOImat(1);
redo1 = ioi_dfg_redo(0);
IOImatCopyChoice = ioi_dfg_IOImatCopyChoice('Mean');

%%%%%%%%%%%

electro_stims         = cfg_branch;
electro_stims.tag     = 'electro_stims';
electro_stims.name    = 'Stimulations from electrophysiology';
electro_stims.val     = {};
electro_stims.help    = {'Electrophysiology information'
    'Taken from IOI.Sess().U, which is generated by GLM module.'
    'Thus to use this option, GLM module must be run first.'}';

default_stims         = cfg_branch;
default_stims.tag     = 'default_stims';
default_stims.name    = 'Default stimulations';
default_stims.val     = {};
default_stims.help    = {'The names, onsets, durations fields in IOI.sess_res will be used'};

stim_choice        = cfg_choice;
stim_choice.name   = 'Choose stimulation selection method';
stim_choice.tag    = 'stim_choice';
stim_choice.values = {default_stims,electro_stims};
stim_choice.val    = {default_stims};
stim_choice.help   = {'Choose stimulation selection method'}';

session_choice = ioi_dfg_session_choice;
[window_before window_after] = ioi_dfg_window(3,20,0);

normalize_choice      = cfg_menu;
normalize_choice.tag  = 'normalize_choice';
normalize_choice.name = 'Normalization choice';
normalize_choice.labels = {'Median over window before','Time zero'};
normalize_choice.values = {1,2};
normalize_choice.val  = {2};
normalize_choice.help = {'Normalization choice.'}';

[generate_figures save_figures] = ioi_dfg_generate_figures;

add_error_bars      = cfg_menu;
add_error_bars.tag  = 'add_error_bars';
add_error_bars.name = 'Add error bars';
add_error_bars.labels = {'Yes','No'};
add_error_bars.values = {1,0};
add_error_bars.val  = {0};
add_error_bars.help = {'Add error bars.'}';

% Executable Branch
explore_stim1      = cfg_exbranch;       % This is the branch that has information about how to run this module
explore_stim1.name = 'Explore stimulations';             % The display name
explore_stim1.tag  = 'explore_stim1'; %Very important: tag is used when calling for execution
explore_stim1.val  = {IOImat redo1 IOImatCopyChoice session_choice stim_choice...
    window_after window_before normalize_choice ...
    generate_figures save_figures add_error_bars};    % The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs
explore_stim1.prog = @ioi_explore_stim_run;  % A function handle that will be called with the harvested job to run the computation
explore_stim1.vout = @ioi_cfg_vout_explore_stim; % A function handle that will be called with the harvested job to determine virtual outputs
explore_stim1.help = {'Calculate average over stimulations.'};

return

%make IOI.mat available as a dependency
function vout = ioi_cfg_vout_explore_stim(job)
vout = cfg_dep;                     % The dependency object
vout.sname      = 'IOI.mat';       % Displayed dependency name
vout.src_output = substruct('.','IOImat'); %{1}); %,'IOImat');
%substruct('()',{1}); % The output subscript reference. This could be any reference into the output variable created during computation
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
