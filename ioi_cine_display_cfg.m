function cine_display1 = ioi_cine_display_cfg
%_______________________________________________________________________
% Copyright (C) 2011 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal

% Executable Branch
cine_display1      = cfg_exbranch;       % This is the branch that has information about how to run this module
cine_display1.name = 'Cine Display';             % The display name
cine_display1.tag  = 'cine_display1'; %Very important: tag is used when calling for execution
cine_display1.val  = {};    % The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs
cine_display1.prog = @ioi_cine_display_run;  % A function handle that will be called with the harvested job to run the computation
cine_display1.vout = @ioi_cfg_vout_cine_display; % A function handle that will be called with the harvested job to determine virtual outputs
cine_display1.help = {'Display and explore previously generated 2D movies'}';

return

%make IOI.mat available as a dependency
function vout = ioi_cfg_vout_cine_display(job)
vout = cfg_dep;                     % The dependency object
vout.sname      = 'IOI.mat';       % Displayed dependency name
vout.src_output = substruct('.','IOImat'); %{1}); %,'IOImat');
%substruct('()',{1}); % The output subscript reference. This could be any reference into the output variable created during computation
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
