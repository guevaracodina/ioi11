function seeds1 = ioi_seeds_cfg
% Graphical interface configuration function for ioi_seeds_run
%_______________________________________________________________________
% Copyright (C) 2010 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%______________________________________________________________________


IOImat = ioi_dfg_IOImat(1);
redo1 = ioi_dfg_redo(0);
IOImatCopyChoice = ioi_dfg_IOImatCopyChoice('Seeds');

RemovePreviousSeed.tag      = 'RemovePreviousSeed';
RemovePreviousSeed.name     = 'Treatment of previous seeds';
RemovePreviousSeed.labels   = {'Keep','Remove'};
RemovePreviousSeed.values   = {0,1};
RemovePreviousSeed.val      = {0};
RemovePreviousSeed.help     = {'If previous seeds had been selected, choose whether to keep'
    'or to remove this information.'}';

select_names        = cfg_menu;
select_names.tag    = 'select_names';
select_names.name   = 'Select names for seeds';
select_names.labels = {'Yes','No'};
select_names.values = {1,0};
select_names.val    = {1};
select_names.help   = {'Option for user to manually enter names of seeds.'
    'If No is selected, seed names will be a number (enumeration).'}';

ArraySeed           = cfg_entry;
ArraySeed.name      = 'Size of array of seeds';
ArraySeed.tag       = 'ArraySeed';       
ArraySeed.strtype   = 'r';
ArraySeed.val{1}    = [2 3];
ArraySeed.num       = [1 2];     
ArraySeed.help      = {'Enter array (number of rows by number of columns) '
    'for automatically generated array of seeds.'
    'Make sure the product of rows by columns (total number of seeds)'
    'is less than 100.'
    'Names of seeds will include array positions, increasing from down to up'
    'vertically and from left to right horizontally.'}'; 

AutoSeed        = cfg_branch;
AutoSeed.tag    = 'AutoSeed';
AutoSeed.name   = 'Automatic seed selection'; 
AutoSeed.val    = {ArraySeed};
AutoSeed.help   = {'Automatic seed selection.'}';
        
ManualSeed      = cfg_branch;
ManualSeed.tag  = 'ManualSeed';
ManualSeed.name = 'Manual seed selection: graphical tool'; 
ManualSeed.val  = {};
ManualSeed.help = {'Manual seed selection: graphical tool.'}';

ManualEnterSeed         = cfg_branch;
ManualEnterSeed.tag     = 'ManualEnterSeed';
ManualEnterSeed.name    = 'Manual seed selection: coordinate entry'; 
ManualEnterSeed.val     = {};
ManualEnterSeed.help    = {'Manual seed selection: coordinate entry.'}';

AutoSeedChoice          = cfg_choice;
AutoSeedChoice.name     = 'Choose seed generation method';
AutoSeedChoice.tag      = 'AutoSeedChoice';
AutoSeedChoice.values   = {ManualSeed ManualEnterSeed AutoSeed}; 
AutoSeedChoice.val      = {ManualSeed}; 
AutoSeedChoice.help     = {'Choose whether to generate seed manually or'
        'automatically'}'; 
        
% Executable Branch
seeds1      = cfg_exbranch;       % This is the branch that has information about how to run this module
seeds1.name = 'Create seeds';             % The display name
seeds1.tag  = 'seeds1'; % Very important: tag is used when calling for execution
seeds1.val  = {IOImat redo1 RemovePreviousSeed IOImatCopyChoice select_names AutoSeedChoice};    % The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs
seeds1.prog = @ioi_seeds_run;  % A function handle that will be called with the harvested job to run the computation
seeds1.vout = @ioi_cfg_vout_seeds; % A function handle that will be called with the harvested job to determine virtual outputs
seeds1.help = {'Create circular ROI''s called seeds'};

return

% make IOI.mat available as a dependency
function vout = ioi_cfg_vout_seeds(~)
vout = cfg_dep;                     % The dependency object
vout.sname      = 'IOI.mat';       % Displayed dependency name
vout.src_output = substruct('.','IOImat'); %{1}); %,'IOImat');
%substruct('()',{1}); % The output subscript reference. This could be any reference into the output variable created during computation
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
