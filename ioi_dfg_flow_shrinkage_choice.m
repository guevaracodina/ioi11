function flow_shrinkage_choice = ioi_dfg_flow_shrinkage_choice

shrink_x      = cfg_entry;
shrink_x.tag  = 'shrink_x';
shrink_x.name = 'Shrink factor for x dimension';
shrink_x.strtype  = 'i';
shrink_x.num = [1 1];
shrink_x.val  = {4};
shrink_x.help = {'Data reduction factor in x.'};

shrink_y      = cfg_entry;
shrink_y.tag  = 'shrink_y';
shrink_y.name = 'Shrink factor for y dimension';
shrink_y.strtype  = 'i';
shrink_y.num = [1 1];
shrink_y.val = {4};
shrink_y.help = {'Data reduction factor in y.'};

force_flow_shrink_recompute      = cfg_menu;
force_flow_shrink_recompute.tag  = 'force_flow_shrink_recompute';
force_flow_shrink_recompute.name = 'Force flow shrink recompute';
force_flow_shrink_recompute.labels = {'Yes','No'};
force_flow_shrink_recompute.values = {1,0};
force_flow_shrink_recompute.val  = {0};
force_flow_shrink_recompute.help = {'This option is used when images of one shrunk size are present'
    'But a different size is desired. This forces generating new shrunk images. The old ones'
    'will be kept too, but will not be available from the new IOI.mat. Therefore, one should'
    'usually select the option to place the new IOI.mat in a new folder, so that the old IOI.mat'
    'can still be used to access the first set of shrunk images.'}';

Configuration_flow_shrinkage         = cfg_branch;
Configuration_flow_shrinkage.tag     = 'configuration_flow_shrink';
Configuration_flow_shrinkage.name    = 'Configuration flow shrinkage';
Configuration_flow_shrinkage.val     = {shrink_x shrink_y force_flow_shrink_recompute};
Configuration_flow_shrinkage.help    = {'Select values.'};

no_shrinkage         = cfg_branch;
no_shrinkage.tag     = 'no_shrinkage';
no_shrinkage.name    = 'No shrinkage';
no_shrinkage.val     = {};
no_shrinkage.help    = {};

flow_shrinkage_choice        = cfg_choice;
flow_shrinkage_choice.name   = 'Choose flow shrinkage method';
flow_shrinkage_choice.tag    = 'flow_shrinkage_choice';
flow_shrinkage_choice.values = {no_shrinkage,Configuration_flow_shrinkage};
flow_shrinkage_choice.val    = {Configuration_flow_shrinkage};
flow_shrinkage_choice.help   = {'Choose whether to shrink the data. Images will then be stored. And could be reused later.'}';
