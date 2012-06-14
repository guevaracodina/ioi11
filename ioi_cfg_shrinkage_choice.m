function shrinkage_choice = ioi_cfg_shrinkage_choice

shrink_x      = cfg_entry;
shrink_x.tag  = 'shrink_x';
shrink_x.name = 'Shrink factor for x dimension';
shrink_x.strtype  = 'i';
shrink_x.num = [1 1];
shrink_x.val  = {2};
shrink_x.help = {'Data reduction factor in x.'};

shrink_y      = cfg_entry;
shrink_y.tag  = 'shrink_y';
shrink_y.name = 'Shrink factor for y dimension';
shrink_y.strtype  = 'i';
shrink_y.num = [1 1];
shrink_y.val = {2};
shrink_y.help = {'Data reduction factor in y.'};

force_shrink_recompute      = cfg_menu;
force_shrink_recompute.tag  = 'force_shrink_recompute';
force_shrink_recompute.name = 'Force shrink recompute';
force_shrink_recompute.labels = {'Yes','No'};
force_shrink_recompute.values = {1,0};
force_shrink_recompute.val  = {0};
force_shrink_recompute.help = {'This option is used when images of one shrunk size are present'
    'But a different size is desired. This forces generating new shrunk images. The old ones'
    'will be kept too, but will not be available from the new IOI.mat. Therefore, one should'
    'usually select the option to place the new IOI.mat in a new folder, so that the old IOI.mat'
    'can still be used to access the first set of shrunk images.'}';


configuration_shrink         = cfg_branch;
configuration_shrink.tag     = 'configuration_shrink';
configuration_shrink.name    = 'Configuration shrinkage';
configuration_shrink.val     = {shrink_x shrink_y force_shrink_recompute};
configuration_shrink.help    = {'Select values.'};

no_shrinkage         = cfg_branch;
no_shrinkage.tag     = 'no_shrinkage';
no_shrinkage.name    = 'No shrinkage';
no_shrinkage.val     = {};
no_shrinkage.help    = {};

shrinkage_choice        = cfg_choice;
shrinkage_choice.name   = 'Choose shrinkage method';
shrinkage_choice.tag    = 'shrinkage_choice';
shrinkage_choice.values = {no_shrinkage,configuration_shrink};
shrinkage_choice.val    = {configuration_shrink};
shrinkage_choice.help   = {'Choose whether to shrink the data. Images will then be stored. And could be reused later.'}';
