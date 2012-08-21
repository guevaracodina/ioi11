function superpose_anatomical = ioi_dfg_superpose_anatomical
threshold      = cfg_entry;
threshold.tag  = 'threshold';
threshold.name = 'Threshold on t-stats';
threshold.strtype  = 'r';
threshold.num = [1 1];
threshold.val  = {2};
threshold.help = {'Threshold on t-statistic images.'};

SuperposeOn         = cfg_branch;
SuperposeOn.tag     = 'SuperposeOn';
SuperposeOn.name    = 'Superpose on anatomical image';
SuperposeOn.val     = {threshold};
SuperposeOn.help    = {'Superpose on anatomical image.'};

SuperposeOff         = cfg_branch;
SuperposeOff.tag     = 'SuperposeOff';
SuperposeOff.name    = 'Do not superpose';
SuperposeOff.val     = {};
SuperposeOff.help    = {'Do not superpose on anatomical image'};

superpose_anatomical      = cfg_choice;
superpose_anatomical.tag  = 'superpose_anatomical';
superpose_anatomical.name = 'Superpose on anatomical image';
superpose_anatomical.values = {SuperposeOn,SuperposeOff};
superpose_anatomical.val  = {SuperposeOn};
superpose_anatomical.help = {'Superpose on anatomical image'}';
