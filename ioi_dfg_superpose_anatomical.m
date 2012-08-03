function superpose_anatomical = ioi_dfg_superpose_anatomical
superpose_anatomical      = cfg_menu;
superpose_anatomical.tag  = 'superpose_anatomical';
superpose_anatomical.name = 'Superpose on anatomical image';
superpose_anatomical.labels = {'Yes','No'};
superpose_anatomical.values = {1,0};
superpose_anatomical.val  = {1};
superpose_anatomical.help = {'Superpose on anatomical image'}';
