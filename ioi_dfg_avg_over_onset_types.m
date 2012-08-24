function avg_over_onset_types = ioi_dfg_avg_over_onset_types
avg_over_onset_types      = cfg_menu;
avg_over_onset_types.tag  = 'avg_over_onset_types';
avg_over_onset_types.name = 'Average over selected onset types';
avg_over_onset_types.labels = {'Yes','No'};
avg_over_onset_types.values = {1,0};
avg_over_onset_types.val  = {0};
avg_over_onset_types.help = {'Average over selected onset types.'
    'With this option, a global average over stims of the different'
    'types that are selected will be produced.'}';
