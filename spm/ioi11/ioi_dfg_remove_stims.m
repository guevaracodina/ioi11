function remove_stims = ioi_dfg_remove_stims
remove_stims      = cfg_entry;
remove_stims.tag  = 'remove_stims';
remove_stims.name = 'Enter an array of time points (in seconds) from which to exclude onsets';
remove_stims.strtype  = 'r';
remove_stims.num = [0 Inf];
remove_stims.val{1} = '';
remove_stims.help = {'Onsets occuring within 1 second of any specified' 
    'time point will be removed from the averaging.'
    'The list of onsets kept will be written to IOI.mat'}';