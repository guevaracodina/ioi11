function which_onset_type = ioi_cfg_which_onset_type
which_onset_type         = cfg_entry; 
which_onset_type.name    = 'Enter onset type(s) to use';
which_onset_type.tag     = 'which_onset_type';       
which_onset_type.strtype = 'r';
which_onset_type.num     = [1 Inf];     
which_onset_type.val     = {1};
which_onset_type.help    = {'Enter which onset type(s) (relevant if there are'
    'several onset types.'
    'Enter (a list of) ordinal number(s) indicating the desired onset type(s) in the onset type list.'}';
