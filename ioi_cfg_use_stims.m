function use_stims = ioi_cfg_use_stims
use_stims      = cfg_entry;
use_stims.tag  = 'use_stims';
use_stims.name = 'Enter an array of stim onset numbers to use';
use_stims.strtype  = 'r';
use_stims.num = [0 Inf];
use_stims.val{1} = '';
use_stims.help = {'Use this option to specify, for example,'
    'that the first 10 onsets should be used, by entering the array 1:10.'
    'These onset numbers will be used for each onset type.'
    'Leave array empty in order to use all available onsets, except possibly '
    'some that are excluded by other mechanisms.'}';
