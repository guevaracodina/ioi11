function [window_start_delay] = ioi_dfg_window_start_delay(wsd)
window_start_delay      = cfg_entry;
window_start_delay.tag  = 'window_start_delay';
window_start_delay.name = 'Window start delay';
window_start_delay.strtype  = 'r';
window_start_delay.num = [1 1];
window_start_delay.val  = {wsd};
window_start_delay.help = {'Start time of the array to average, after time zero, in seconds.'};