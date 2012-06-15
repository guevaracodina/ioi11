function [window_before window_after window_offset] = ioi_cfg_window(bef,aft,ofst)
window_before      = cfg_entry;
window_before.tag  = 'window_before';
window_before.name = 'Window before';
window_before.strtype  = 'r';
window_before.num = [1 1];
window_before.val  = {bef};
window_before.help = {'Size of window to keep prior to each stimulation onset, in seconds.'};

window_after      = cfg_entry;
window_after.tag  = 'window_after';
window_after.name = 'Window after';
window_after.strtype  = 'r';
window_after.num = [1 1];
window_after.val  = {aft};
window_after.help = {'Size of window to keep after each stimulation onset, in seconds.'};

window_offset      = cfg_entry;
window_offset.tag  = 'window_offset';
window_offset.name = 'Window offset';
window_offset.strtype  = 'r';
window_offset.num = [1 1];
window_offset.val  = {ofst};
window_offset.help = {'To look back in time, include an offset in seconds. '
    'A positive number corresponds to a shift back in time. '
    'This works in combination with variables window_after and window_before:'
    'Time 0 will be negative window_offset, window_after starts at time 0, and window_before ends at time 0.'};
