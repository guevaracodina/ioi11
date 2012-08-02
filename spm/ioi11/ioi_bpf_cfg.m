function bpf        = ioi_bpf_cfg(choose_on, freq, order, type)
% Configuration unit for band-pass filter in Matlab batch gui.
% SYNTAX
% bpf               = ioi_dfg_bpf(choose_on, freq, order, type)
% INPUTS
% choose_on         Determines if filter is enabled (1) or not (0)
% freq              2-element vector with the cut-off frequencies
% order             integer defining the order of the filter
% type              @butter, @cheby1, @cheby2, @ellip, @yulewalk
% OUTPUT
% bpf               Executable branch passed to graphical interface
%                   configuration function
%_______________________________________________________________________________
% Copyright (C) 2010 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________
bpf_freq            = cfg_entry; 
bpf_freq.tag        = 'bpf_freq';       
bpf_freq.name       = 'Cutoff frequencies for BPF';
bpf_freq.strtype    = 'r';
bpf_freq.num        = [1 2];     
bpf_freq.val        = {freq};
bpf_freq.help       = {'Enter band-pass frequencies in Hz for BPF.'};

bpf_order           = cfg_entry; 
bpf_order.tag       = 'bpf_order';       
bpf_order.name      = 'Order of BPF';
bpf_order.strtype   = 'r';
bpf_order.num       = [1 1];     
bpf_order.val       = {order};
bpf_order.help      = {'Enter order of BPF (preferred value = 8).'};

bpf_type            = cfg_menu;
bpf_type.tag        = 'bpf_type';
bpf_type.name       = 'Type of band-pass filter';
bpf_type.labels     = {'Butterworth', 'Chebyshev I', 'Chebyshev II', 'Elliptic', 'Yulewalk'};
bpf_type.values     = {'butter', 'cheby1', 'cheby2', 'ellip', 'yulewalk'};
bpf_type.val        = {type};
bpf_type.help       = {'Filter and downsample whole image time-series. It creates a new sub-folder for each session'};

bpf_On              = cfg_branch;
bpf_On.tag          = 'bpf_On';
bpf_On.name         = 'Enable BP filter';
bpf_On.val          = {bpf_freq, bpf_order, bpf_type}; 
bpf_On.help         = {'Band-pass filter.'};

bpf_Off             = cfg_branch;
bpf_Off.tag         = 'bpf_Off';
bpf_Off.name        = 'BP filter off';
bpf_Off.val         = {}; 
bpf_Off.help        = {'Band-pass filter turned off.'};

bpf                 = cfg_choice;
bpf.tag             = 'bpf';
bpf.name            = 'Band Pass Filter';
bpf.values          = {bpf_On bpf_Off};
if choose_on
    bpf.val         = {bpf_On};
else
    bpf.val         = {bpf_Off};
end
bpf.help            = {'Choose whether to include a Band-Pass Filter. Parameters are: order (e.g. 4) and frequency (e.g. [0.009 0.08] Hz)'}';
 
