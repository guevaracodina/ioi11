function [E stim_choice] = ioi_get_E_for_electrophysiology(job)
E = [];
stim_choice=0;
if isfield(job.stim_choice,'electro_stims')
    %choose first or second electrode
    if isfield(job.stim_choice.electro_stims,'seizure_choice')
        E.szOn = 1; E.spkOn = 0;
        electro_stims = job.stim_choice.electro_stims.seizure_detection;
        E.tb = electro_stims.seizure_tb;
        E.ta = electro_stims.seizure_ta;
        %minimal distance between peaks in seconds
        E.dP = electro_stims.seizure_dP;
    end
    if isfield(job.stim_choice.electro_stims,'spike_choice')
        E.szOn = 0; E.spkOn = 1;
        electro_stims = job.stim_choice.electro_stims.spike_detection;
        E.tb = electro_stims.spike_tb;
        E.ta = electro_stims.spike_ta;
        %minimal distance between peaks in seconds
        E.dP = electro_stims.spike_dP/1000;
    end
    E.el_choice = electro_stims.electrophysiology_choice;
    stim_choice = 1;
    %sampling frequency
    E.sf = electro_stims.sf;
    %number of standard deviations
    E.nSD = electro_stims.nSD;
    %minimum standard deviation imposed
    E.mbSD = electro_stims.mbSD;
    
    %HPF
    if isfield(electro_stims.hpf_butter,'electro_hpf_butter_On')
        E.electro_hpf_butter_On = 1;
        E.hpf_butter_freq = electro_stims.hpf_butter.electro_hpf_butter_On.electro_hpf_butter_freq;
        E.hpf_butter_order = electro_stims.hpf_butter.electro_hpf_butter_On.electro_hpf_butter_order;
    else
        E.electro_hpf_butter_On = 0;
    end
    %LPF
    if isfield(electro_stims.electro_lpf_butter,'electro_lpf_butter_On')
        E.electro_lpf_butter_On = 1;
        E.lpf_butter_freq = electro_stims.electro_lpf_butter.electro_lpf_butter_On.electro_lpf_butter_freq;
        E.lpf_butter_order = electro_stims.electro_lpf_butter.electro_lpf_butter_On.electro_lpf_butter_order;
    else
        E.electro_lpf_butter_On = 0;
    end
    E.write_pictures = electro_stims.write_pictures;
    E.use_epilepsy_convention = electro_stims.use_epilepsy_convention;
    E.electrophysiology_onset_name = electro_stims.electrophysiology_onset_name;  
end