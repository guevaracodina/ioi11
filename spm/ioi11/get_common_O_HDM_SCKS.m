function [O HPF LPF] = get_common_O_HDM_SCKS(job)
%select a subset of sessions
if isfield(job.session_choice,'select_sessions')
    O.all_sessions = 0;
    O.selected_sessions = job.session_choice.select_sessions.selected_sessions;
else
    O.all_sessions = 1;
end
%select a subset of ROIs
if isfield(job.ROI_choice,'select_ROIs')
    O.all_ROIs = 0;
    O.selected_ROIs = job.ROI_choice.select_ROIs.selected_ROIs;
else
    O.all_ROIs = 1;
end
%Hemodynamic model choice
O.PhysioModel_Choice = job.PhysioModel_Choice;
%modalities to include
O.includeHbR = job.includeHbR;
O.includeHbO = job.includeHbO;
O.includeHbT = job.includeHbT;
O.includeFlow = job.includeFlow;
%Baseline choice
if isfield(job.baseline_choice,'baseline_percentile_choice')
    O.baseline_choice = 1;
    O.baseline_correction{1} = job.baseline_choice.baseline_percentile_choice.baseline_percentile_HbR/100;
    O.baseline_correction{2} = job.baseline_choice.baseline_percentile_choice.baseline_percentile_HbT/100;
    O.baseline_correction{3} = job.baseline_choice.baseline_percentile_choice.baseline_percentile_flow/100;
else
    if isfield(job.baseline_choice,'baseline_offset_choice')
        O.baseline_choice = 2;
        O.baseline_correction{1} = job.baseline_choice.baseline_offset_choice.baseline_offset_HbR;
        O.baseline_correction{2} = job.baseline_choice.baseline_offset_choice.baseline_offset_HbT;
        O.baseline_correction{3} = job.baseline_choice.baseline_offset_choice.baseline_offset_flow;
    else
        O.baseline_choice = 0;
        O.baseline_correction = [];
    end
end
%Temporal filters:
%HPF
if isfield(job.hpf_butter,'hpf_butter_On')
    HPF.hpf_butter_On = 1;
    HPF.hpf_butter_freq = job.hpf_butter.hpf_butter_On.hpf_butter_freq;
    HPF.hpf_butter_order = job.hpf_butter.hpf_butter_On.hpf_butter_order;
else
    HPF.hpf_butter_On = 0;
end
%LPF
if isfield(job.lpf_choice,'lpf_gauss_On')
    LPF.lpf_gauss_On = 1;
    LPF.fwhm1 = job.lpf_choice.lpf_gauss_On.fwhm1;
else
    LPF.lpf_gauss_On = 0;
end