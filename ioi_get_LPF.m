function LPF = ioi_get_LPF(job)
%LPF
if isfield(job.lpf_choice,'lpf_gauss_On')
    LPF.lpf_gauss_On = 1;
    LPF.fwhm1 = job.lpf_choice.lpf_gauss_On.fwhm1;
    LPF.apply_lpf_on_flow_only = job.lpf_choice.lpf_gauss_On.apply_lpf_on_flow_only;
else
    LPF.lpf_gauss_On = 0;
end