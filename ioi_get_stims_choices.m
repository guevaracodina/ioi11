function [rmi ust remove_stims_SD] = ioi_get_stims_choices(job)
if isfield(job,'remove_stims')
    rmi = job.remove_stims;
else
    rmi = '';
end
if isfield(job,'use_stims')
    ust = job.use_stims;
else
    ust = '';
end
if isfield(job,'remove_stims_SD')
    remove_stims_SD = job.remove_stims_SD;
else
    remove_stims_SD = 1;
end