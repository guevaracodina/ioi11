function [shrinkage_choice SH] = ioi_get_shrinkage_choice(job)
if isfield(job.shrinkage_choice,'configuration_shrink')
    SH.shrink_x = job.shrinkage_choice.configuration_shrink.shrink_x;
    SH.shrink_y = job.shrinkage_choice.configuration_shrink.shrink_y;
    try
        SH.force_shrink_recompute = job.shrinkage_choice.configuration_shrink.force_shrink_recompute;
    catch
        SH.force_shrink_recompute = 0;
    end
    shrinkage_choice = 1;
else
    shrinkage_choice = 0;
end
