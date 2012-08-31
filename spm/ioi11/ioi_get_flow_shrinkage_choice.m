function [flow_shrinkage_choice SH] = ioi_get_flow_shrinkage_choice(job)
if isfield(job.flow_shrinkage_choice,'configuration_flow_shrink')
    SH.shrink_x = job.flow_shrinkage_choice.configuration_flow_shrink.shrink_x;
    SH.shrink_y = job.flow_shrinkage_choice.configuration_flow_shrink.shrink_y;
    try
        SH.force_flow_shrink_recompute = job.flow_shrinkage_choice.configuration_flow_shrink.force_flow_shrink_recompute;
    catch
        SH.force_flow_shrink_recompute = 0;
    end
    flow_shrinkage_choice = 1;
else
    flow_shrinkage_choice = 0;
    SH = [];
end
