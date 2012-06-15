function [all_ROIs selected_ROIs] = ioi_get_ROIs(job)
if isfield(job.ROI_choice,'select_ROIs')
    all_ROIs = 0;
    selected_ROIs = job.ROI_choice.select_ROIs.selected_ROIs;
else
    all_ROIs = 1;
    selected_ROIs = [];
end