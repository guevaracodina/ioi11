%% script that cleans up SPM generated files during GLM regression of ROIs/whole
% image
topDir = 'F:\Edgar\Data\IOS_Results\';
clear IOI
% Prompts user to select IOI.mat files
[IOImatList, sts] = cfg_getfile(Inf,'mat','Select IOI mat',[], topDir, '^IOI.mat$');

%% Do the cleanup
if sts
    cleanupOK = false(size(IOImatList));
    for iSubjects = 1:numel(IOImatList)
        % Loads IOI structure
        load(IOImatList{iSubjects});
        % Removes SPM generated files
        cleanupOK(iSubjects) = ioi_fc_GLM_on_ROI_cleanup(IOI);
    end
else
    disp('User cancelled input!')
end

% EOF
