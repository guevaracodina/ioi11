function cleanupOK = ioi_fc_GLM_on_ROI_cleanup(IOI)
% Cleanup SPM generated files during the GLM regression. Keeps only the .nii
% files of the regressed ROI/whole image time course.
% SYNTAX
% cleanupOK = ioi_fc_GLM_on_ROI_cleanup(IOI)
% INPUT
% IOI       IOI structure
% OUTPUT
% cleanupOK True if cleanup was succesful.
%_______________________________________________________________________________
% Copyright (C) 2012 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________

try
    tic
    cleanupOK = false;
    % Cell with names of files to be deleted
    files2Delete = {    
                    'beta_0001.hdr'
                    'beta_0001.img'
                    'mask.hdr'
                    'mask.img'
                    'ResMS.hdr'
                    'ResMS.img'
                    'RPV.hdr'
                    'RPV.img'
                    'SPM.mat'
                    };
    % --------------------------------------------------------------------------
    % Clean-up files from whole image regression
    % --------------------------------------------------------------------------
    % Loop over sessions
    for iSessions = 1:size(IOI.fcIOS.SPM.fnameSPM, 1)
        % Loop over colors
        for iColors = 1:size(IOI.fcIOS.SPM.fnameSPM, 2)
            % Only remove SPM files if regression was succesful
            if IOI.fcIOS.SPM.wholeImageRegressOK{iSessions, iColors}
                % Loop over files to delete
                for iFiles = 1:numel(files2Delete)
                    fName = fullfile(IOI.fcIOS.SPM.fnameSPM{iSessions, iColors},...
                        files2Delete{iFiles});
                    % Undocumented MATLAB feature to delete files with Java
                    java.io.File(fName).delete();
                end
            end
        end
    end
    % --------------------------------------------------------------------------
    
    % --------------------------------------------------------------------------
    % Clean-up files from ROIs regression
    % --------------------------------------------------------------------------
    % Loop over ROIs
    for iROI = 1:numel(IOI.fcIOS.SPM.fnameROISPM)
        % Loop over sessions
        for iSessions = 1:size(IOI.fcIOS.SPM.fnameROISPM{iROI}, 1)
            % Loop over colors
            for iColors = 1:size(IOI.fcIOS.SPM.fnameROISPM{iROI}, 2)
                % Only remove SPM files if regression was succesful
                if IOI.fcIOS.SPM.ROIregressOK{iROI}{iSessions, iColors}
                    % Loop over files to delete
                    for iFiles = 1:numel(files2Delete)
                        fName = fullfile(IOI.fcIOS.SPM.fnameROISPM{iROI}{iSessions, iColors},...
                            files2Delete{iFiles});
                        % Undocumented MATLAB feature to delete files with Java
                        java.io.File(fName).delete();
                    end
                end
            end
        end
    end
    % --------------------------------------------------------------------------
    
    % Clean up succesful
    cleanupOK = true;
    fprintf('Cleanup GLM files for %s done! Elapsed time: %s\n',...
        IOI.subj_name, datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS'));
catch exception
    cleanupOK = false;
    disp(exception.identifier)
    disp(exception.stack(1))
end

% EOF
