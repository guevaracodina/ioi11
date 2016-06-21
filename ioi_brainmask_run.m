function [out, BW_mask] = ioi_brainmask_run(job)
% Manual segmentation of the brain to provide a mask for resting-state
% functional connectivity mapping with optical intrinsic signal imaging (fcOIS)
% analysis. User should only select those pixels belonging to the brain.
% SYNTAX:
% out = ioi_brainmask_run(job)
% INPUTS:
% job       Matlab batch job structure
% OUTPUTS:
% out       Structure containing the names of IOI matrices
% BW_mask   Binary mask that contains brain pixels
%_______________________________________________________________________________
% Copyright (C) 2012 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________

try
    SelectPreviousMask = job.SelectPreviousMask;
catch
    SelectPreviousMask = false;
end

% Loop over subjects
for SubjIdx = 1:length(job.IOImat)
    try
        tic
        clear IOI
        % Load IOI.mat information
        [IOI IOImat dir_ioimat] = ioi_get_IOI(job,SubjIdx);
        % Load IOI.mat information
        % IOImat = job.IOImat{SubjIdx};
        % load(IOImat);
        
        if ~isfield(IOI, 'fcIOS')
            % Create fcIOS field to contain the whole structure of fcIOS
            % utilities
            IOI.fcIOS = struct([]);
            IOI.fcIOS(1).mask = struct([]);
            IOI.fcIOS(1).LPF  = struct([]);
            IOI.fcIOS(1).filtNdown = struct([]);
            IOI.fcIOS(1).SPM = struct([]);
            IOI.fcIOS(1).corr = struct([]);
        end
        if ~isfield(IOI.fcIOS.mask,'maskOK') || job.force_redo
            
            if SelectPreviousMask
                % Choose a previous brain mask.
                % Goto interactive window
                h2 = spm_figure('GetWin', 'Interactive');
                spm_input(['Subject ' int2str(SubjIdx) 'of ' int2str(length(job.IOImat)) ' (' IOI.subj_name ')' ],'-1','d');
                SelPrevMask = spm_input('Select a previous brain mask?','+1','y/n');
                if SelPrevMask == 'y'
                    [tPrevMask stsPrevMask] = spm_select(1,'mat','Select IOI.mat structure containing information on desired brain mask','',dir_ioimat,'IOI.mat',1);
                    if stsPrevMask
                        IOI0 = IOI; %Store current IOI
                        try
                            load(tPrevMask);
                            IOI_withMask = IOI;
                            IOI = IOI0;
                            try
                                IOI.fcIOS.mask(1).fname = IOI_withMask.fcIOS.mask.fname;
                                clear IOI_withMask
                            catch
                                IOI = IOI0;
                                disp('Specified IOI.mat structure does not contain valid brain mask information')
                            end
                        catch
                            disp('Could not load IOI.mat structure containing desired brain mask')
                        end
                    end
                end
            else % Choose manually a new brain mask. 
                % 1-Lis une image nifti de type anatomique (dans notre cas c'est l'image
                % verte de l'IOI) qui n'est qu'une tranche (sans temps)
                [dirName fileName fileExt] = fileparts(IOI.res.file_anat);
                fileName = strcat(fileName, fileExt);
                
                % Read anatomical NIFTI file
                im_anat = ioi_read_vol(fullfile(dirName, fileName));
                
                % 2-Afficher et utiliser imroi pour que l'usager puisse manuellement
                % identifier la zone d'intérêt contenant le cerveau (on click un
                % polygone)
                
                % Display anatomical image on SPM graphics window
                spm_figure('GetWin', 'Graphics');
                spm_figure('Clear', 'Graphics');
                
                % Start interactive ROI tool to choose spline brain mask
                % ------------------------------------------------------------------
                % BW_mask = ioi_roi_spline(im_anat);
                % ioi_roi_spline disabled beause of unactive spline toolbox // EGC
                imagesc(im_anat);
                axis image
                colormap gray
                title('Choose a mask containing only brain pixels','FontSize',13)
                BW_mask = roipoly;
                % ------------------------------------------------------------------
                
                % Display masked image on SPM graphics window
                minVal = min(im_anat(:));
                maxVal = max(im_anat(:));
                imagesc(im_anat .* BW_mask, [minVal maxVal]);
                axis image
                title(sprintf('Mask for subject %s; %d of %d', IOI.subj_name, SubjIdx, length(job.IOImat)),'FontSize',13)
                
                % 3-Sauver en sortie une image nifti qui s'appelle brainmask.nii qui
                % vaut 1 dans le cerveau et 0 en dehors.
                
                % Create filename according the existing nomenclature at subject level
                brainMaskName = [IOI.subj_name '_anat_brainmask.nii'];
                
                % Create and write a NIFTI file in the subject folder
                hdr = spm_vol(fullfile(dirName,fileName));
                ioi_create_vol(fullfile(dirName, brainMaskName), ...
                    hdr.dim, hdr.dt, hdr.pinfo, hdr.mat, hdr.n, BW_mask);
                if isempty(IOI.fcIOS.mask)
                    % To avoid emptyDotAssignment we create a field
                    IOI.fcIOS.mask = struct('fname', []);
                end
                
                % Identifier dans IOI le nom du fichier masque
                IOI.fcIOS.mask.fname = fullfile(dirName, brainMaskName);
            end % SelectPreviousMask
            % Mask created succesfully!
            IOI.fcIOS.mask.maskOK = true;
            save(IOImat,'IOI');
        end
        out.IOImat{SubjIdx} = IOImat;
        disp(['Elapsed time: ' datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS')]);
        disp(['Subject ' int2str(SubjIdx) ' (' IOI.subj_name ')' ' complete']);
    catch exception
        disp(exception.identifier)
        disp(exception.stack(1))
        out.IOImat{SubjIdx} = job.IOImat{SubjIdx};
    end % End of try
end % End of main for (loop over subjects)
end % End of function

%_______________________________________________________________________________
% Auxiliary functions
% EGC --> Moved to independent files, so they can be used outside
% ioi_brainmask_run
%_______________________________________________________________________________
% function hdr = ioi_create_vol(fname, dim, dt, pinfo, mat, n, data)
% hdr = struct('fname',fname,...
%     'dim', dim,...
%     'dt',   dt,...
%     'pinfo',pinfo,...
%     'mat',  mat,...
%     'n', n);
% hdr = spm_create_vol(hdr);
% spm_write_vol(hdr, data);
% end % End ioi_create_vol
% 
% function im_anat = ioi_read_vol(nifti_anat_file)
% % Read anatomical image from nifti file
% vol = spm_vol(nifti_anat_file);
% im_anat = spm_read_vols(vol);
% end % End ioi_read_vol

% EOF
