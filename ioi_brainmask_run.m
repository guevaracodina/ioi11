function [out, BW_mask] = ioi_brainmask_run(job)
% Manual segmentation of the brain to provide a mask for resting-state
% functional connectivity mapping with optical intrinsic signal imaging (fcOIS)
% analysis. User should only select those pixels belonging to the brain.
% SYNTAX:
% out = ioi_brainmask_run(job)
% INPUTS:
% job
% OUTPUTS:
% out       Structure containing the names of IOI matrices
% BW_mask   Binary mask that contains brain pixels
%_______________________________________________________________________________
% Copyright (C) 2012 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________

for SubjIdx=1:length(job.IOImat)
    try
        tic
        clear IOI
        % Load IOI.mat information
        IOImat = job.IOImat{SubjIdx};
        load(IOImat);
        
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
            BW_mask = ioi_roi_spline(im_anat);
            % ------------------------------------------------------------------
            
            % Display masked image on SPM graphics window
            minVal = min(im_anat(:));
            maxVal = max(im_anat(:));
            imagesc(im_anat .* BW_mask, [minVal maxVal]);
            axis image
            title(['Mask for subject ' IOI.subj_name],'FontSize',13)
            
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
end % End of main for
end % End of function

%_______________________________________________________________________________
% Auxiliary functions
%_______________________________________________________________________________
function hdr = ioi_create_vol(fname, dim, dt, pinfo, mat, n, data)
hdr = struct('fname',fname,...
    'dim', dim,...
    'dt',   dt,...
    'pinfo',pinfo,...
    'mat',  mat,...
    'n', n);
hdr = spm_create_vol(hdr);
spm_write_vol(hdr, data);
end % End ioi_create_vol

function im_anat = ioi_read_vol(nifti_anat_file)
% Read anatomical image from nifti file
vol = spm_vol(nifti_anat_file);
im_anat = spm_read_vols(vol);
end % End ioi_read_vol
