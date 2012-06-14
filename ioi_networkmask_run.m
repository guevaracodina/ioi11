function out = ioi_networkmask_run(job)
% Manual segmentation of the brain to provide a mask for resting-state
% functional connectivity mapping with optical intrinsic signal imaging (fcOIS)
% analysis. User should only select those pixels belonging to the brain.
% SYNTAX:
% out = ioi_networkmask_run(job)
% INPUTS:
% job
% OUTPUTS:
% out       
%_______________________________________________________________________
% Copyright (C) 2011 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%______________________________________________________________________

for SubjIdx=1:length(job.IOImat)
    try
        tic
        clear IOI
        
        %Load IOI.mat information
        IOImat = job.IOImat{SubjIdx};
        load(IOImat);
        if ~isfield(IOI, 'fcIOS')
            % Create fcIOS field to contain the whole structure of fcIOS
            % utilities
            IOI.fcIOS = struct([]);
            IOI.fcIOS.mask = struct([]);
            IOI.fcIOS.seeds =[];
            IOI.fcIOS.filtNdown = struct([]);
            IOI.fcIOS.filtNdown = struct([]);
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
            
            % Grayscale colormap for contrast enhancement
            cmap = contrast(im_anat);
            h_im = imagesc(im_anat); colormap(cmap);
            axis image
            title('Choose a mask containing only brain pixels','FontSize',13)
            
            % Start interactive ROI tool to choose brain mask
            h_ROI = impoly(gca);
            
            % Set color of ROI brain mask
            setColor(h_ROI, 'y')
            
            % Keep the polygon inside the original dimensions
            fcn = makeConstrainToRectFcn('impoly',get(gca,'XLim'), get(gca,'YLim'));
            setPositionConstraintFcn(h_ROI,fcn)
            
            % Create a mask
            BW_mask = createMask(h_ROI,h_im);
            
            % Display masked image on SPM graphics window
            imagesc(im_anat .* BW_mask);
            axis image
            title(['Mask for subject ' IOI.subj_name],'FontSize',13)
            
            % 3-Sauver en sortie une image nifti qui s'appelle brainmask.nii qui
            % vaut 1 dans le cerveau et 0 en dehors.
            
            % Create filename according the existing nomenclature at subject level
            maskName = [IOI.subj_name '_anat_brainmask.nii'];
            
            % Create and write a NIFTI file in the subject folder
            hdr = spm_vol(fullfile(dirName,fileName));
            ioi_create_vol(fullfile(dirName, maskName), ...
                hdr.dim, hdr.dt, hdr.pinfo, hdr.mat, hdr.n, BW_mask);
            
            % Identifier dans IOI le nom du fichier masque
            IOI.fcIOS.mask.fname = fullfile(dirName, maskName);
            
            % Mask created succesfully!
            IOI.fcIOS.mask.maskOK = true;
            save(IOImat,'IOI');
            toc
        end
        out.IOImat{SubjIdx} = IOImat;
        disp(['Subject ' int2str(SubjIdx) ' (' IOI.subj_name ')' ' complete']);
    catch exception
        disp(exception.identifier)
        disp(exception.stack(1))
        out.IOImat{SubjIdx} = job.IOImat{SubjIdx};
    end % End of try
end % End of main for
end % End of function

%_______________________________________________________________________
% Auxiliary functions
%_______________________________________________________________________
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
