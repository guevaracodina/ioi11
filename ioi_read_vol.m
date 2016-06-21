function im_anat = ioi_read_vol(nifti_anat_file)
% Read anatomical image from nifti file
vol = spm_vol(nifti_anat_file);
im_anat = spm_read_vols(vol);
end % End ioi_read_vol