% Save the data in nifti format
function ioi_save_nifti(image, fname, voxel_size)
% At this point we have our images, need to transform in nifti for further
% SPM processing
% Should we hardcode all of this or should it be the user's choice?
origin=[0 0 0];
datatype=16; % 'single'
%datatype = 4; %int16
description='IOI Image, created with ioi_msioi_run.';
ioi_write_nifti(image,voxel_size,origin,datatype,description,fname);
