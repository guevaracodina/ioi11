% No check on dimensions since it is only called internally,
% Returns a temporal slice of a volume time series
function slice=ioi_read_time_vol(V,n)

spm_check_orientations(V);
%-Read in image data
%--------------------------------------------------------------------------
slice = zeros([V(1).dim(1:3)]);       %-image data matrix

for p=1:V(1).dim(3)
    slice(:,:,p) = spm_slice_vol(V(n),spm_matrix([0 0 p]),V(n).dim(1:2),0);
end
