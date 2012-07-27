function [IOI mask] = ioi_get_brain_mask(IOI,job)
% Extracts time series for the brain mask used in functional connectivity
% mapping with Intrinsic Optical Signal (fcIOS)
%_______________________________________________________________________________
% Copyright (C) 2012 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________

[all_sessions selected_sessions] = ioi_get_sessions(job);
% [all_ROIs selected_ROIs] = ioi_get_ROIs(job);
%loop over ROIs
% first_pass = 1;
% for r1=1:length(IOI.res.ROI)
%     if all_ROIs || sum(r1==selected_ROIs)
r1=1; % Only one brain mask per subject
vol = spm_vol(IOI.fcIOS.mask.fname);
tmp_mask = logical(spm_read_vols(vol));

%shrink mask to voxel size
if IOI.res.shrinkageOn
    sz = size(tmp_mask);
    %careful, this floor might not lead to the
    %correct size - better to do a check on image
    %size
    mask{r1} = imresize(tmp_mask,[floor(sz(1)/IOI.res.shrink_x) floor(sz(2)/IOI.res.shrink_y)],'bicubic');
else
    mask{r1} = tmp_mask;
end
%         if first_pass

%check size of mask is OK
if all_sessions %first available session
    s1 = 1;
else
    s1 = selected_sessions(1);
end

%find non-deleted file
foundfile = 0;

for c0=1:length(IOI.sess_res{s1}.fname)
    if ~foundfile
        % Check if laser is recorded
        if ~isempty(IOI.sess_res{s1}.fname{c0})
            fname_list = IOI.sess_res{s1}.fname{c0};
            fname = fname_list{1}; %1st file in list
            try
                vols = spm_vol(fname);
                d = spm_read_vols(vols);
                foundfile = 1;
            end
        end
    end
end
[d1 d2 d3 d4] = size(d);
%             first_pass = 0;
%         end
if ~(d1 == size(mask{1},1) && d2 == size(mask{1},2))
    mask{r1} = imresize(mask{r1},[d1 d2]);
end
%     end
% end
