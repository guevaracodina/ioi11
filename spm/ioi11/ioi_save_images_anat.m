function ioi_save_images_anat(image, fname, voxel_size,Im,tit,A)
anat = A.anat;
try
    mask = A.mask;
end
try
    Im.save_figures;
catch
    Im.save_figures = 1;
    Im.generate_figures = 0;
end
if ~isfield(Im,'cbar')
    Im.cbar.colorbar_override = 0;
end
% Colormap to enhance contrast
cmap = contrast(anat);

th_z = 1; %set low threshold
%must be a t-statistic
s_map = image;

min_T = min(s_map(index_over2));
max_T = max(s_map(index_over2));
if min_T == max_T %for example if only one data point in index_over
    if max_T > 0
        min_T = -min_T; % - min_max_gap; %careful -- need to include grey
    else
        max_T = -min_T;
    end
end

%picking lower half of jet colormap;
load(split);
% fcool = 64;
% fgray = 64;
% fhot = 64;
cmap_res = 3*64;
sbar = linspace(min_T, max_T, cmap_res); %for colormap in 1-192 range
%fraction of colorbar in each of cool, gray and hot areas:
fcool = round((-th_z-min_T)/(max_T-min_T)*cmap_res);
fgray = round(2*th_z/(max_T-min_T)*cmap_res);
fhot = round((max_T-th_z)/(max_T-min_T)*cmap_res);

%Choose location of background brain as a function of th_z:
%want it to be in between -th_z and +th_z
%To ensure no overlap between background brain and
%(de)activations, reduce th_z by a factor
th_z_eff = th_z*th_z_shrink;
maxb = max(brain(:)); %minb = min(brain(:)); %expect =0
T_brain_over = (2*th_z_eff/maxb).*brain -th_z_eff;
%add both negative and positive regions of interest
if ~isempty(index_over)
    T_brain_over(index_over) = s_map(index_over);
end
if combinedfig && ~isempty(index_over2)
    T_brain_over(index_over2) = s_map(index_over2);
end
jet2 = colormap(jet(2*fcool)); %doubling resolution of jet colormap            
split1 = [jet2(1:fcool,:);gray(fgray);hot(fhot)];
                    
                    
%ioi_save_nifti(image, fname, voxel_size)
h = figure;
%save also as figures
if Im.cbar.colorbar_override
    image(image>Im.cbar.c_max) = Im.cbar.c_max;
    image(image<Im.cbar.c_min) = Im.cbar.c_min;
    image(1,1) = Im.cbar.c_max;
    image(end,end) = Im.cbar.c_min;
    %hc = colorbar;
    %set(hc, 'YLim', [Im.cbar.c_min Im.cbar.c_max]);
end
imagesc(image);
colorbar;
title(tit);
if Im.save_figures
    print(h, '-dtiffn', [fname '.tiff']);
    saveas(h,[fname '.fig']);
end
if ~Im.generate_figures, close(h); end
end
