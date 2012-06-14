function ioi_save_images(image, fname, voxel_size,Im,tit)
try
    Im.save_figures;
catch
    Im.save_figures = 1;
    Im.generate_figures = 0;
end
if ~isfield(Im,'cbar')
    Im.cbar.colorbar_override = 0;
end
ioi_save_nifti(image, fname, voxel_size)
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
