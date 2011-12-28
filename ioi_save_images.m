function ioi_save_images(image, fname, voxel_size,Im,tit)
try 
    Im.save_figures;
catch
    Im.save_figures = 1;
    Im.generate_figures = 0;
end
ioi_save_nifti(image, fname, voxel_size)
h = figure;
%save also as figures
imagesc(image);
title(tit);
colorbar;
if Im.save_figures
    print(h, '-dtiffn', [fname '.tiff']);
    saveas(h,[fname '.fig']);
end
if ~Im.generate_figures, close(h); end
end
