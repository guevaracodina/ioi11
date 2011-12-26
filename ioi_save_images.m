function ioi_save_images(image, fname, voxel_size,Im,tit)
ioi_save_nifti(image, fname, voxel_size)
h = figure;
%save also as figures
imagesc(image);
title(tit);
colorbar;
if Im.save_figures
    print(h, '-dtiffn', [fname '.tif']);
    saveas(h,[fname '.fig']);
end
if ~Im.generate_figures, close(h); end
end
