function ioi_save_masked_images(image, fname, voxel_size,Im,tit,sgn,thold)
%also generate masked figures
[nx ny] = size(image);
image = image(:);
% if sgn > 0
%     image(image<thold) = 0;
% else
%     image(image>-thold) = 0;
% end
image(image>-thold & image<thold) = 0;
image(isnan(image)) = 0;
image = reshape(image,[nx ny]);
fname = [fname 'masked'];
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
