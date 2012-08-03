function ioi_save_images(image, fname, voxel_size,Z,tit)
try
    Z.save_figures;
catch
    Z.save_figures = 1;
    Z.generate_figures = 0;
end
if ~isfield(Z,'cbar')
    Z.cbar.colorbar_override = 0;
end
if Z.superpose_anatomical
    V = spm_vol(Z.file_anat); 
    Y = spm_read_vols(V); %this could be stored instead of reread each time
    %rescale values of anatomical image
    load split
    fcool = 64;
fgray = 64;
fhot = 64;
    thz = 1;
    %set a threshold
    index_over = image > thz;
    index_under = image < -thz;
    min_T = min(image);
    max_T = max(image);
    smin_T = max_T - ((max_T - min_T)./63) * 127;
    sbar = linspace(smin_T, max_T, 128);
    th_image = ((-sbar(1) + sbar(64))/(0.5)).*image + sbar(1);
    th_image(index_over) = image(index_over);
    th_image(index_under) = image(index_under);
end
%include ROIs
if Z.superpose_ROIs
    for r1=1:length(Z.ROIname)
        Vr{r1} = spm_vol(Z.ROIname{r1}); 
        Yr{r1} = spm_read_vols(Vr{r1}); %this could be stored instead of reread each time
        %extract border of ROI
        tmp = diff(Yr{r1},1,1) + diff(Yr{r1},1,2);
        %add to the image
        image(tmp>0) = 1;
    end
end
    
ioi_save_nifti(image, [fname '.nii'], voxel_size)
h = figure;
%save also as figures
if Z.cbar.colorbar_override
    image(image>Z.cbar.c_max) = Z.cbar.c_max;
    image(image<Z.cbar.c_min) = Z.cbar.c_min;
    image(1,1) = Z.cbar.c_max;
    image(end,end) = Z.cbar.c_min;
    %hc = colorbar;
    %set(hc, 'YLim', [Im.cbar.c_min Im.cbar.c_max]);
end
imagesc(image);
colorbar;
% Here I prevent title to display characters preceded by _ as subscripts
title(tit, 'interpreter', 'none');
if Z.save_figures
    %print(h, '-dtiffn', [fname '.tiff']);
    print(h, '-dpng', [fname '.png'], '-r300');
    saveas(h,[fname '.fig']);
    % Save as PNG: ~10x smaller file size and 2x the resolution //EGC
    % It is necessary to add file extension //EGC

end
if ~Z.generate_figures, close(h); end