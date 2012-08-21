function ioi_save_images(image, fname, voxel_size,Z,tit)
%save nifti first
ioi_save_nifti(image, [fname '.nii'], voxel_size)
if ~isfield(Z,'save_figures')
    Z.save_figures = 1;
    Z.generate_figures = 0;
end
if ~isfield(Z,'cbar')
    Z.cbar.colorbar_override = 0;
end
if ~isfield(Z,'do_superpose')
    Z.do_superpose = 0;
end
h = figure;
try
    %include ROIs
    if isfield(Z,'superpose_ROIs') && Z.superpose_ROIs && Z.do_superpose
        max_T = max(image(:));
        for r1=1:length(Z.ROI)
            Vr{r1} = spm_vol(Z.ROI{r1}.fname);
            Yr{r1} = spm_read_vols(Vr{r1}); %this could be stored instead of reread each time
            if r1 == 1
                [nx ny] = size(Yr{r1});
                image = imresize(image,[nx ny]);
            end
            %extract border of ROI
            a0 = [diff(Yr{r1},1,1); zeros(1,ny)];
            if isempty(a0(:)>0) || isempty(a0(:)<0)
                %make line wider
                a0 = [diff(Yr{r1},1,1); zeros(1,ny)] + [zeros(1,ny); diff(Yr{r1},1,1)];
            end
            b0 = [diff(Yr{r1},1,2) zeros(nx,1)];
            if isempty(b0(:)>0) || isempty(b0(:)<0)
                %make line wider
                b0 = [diff(Yr{r1},1,2) zeros(nx,1)] + [zeros(nx,1) diff(Yr{r1},1,2)];
            end
            tmp = a0 + b0;
            %add to the image
            image(tmp ~= 0) = max_T;
            %image(tmp>0) = max_T;
        end
    end
catch
    disp('Unable to superpose ROIs on image');
end

try
    if isfield(Z,'superpose_anatomical') && Z.superpose_anatomical && Z.do_superpose
        adj = 1.01;
        prec = 3;
        min_T = min(image(:));
        max_T = max(image(:));
        cmap = 64*prec;
        thz = Z.superpose_threshold;
        if min_T > -thz, min_T = -adj*thz; end
        if max_T < thz, max_T = adj*thz; end
        tick_number = 7;
        fontsize_choice = 10;
        V = spm_vol(Z.file_anat);
        A = spm_read_vols(V); %this could be stored instead of reread each time
        %rescale values of anatomical image
        mA = min(A(:));
        MA = max(A(:));
        A = (max_T-min_T)/(6*adj)*(2*A-mA-MA)/(MA-mA);
        [nx ny] = size(A);
        sbar1 = linspace(min_T, -thz, cmap);
        %sbar2 = linspace(-thz, thz, 64);
        sbar3 = linspace(thz, max_T, cmap);
        
        B = imresize(image,[nx ny]);
        index_over = B > thz;
        index_under = B < -thz;
        A(index_over) = B(index_over);
        A(index_under) = B(index_under);
        image = A;
        %figure; 
        imagesc(image);
        axis off
        load split
        fcool = round((-thz-min_T)/(max_T-min_T)*cmap);
        fgray = round(2*thz/(max_T-min_T)*cmap);
        fhot = round((max_T-thz)/(max_T-min_T)*cmap);
        %fcool = cmap; fgray = cmap; fhot = cmap;
        jet2 = colormap(jet(2*fcool)); %doubling resolution of jet colormap,
        %then remove colors already used in colormap hot:
        split1 = [jet2(1:fcool,:);gray(fgray);hot(fhot)];
        colormap(split1);
        hc1 = colorbar('EastOutside');
        hc2 = colorbar('WestOutside');
        hc1_min = sbar3(1);
        hc1_max = sbar3(cmap);
        hc2_min = sbar1(1);
        hc2_max = sbar1(cmap);
        hc1 = ioi_set_colorbar_limits_fontsize(hc1,hc1_min,hc1_max,tick_number,fontsize_choice);
        hc2 = ioi_set_colorbar_limits_fontsize(hc2,hc2_min,hc2_max,tick_number,fontsize_choice);
    else
        imagesc(image); colorbar;
    end
catch
    disp('Unable to superpose anatomical image');
    imagesc(image); colorbar;
end

% %save also as figures
if Z.cbar.colorbar_override
    image(image>Z.cbar.c_max) = Z.cbar.c_max;
    image(image<Z.cbar.c_min) = Z.cbar.c_min;
    image(1,1) = Z.cbar.c_max;
    image(end,end) = Z.cbar.c_min;
    %hc = colorbar;
    %set(hc, 'YLim', [Im.cbar.c_min Im.cbar.c_max]);
end

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