%% Load data
fprintf('Loading data...\n');
ioimat = 'D:\Edgar\Data\IOS_Res\12_06_20,CO04\ROIseeds\Series\FiltNDown\GLMfcIOS\IOI.mat';
load(ioimat);
load(IOI.fcIOS.SPM.fnameROIregress)
colorNames = fieldnames(IOI.color);

for s1=1,
    figure;
    set(gcf,'color','w')
    for c1=5:7,
        for r1=1:2,
            tic
            vol = spm_vol(IOI.fcIOS.SPM.fname{s1, c1});
            y = spm_read_vols(vol);
            ROIvol = spm_vol(IOI.fcIOS.SPM.fnameROInifti{r1}{s1, c1});
            ROI = spm_read_vols(ROIvol);
            maskVol = spm_vol(IOI.fcIOS.mask.fname);
            mask = spm_read_vols(maskVol);
            ROImaskVol = spm_vol(IOI.res.ROI{r1}.fname);
            ROImask = spm_read_vols(ROImaskVol);
            % Preallocate
            tempCorrMap = zeros([size(y,1) size(y,2)]);
            fprintf('Data loaded, %0.1f sec. elapsed\n',toc);
            % Find Pearson's correlation coefficient
            tic
            fprintf('Computing correlation map...\n');
            for iX = 1:size(y,1),
                for iY = 1:size(y,2),
                    tempCorrMap(iX, iY) = corr(squeeze(ROI), squeeze(y(iX, iY, 1, :)));
                end
            end
            corrMap{r1}{s1, c1} = tempCorrMap;
            fprintf('Pearson''s correlation coefficient computed. Seed %d S%d C%d (%s) %0.1f sec. elapsed\n',r1,s1,c1,colorNames{1+c1},toc);
        end % ROI/seeds loop
    end % Colors loop
end % Sessions loop



%% Display results
figure;
s1 = 1;
r1 = 2;
c1 = 7;
set(gcf,'Color','w')
subplot(221); imagesc(corrMap{r1}{s1, c1})
axis image
colorbar
title(sprintf('Correlation Map. ROI %d Session %d Color %s',r1,s1,IOI.color.eng(c1)))

ROImaskVol = spm_vol(IOI.res.ROI{r1}.fname);
ROImask = spm_read_vols(ROImaskVol);
subplot(222); imagesc(ROImask)
axis image
title('ROI mask')

subplot(223); imagesc(mask)
axis image
title('Brain mask')
