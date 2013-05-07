%% Resolution of intrinsic imaging system
% Using R3L3S1N - Negative 1951 USAF Test Target, 3" x 3" 
% http://www.thorlabs.com/catalogpages/V21/1786.PDF

load('E:\Edgar\Data\IOS_Resolution\12_09_28,MT10\S01\all_images.mat')

%% Plot mean intensities
avg_Int = zeros([numel(images) 1]);
for i=1:numel(images)
    avg_Int(i) = mean(images{i}(:));
end

%% Plot intensities
figure; set(gcf,'color','w')
plot(avg_Int,'k-o')

%% Separate images acording to intensities
Rmean = mean(images{23}(:));
Gmean = mean(images{24}(:));
Ymean = mean(images{25}(:));
Lmean = mean(images{26}(:));
imBound = 5;
clear Rimages Gimages Yimages Limages
Rcount = 1; Gcount = 1; Ycount = 1; Lcount = 1;
for i=1:numel(images)
    avg_Int = mean(images{i}(:));
    if (avg_Int > Rmean-imBound) && (avg_Int < Rmean+imBound)
        Rimages(:,:,Rcount) = images{i};
        Rcount = Rcount + 1;
    elseif (avg_Int > Gmean-imBound) && (avg_Int < Gmean+imBound)
        Gimages(:,:,Gcount) = images{i};
        Gcount = Gcount + 1;
    elseif (avg_Int > Ymean-imBound) && (avg_Int < Ymean+imBound)
        Yimages(:,:,Ycount) = images{i};
        Ycount = Ycount + 1;
    else
        Limages(:,:,Lcount) = images{i};
        Lcount = Lcount + 1;
    end
end

%% Average images
Ravg = mean(Rimages,3);
Gavg = mean(Rimages,3);
Yavg = mean(Rimages,3);
Lavg = mean(Rimages,3);

%% Plot single images
imLims = [0 4095];
figure; set(gcf,'color','w')
colormap (gray(255))

subplot(221); imagesc(images{23},imLims); axis image; 
set(gca,'Xtick',[]); set(gca,'yTick',[]); colorbar; 
title('Red','FontSize',14);

subplot(222); imagesc(images{24},imLims); axis image;
set(gca,'Xtick',[]); set(gca,'yTick',[]); colorbar; 
title('Green','FontSize',14);

subplot(223); imagesc(images{25},imLims); axis image; 
set(gca,'Xtick',[]); set(gca,'yTick',[]); colorbar; 
title('Yellow','FontSize',14); 

subplot(224); imagesc(images{26},imLims); axis image;
set(gca,'Xtick',[]); set(gca,'yTick',[]); colorbar; 
title('Laser','FontSize',14); 


%% Plot average images
imLims = [0 4095];
figure; set(gcf,'color','w')
colormap (gray(255))

subplot(221); imagesc(Ravg,imLims); axis image; 
set(gca,'Xtick',[]); set(gca,'yTick',[]); colorbar; 
title('<Red>','FontSize',14);

subplot(222); imagesc(Gavg,imLims); axis image;
set(gca,'Xtick',[]); set(gca,'yTick',[]); colorbar; 
title('<Green>','FontSize',14);

subplot(223); imagesc(Yavg,imLims); axis image; 
set(gca,'Xtick',[]); set(gca,'yTick',[]); colorbar; 
title('<Yellow>','FontSize',14); 

subplot(224); imagesc(Lavg,imLims); axis image;
set(gca,'Xtick',[]); set(gca,'yTick',[]); colorbar; 
title('<Laser>','FontSize',14); 

%% Plot single image
close all
resTarget = Gavg(45:755,85:795);
imLims = [118 1217];
h1 = figure; set(gcf,'color','w')
colormap (gray(255))
imagesc(resTarget, imLims); axis image; 
set(gca,'Xtick',[]); set(gca,'yTick',[]);  

% Zoom on smaller elements
resTargetZoom = resTarget(284:468, 322:506);
h2 = figure; set(gcf,'color','w')
colormap (gray(255))
imagesc(resTargetZoom, imLims); axis image; 
set(gca,'Xtick',[]); set(gca,'yTick',[]); 

figure(h2)
% Profile across Group 2, elements 1:4 (4, 4.49, 5.04, 5.66 line pairs/mm)
xi = [163.5921; 163.5921];
yi = [29.4401;  101.3845];
[cx,cy,resTargetProfile,xi,yi] = improfile(resTargetZoom, xi, yi);

resTargetProfile = resTargetProfile - min(resTargetProfile(:));
resTargetProfile = resTargetProfile ./ max(resTargetProfile);
% Group 0, element 1 (1 line pair/mm) => 35.5 - 18 pixels = 500 um
umPerPixel = 500/(35.5-18);
distance = (0:numel(resTargetProfile)-1)*umPerPixel;

h3 = figure; set(gcf,'color','w')
plot(distance, resTargetProfile, 'k-', 'LineWidth', 2)
axis tight
set(gca, 'FontSize',10)
set(gca, 'YTick', [0 0.5 1])
xlabel('length (\mum)','FontSize',10); 
ylabel('Intensity [a.u.]','FontSize',10); 
res = 1/(5.66*2); % FWHM ~70um

%% Save images
addpath(genpath('D:\Edgar\ssoct\Matlab'))
figSize = [6.5/3 6.5/3];

figure(h1);
% Specify window units
set(h1, 'units', 'inches')
% Change figure and paper size
set(h1, 'Position', [0.1 0.1 figSize])
set(h1, 'PaperPosition', [0.1 0.1 figSize])
export_fig(fullfile('D:\Edgar\Documents\Dropbox\Docs\Thesis\Figures\OIS_resolution','OIS_resTarget'),'-png',h1)

figure(h2);
% Specify window units
set(h2, 'units', 'inches')
% Change figure and paper size
set(h2, 'Position', [0.1 0.1 figSize])
set(h2, 'PaperPosition', [0.1 0.1 figSize])
export_fig(fullfile('D:\Edgar\Documents\Dropbox\Docs\Thesis\Figures\OIS_resolution','OIS_resTargetZoom'),'-png',h2)

figure(h3);
% Specify window units
set(h3, 'units', 'inches')
% Change figure and paper size
set(h3, 'Position', [0.1 0.1 figSize])
set(h3, 'PaperPosition', [0.1 0.1 figSize])
export_fig(fullfile('D:\Edgar\Documents\Dropbox\Docs\Thesis\Figures\OIS_resolution','OIS_resTargetProfile'),'-png',h3)

