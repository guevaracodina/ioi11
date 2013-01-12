%% Resolution of intrinsic imaging system

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


%% Get a 1x1 binning anatomical image
% 
% addpath(genpath('D:\Edgar\ssoct\Matlab'))
% export_fig(fullfile('D:\Edgar\Documents\Dropbox\Docs\Epilepsy\figs','LFP_4AP'),'-png',gcf)

