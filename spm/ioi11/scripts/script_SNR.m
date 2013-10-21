%% script_SNR
load('E:\Edgar\Data\IOS_Resolution\12_09_28,MT11\IOI.mat')
imagePath = 'D:\Edgar\Documents\Dropbox\Docs\Thesis\Figures\OIS_SNR';
imSize = [0.1 0.1 3 3];

%% NIfTI volumes
% Only 1 session
s1 = 1;
% IOI.color.eng = 'RGYL'
c1 = 1;
R = ioi_get_images(IOI,1:IOI.sess_res{s1}.n_frames,c1,s1,'E:\Edgar\Data\IOS_Resolution\12_09_28,MT11\',false);
c1 = 2;
G = ioi_get_images(IOI,1:IOI.sess_res{s1}.n_frames,c1,s1,'E:\Edgar\Data\IOS_Resolution\12_09_28,MT11\',false);
c1 = 3;
Y = ioi_get_images(IOI,1:IOI.sess_res{s1}.n_frames,c1,s1,'E:\Edgar\Data\IOS_Resolution\12_09_28,MT11\',false);
c1 = 4;
L = ioi_get_images(IOI,1:IOI.sess_res{s1}.n_frames,c1,s1,'E:\Edgar\Data\IOS_Resolution\12_09_28,MT11\',false);

%% SNR computation
h = figure();  set(gcf, 'color', 'w')
% Choose color IOI.color.eng = 'RGYL'
c1 = 3;
switch(c1)
    case 1
        OISimg = squeeze(R(:,:,randi(size(R,3), 1)));
        imagesc(OISimg)
    case 2
        OISimg = squeeze(G(:,:,randi(size(G,3), 1)));
        imagesc(OISimg)
    case 3
        OISimg = squeeze(Y(:,:,randi(size(Y,3), 1)));
        imagesc(OISimg)
    case 4
        OISimg = squeeze(L(:,:,randi(size(L,3), 1)));
        imagesc(OISimg)
end
colormap(gray(256))
set(gca,'XTick',[])
set(gca,'YTick',[])
axis image
% Rectangular ROI
title('Choose signal ROI')
roiPos = round(wait(imrect));
% ROI coordinates
y1 = roiPos(2);
y2 = roiPos(2) + roiPos(4);
x1 = roiPos(1);
x2 = roiPos(1) + roiPos(3);
% ROI signal
OISsignal = OISimg(y1:y2,x1:x2);
OISsignal = mean2(OISsignal);
title('Choose background ROI')
roiPos = round(wait(imrect));
% ROI coordinates
y1 = roiPos(2);
y2 = roiPos(2) + roiPos(4);
x1 = roiPos(1);
x2 = roiPos(1) + roiPos(3);
% ROI signal
OISnoise = OISimg(y1:y2,x1:x2);
OISnoise = std2(OISnoise);
title(sprintf('SNR(%c) = %0.2f dB',IOI.color.eng(c1), 20*log10(OISsignal / OISnoise)),'FontSize',12); 
fprintf('SNR(%c) = %f dB\n',IOI.color.eng(c1), 20*log10(OISsignal / OISnoise));

% Specify window units
set(h, 'units', 'inches')
% Change figure and paper size
set(h, 'Position', imSize)
set(h, 'PaperPosition', imSize)

%% Print images
% Save as PNG
print(h, '-dpng', fullfile(imagePath, sprintf('OIS_SNR_%c', IOI.color.eng(c1))), '-r300');
% Save as a figure
saveas(h, fullfile(imagePath, sprintf('OIS_SNR_%c', IOI.color.eng(c1))), 'fig');

%% Temporal SNR map
% Choose color IOI.color.eng = 'RGYL'
switch(c1)
    case 1
        OISmean = mean(R, 3);
        OISstd = std(R, 0, 3);
    case 2
        OISmean = mean(G, 3);
        OISstd = std(G, 0, 3);
    case 3
        OISmean = mean(Y, 3);
        OISstd = std(Y, 0, 3);
    case 4
        OISmean = mean(L, 3);
        OISstd = std(L, 0, 3);
end

% Temporal SNR
OIS_SNR = 20*log10(OISmean ./ OISstd);

%% Temporal SNR figure
h = figure;  set(gcf, 'color', 'w')
imagesc(OIS_SNR)
colormap(gray(256));
set(gca,'XTick',[])
set(gca,'YTick',[])
t = colorbar; 
set(get(t,'title'),'String', 'dB', 'FontSize', 12);
title(sprintf('SNR map(%c)',IOI.color.eng(c1)),'FontSize',12);
axis image

% Specify window units
set(h, 'units', 'inches')
% Change figure and paper size
set(h, 'Position', imSize)
set(h, 'PaperPosition', imSize)

%% Print images
% Save as PNG
print(h, '-dpng', fullfile(imagePath, sprintf('OIS_SNR_map_%c', IOI.color.eng(c1))), '-r300');
% Save as a figure
saveas(h, fullfile(imagePath, sprintf('OIS_SNR_map_%c', IOI.color.eng(c1))), 'fig');

% EOF
