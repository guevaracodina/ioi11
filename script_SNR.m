%% script_SNR
load('E:\Edgar\Data\IOS_Resolution\12_09_28,MT11\IOI.mat')

%% NIfTI volumes
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
c1 = 2;
switch(c1)
    case 1
        imagesc(R)
    case 2
        Gimg = squeeze(G(:,:,400));
        imagesc(Gimg)
    case 3
        imagesc(Y)
    case 4
        imagesc(L)
end
colormap(jet(256))
set(gca,'XTick',[])
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
Gsignal = Gimg(y1:y2,x1:x2);
Gsignal = mean2(Gsignal);
title('Choose background ROI')
roiPos = round(wait(imrect));
% ROI coordinates
y1 = roiPos(2);
y2 = roiPos(2) + roiPos(4);
x1 = roiPos(1);
x2 = roiPos(1) + roiPos(3);
% ROI signal
Gnoise = Gimg(y1:y2,x1:x2);
Gnoise = std2(Gnoise);
title(sprintf('SNR(G) = %0.2f dB',20*log10(Gsignal / Gnoise)),'FontSize',12); 
fprintf('SNR(%c) = %f dB\n',IOI.color.eng(c1), 20*log10(Gsignal / Gnoise));

% Specify window units
set(h, 'units', 'inches')
% Change figure and paper size
set(h, 'Position', [0.1 0.1 3 3])
set(h, 'PaperPosition', [0.1 0.1 3 3])

%% Print images
% Save as PNG
print(h, '-dpng', fullfile('D:\Edgar\Documents\Dropbox\Docs\Thesis\Figures\OIS_SNR', 'OIS_SNR'), '-r300');
% Save as a figure
saveas(h, fullfile('D:\Edgar\Documents\Dropbox\Docs\Thesis\Figures\OIS_SNR', 'OIS_SNR'), 'fig');


%% Rectangular ROI (OIS)
title('Choose signal ROI')
roiPos = round(wait(imrect));
% ROI coordinates
y1 = roiPos(2);
y2 = roiPos(2) + roiPos(4);
x1 = roiPos(1);
x2 = roiPos(1) + roiPos(3);
% ROI signal
OISsignal = OISimg(y1:y2,x1:x2);
OISsignal = mean(OISsignal(:));
title('Choose background ROI')
roiPos = round(wait(imrect));
% ROI coordinates
y1 = roiPos(2);
y2 = roiPos(2) + roiPos(4);
x1 = roiPos(1);
x2 = roiPos(1) + roiPos(3);
% ROI signal
OISnoise = OISimg(y1:y2,x1:x2);
OISnoise = std(OISnoise(:));
title('')
fprintf('SNR(OIS) = %f dB\n',20*log10(OISsignal / OISnoise));

%% SNR map
HbTmean = mean(squeeze(HbT(:,:,1,:)),3);
HbTstd = std(squeeze(HbT(:,:,1,:)), 0, 3);
HbTSNR = 20*log10(HbTmean ./ HbTstd);

SO2mean = mean(squeeze(SO2(:,:,1,:)),3);
SO2std = std(squeeze(SO2(:,:,1,:)), 0, 3);
SO2SNR = 20*log10(SO2mean ./ SO2std);

h = figure;  set(gcf, 'color', 'w')
subplot(121)
imagesc(PAT.PAparam.WidthAxis, PAT.PAparam.DepthAxis, HbTSNR)
colormap(gray); colorbar
xlabel('[mm]','FontSize',12); ylabel('[mm]','FontSize',12); 
title('HbT','FontSize',12);
axis image

subplot(122)
imagesc(PAT.PAparam.WidthAxis, PAT.PAparam.DepthAxis, SO2SNR)
colormap(gray); colorbar
xlabel('[mm]','FontSize',12); ylabel('[mm]','FontSize',12); 
title('SO_2','FontSize',12);
axis image

% Specify window units
set(h, 'units', 'inches')
% Change figure and paper size
set(h, 'Position', [0.1 0.1 6 3])
set(h, 'PaperPosition', [0.1 0.1 6 3])

%% Print images
% Save as PNG
print(h, '-dpng', fullfile('D:\Edgar\Documents\Dropbox\Docs\Thesis\Figures\OIS_SNR', 'OIS_SNR_map'), '-r300');
% Save as a figure
saveas(h, fullfile('D:\Edgar\Documents\Dropbox\Docs\Thesis\Figures\OIS_SNR', 'OIS_SNR_map'), 'fig');


% EOF
