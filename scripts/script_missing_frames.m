%% Load data
BT01 = load('F:\Edgar\Data\IOS_Test_Results\12_08_07,BT01\S01\missing_frames.mat');
BT02 = load('F:\Edgar\Data\IOS_Test_Results\12_08_07,BT02\S01\missing_frames.mat');
BT03 = load('F:\Edgar\Data\IOS_Test_Results\12_08_07,BT03\S01\missing_frames.mat');
BT01IOI = load('F:\Edgar\Data\IOS_Test_Results\12_08_07,BT01\IOI.mat');
BT02IOI = load('F:\Edgar\Data\IOS_Test_Results\12_08_07,BT02\IOI.mat');
BT03IOI = load('F:\Edgar\Data\IOS_Test_Results\12_08_07,BT03\IOI.mat');

%% Plot missing frames stats

% ------------------------------------------------------------------------------
BT01_total_frames = zeros([BT01.nColors * BT01.n_frames 1]);
BT01t = 0:BT01IOI.IOI.dev.TR:BT01IOI.IOI.dev.TR*(numel(BT01_total_frames)-1);
BT01_total_frames(BT01.missing_frames) = 1;
frameRange = 1e3:2e3;
fprintf('Subject %s %d missing frames out of %d frames (%0.2f%%)\n',...
                        BT01IOI.IOI.subj_name,...
                        numel(BT01.missing_frames), BT01.n_frames*BT01.nColors, ...
                        100*numel(BT01.missing_frames)/(BT01.n_frames*BT01.nColors));
figure; set(gcf,'color','w')
subplot(312)
imagesc(BT01t(frameRange), 1, BT01_total_frames(frameRange)', [0 1]);
colormap(flipud(gray))
% stem(BT01t(frameRange), BT01_total_frames(frameRange),'k-','MarkerSize',1);
title(sprintf('%s 4x4bin miss %d/%d(%0.2f%%)',...
                        BT01IOI.IOI.subj_name,...
                        numel(BT01.missing_frames), BT01.n_frames*BT01.nColors, ...
                        100*numel(BT01.missing_frames)/(BT01.n_frames*BT01.nColors)),...
                        'interpreter','none','FontSize',14);
xlabel('t [s]','FontSize',14)
set(gca,'yTick',0)

% ------------------------------------------------------------------------------
BT02_total_frames = zeros([BT02.nColors * BT02.n_frames 1]);
BT02t = 0:BT02IOI.IOI.dev.TR:BT02IOI.IOI.dev.TR*(numel(BT02_total_frames)-1);
BT02_total_frames(BT02.missing_frames) = 1;
frameRange = 1e3:2e3;
fprintf('Subject %s %d missing frames out of %d frames (%0.2f%%)\n',...
                        BT02IOI.IOI.subj_name,...
                        numel(BT02.missing_frames), BT02.n_frames*BT02.nColors, ...
                        100*numel(BT02.missing_frames)/(BT02.n_frames*BT02.nColors));
subplot(311)
imagesc(BT02t(frameRange), 1, BT02_total_frames(frameRange)', [0 1]);
colormap(flipud(gray))
% stem(BT02t(frameRange), BT02_total_frames(frameRange),'k-','MarkerSize',1);
title(sprintf('%s 8x8bin miss %d/%d(%0.2f%%)',...
                        BT02IOI.IOI.subj_name,...
                        numel(BT02.missing_frames), BT02.n_frames*BT02.nColors, ...
                        100*numel(BT02.missing_frames)/(BT02.n_frames*BT02.nColors)),...
                        'interpreter','none','FontSize',14);
xlabel('t [s]','FontSize',14)
set(gca,'yTick',0)

% ------------------------------------------------------------------------------
BT03_total_frames = zeros([BT03.nColors * BT03.n_frames 1]);
BT03t = 0:BT03IOI.IOI.dev.TR:BT03IOI.IOI.dev.TR*(numel(BT03_total_frames)-1);
BT03_total_frames(BT03.missing_frames) = 1;
frameRange = 1e3:2e3;
fprintf('Subject %s %d missing frames out of %d frames (%0.2f%%)\n',...
                        BT03IOI.IOI.subj_name,...
                        numel(BT03.missing_frames), BT03.n_frames*BT03.nColors, ...
                        100*numel(BT03.missing_frames)/(BT03.n_frames*BT03.nColors));
subplot(313)
imagesc(BT03t(frameRange), 1, BT03_total_frames(frameRange)', [0 1]);
colormap(flipud(gray))
% stem(BT03t(frameRange), BT03_total_frames(frameRange),'k-','MarkerSize',1);
title(sprintf('%s 2x2bin miss %d/%d(%0.2f%%)',...
                        BT03IOI.IOI.subj_name,...
                        numel(BT03.missing_frames), BT03.n_frames*BT03.nColors, ...
                        100*numel(BT03.missing_frames)/(BT03.n_frames*BT03.nColors)),...
                        'interpreter','none','FontSize',14);
xlabel('t [s]','FontSize',14)
set(gca,'yTick',0)
