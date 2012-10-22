%% Script to display temperature control in MouseSTAT
load('D:\Edgar\Data\IOS_Results\12_10_19,TT00\IOI.mat')
load('D:\Edgar\Data\IOS_Results\12_10_19,TT00\elinfo_S01.mat')

%% Retrieve data
tPad = ConvertedData.Data.MeasuredData(1,10).Data;
dt = ConvertedData.Data.MeasuredData(1,10).Property(1,3).Value*2;
tRat = ConvertedData.Data.MeasuredData(1,9).Data;
tVector = (0:ConvertedData.Data.MeasuredData(1,9).Total_Samples-1)'*dt;

%% Plot temperature
close all
figure; set(gcf, 'color', 'w')
h = plot(tVector, tPad, 'kx', ...
    tVector, tRat, 'r.',...
    tVector, 37*ones(size(tVector)), 'b:');
% IOI/Temp Pad ChanNames{1}{10}, IOI/Temp rat ChanNames{1}{11}
legend({ChanNames{1}{10} ChanNames{1}{11}})
set(h(3), 'LineWidth', 3)
set(gca,'FontSize', 12)
ylabel('Temp [\circC]', 'FontSize', 14)
xlabel('time [s]', 'FontSize', 14)
title(IOI.subj_name, 'FontSize', 14, 'FontWeight', 'Bold', 'Interpreter', 'none')


%% Print results
addpath(genpath('D:\Edgar\ssoct\Matlab'))
export_fig(fullfile('D:\Edgar\Documents\Dropbox\Docs\fcOIS\2012_10_22_Report',...
    'tempMouseSTAT'),'-png',gcf)
