%% Breakpoint on line 221 of ioi_fc_GLM_on_ROI_run
% fprintf('\nGLM for %s ROI %d Session %d Color %d (%s) done!\n',IOI.subj_name,r1,s1,c1,colorNames{1+c1})
% Identify in IOI the file name of the time series
IOI.fcIOS.SPM(1).fnameROIregress = fnameROIregress;
% Get mask for each ROI
[~, mask] = ioi_get_ROImask(IOI,job);
Amask = []; % Initialize activation mask
% We are not extracting brain mask here
job.extractingBrainMask = false;
job.extractBrainMask = false;
% Extract ROI from regressed whole image series.
[ROIregressTmp IOI] = ...
    ioi_extract_core(IOI,job,mask,Amask,'regressData');

%% Plot
h = figure; set(h,'color','w')
subplot(311); plot(ROIregress{r1}{s1, c1}, 'k-', 'LineWidth', 1.5)
title('Regression on ROIs','FontSize',10)
ylabel('\DeltaHbO_2 (a.u.)','FontSize',10)
set(gca,'FontSize',8)
xlim([0 863]); ylim([-2.5 3.1]);
subplot(312); plot(ROIregressTmp{r1}{s1, c1}, 'k-', 'LineWidth', 1.5)
title('Regression on pixels','FontSize',10)
ylabel('\DeltaHbO_2 (a.u.)','FontSize',10)
xlim([0 863]); ylim([-2.5 3.1]);
set(gca,'FontSize',8)
subplot(313); plot(100*(ROIregress{r1}{s1, c1} - ROIregressTmp{r1}{s1, c1}) ./ ROIregress{r1}{s1, c1},...
    'k-', 'LineWidth', 1.5)
axis tight; xlim([0 863]); 
title('Error (%)','FontSize',10)
ylabel('Error (%)','FontSize',10)
xlabel('t(s)','FontSize',10)
set(gca,'FontSize',8)


%% Print
% Specify window units
set(h, 'units', 'inches')
% Change figure and paper size
set(h, 'Position', [0.1 0.1 6 3])
set(h, 'PaperPosition', [0.1 0.1 6 3])
% Save as PNG at the user-defined resolution
print(h, '-dpng', ...
    fullfile('D:\Edgar\Documents\Dropbox\Docs\fcOIS', 'GLM_test'),...
    sprintf('-r%d',300));
% Save as a figure
saveas(h, fullfile('D:\Edgar\Documents\Dropbox\Docs\fcOIS', 'GLM_test'), 'fig');
