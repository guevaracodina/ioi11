function ioi_plot_ROI_time_series(filtNdownROI, filtNdownBrain, rN, sN, cN)
% Display plots on SPM graphics window
% spm_figure('GetWin', 'Graphics');
% spm_figure('Clear', 'Graphics');
nCols = numel(sN);
nRows = numel(cN);
for r1 = rN,
    figure;
    nPlot = 1;
    for s1 = sN,
        for c1 = cN,
            subplot(nCols,nRows,nPlot); plot(filtNdownROI{r1}{s1,c1});
            title(sprintf('Filtered ROI %d Session %d Color %d',r1,s1,c1));
%             subplot(nCols,nRows,nPlot); plot(filtNdownBrain{r1}{s1,c1});
%             title(sprintf('Global brain signal,  Session %d Color %d',s1,c1));
            nPlot = nPlot + 1;
        end
    end
end
