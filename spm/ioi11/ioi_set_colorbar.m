function hc1 = ioi_set_colorbar(gcf,len)
figure(gcf)
hc1 = colorbar;
set(hc1, 'YLim', [1 len+1]);
y_tick = linspace(1, len, len)'+0.49;
set(hc1, 'YTick', y_tick);
%set(hc1, 'YTickMode', 'Manual');
set(hc1, 'FontSize', 12);
%Customize here number of decimals
set(hc1,'YTickLabel',sprintf('%.0f |',get(hc1,'YTick')'));