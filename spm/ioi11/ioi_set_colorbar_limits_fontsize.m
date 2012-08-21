function hc = ioi_set_colorbar_limits_fontsize(hc,hc_min,hc_max,tick_number,fontsize_choice)
if hc_min==hc_max %quick fix
    hc_min=0.999*hc1_max; %-0.0001;
end
set(hc, 'YLim', [hc_min hc_max]);
y_tick = linspace(hc_min, hc_max, tick_number)';
set(hc, 'YTick', y_tick);
set(hc, 'FontSize', fontsize_choice);
%Customize here number of decimals
set(hc,'YTickLabel',sprintf('%.1f |',get(hc,'YTick')'));