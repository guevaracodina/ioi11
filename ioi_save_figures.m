function ioi_save_figures(save_figures,generate_figures,h,tit,dir_fig)
try
    fname = fullfile(dir_fig,tit);
    if save_figures
        print(h, '-dpng', [fname '.png'], '-r300');
        saveas(h,[fname '.fig']);
    end
    if ~generate_figures, close(h); end
end
