function ioi_save_figures(save_figures,generate_figures,h,tit,dir_fig)
%figure(h);
%title(tit);
if save_figures
    filen = fullfile(dir_fig,[tit '.tiff']); %save as .tiff
    print(h, '-dtiffn', filen);
    filen2 = fullfile(dir_fig,[tit '.fig']); %save as .fig
    saveas(h,filen2,'fig');
end
if ~generate_figures, close(h); end
end