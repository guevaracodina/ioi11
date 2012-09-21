function ioi_plot_roi_time_series(IOI,ROI,generate_figures,save_figures)
lsp{1} = '.c'; lsp{2} = '.m'; lsp{3} = '.g';
for cC = 6:7;
    for cROI = 1:2
        aROI = ROI{cROI}{cC};
        Ns = length(aROI);
        lp = linspace(0,Ns*TR,Ns);
        %add a HPF
        d = ButterHPF(1/TR,OP.bf,OP.bo,aROI);
        
        h = figure; plot(lp,d); hold on; plot(lp,tmp_d,'y'); plot(lp,aROI,'g')
        if ~isempty(rmi{i0})
            try
                stem(rmi{i0},10*ones(1,length(rmi{i0})),'xk');
                stem(armonsets{i0},20*ones(1,length(armonsets{i0})),'+r');
            end
        end
        for j0=1:3
            stem(onsets{i0,j0},5*ones(1,length(onsets{i0,j0})),lsp{j0});
        end
        tit = [IOI.subj_name ', Color ' IOI.color.eng(cC) ', ROI ' int2str(cROI)];
        title(tit);
        ioi_save_figures(save_figures,generate_figures,h,tit,pathGM);
    end
end