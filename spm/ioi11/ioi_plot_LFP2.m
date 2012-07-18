function ioi_plot_LFP2(IOI,el,s1,el2,DS,save,fname)
%el = load(IOI.res.el{s1});
n0 = length(el);
lp0 = linspace(0,n0/(IOI.res.sfel),n0/DS);
h = figure; plot(lp0,el(1:DS:end)/IOI.res.sfel,'k'); hold on
plot(lp0,-el2(1:DS:end)/IOI.res.sfel,'r'); hold on
xlabel('time (s)')
ylabel('LFP (mV)')
title(['LFP: Session ' int2str(s1)'])
if save
    print(h, '-dpng', [fname '.png'], '-r300');
    saveas(h,[fname '.fig']);
    close(h);
end
