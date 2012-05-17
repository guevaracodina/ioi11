function ioi_plot_LFP(IOI,el,s1)
%el = load(IOI.res.el{s1});
n0 = length(el);
lp0 = linspace(0,n0/(IOI.res.sfel),n0/10);
figure; plot(lp0,el(1:10:end)/IOI.res.sfel,'k');
xlabel('time (s)')
ylabel('LFP (mV)')
title(['LFP: Session ' int2str(s1)'])
