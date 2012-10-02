%% Limitations of IIR filters
fs = 5; n = 5; Wn = 2*[0.009 0.08]/fs;
ftype = 'bandpass';

% Transfer Function design
[b,a] = butter(n,Wn,ftype);
h1=dfilt.df2(b,a);      % This is an unstable filter.

% Zero-Pole-Gain design
[z, p, k] = butter(n,Wn,ftype);
[sos,g]=zp2sos(z,p,k);
h2=dfilt.df2sos(sos,g);

% Plot and compare the results
hfvt=fvtool(h1,h2,'FrequencyScale','log');
legend(hfvt,'TF Design','ZPK Design')
hAx = findobj(get(hfvt,'Children'),'Type','axes');
set(hAx,'FontSize',12)
set(gcf,'Color','w')
set(get(hAx,'YLabel'),'FontSize',12)
set(get(hAx,'XLabel'),'FontSize',12)
set(get(hAx,'Title'),'FontSize',12)
set(get(hAx,'Title'),'String',sprintf('n = %d, f_s = %d Hz, W_n = [%0.3f %0.3f] Hz',n,fs,fs*Wn/2))
print(gcf, '-dpng', fullfile('D:\Edgar\Documents\Dropbox\Docs\fcOIS\2012_10_01_Report','filter_comp'), '-r300');
