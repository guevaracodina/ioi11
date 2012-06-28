function ioi_display_electro(IOI)
%Plot electrophysiology for all sessions
try  
    ds = 10; %downsampling factor
    if isfield(IOI.res,'sfel')
        sfel = IOI.res.sfel;
    else
        sfel = 10000;
    end       
    for i0=1:length(IOI.res.el)
        el = load(IOI.res.el{i0});
        el = el.el;
        figure;
        n = length(el);
        lp = linspace(0,n/sfel,n/ds);
        plot(lp,el(1:ds:end),'k')
        tit = ['Session: ' int2str(i0)];
        title(tit);
        xlabel('Time (s)')
        ylabel('LFP (mV)')
    end
catch exception
    disp(exception.identifier)
    disp(exception.stack(1))
end