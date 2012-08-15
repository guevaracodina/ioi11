function [pkh pk dur] = ioi_get_onsets_from_electrophysiology(IOI,s1,E,cdir,elDir)
try
    dur = [];
    %sampling frequency
    sf = E.sf;
    nP = E.dP;
    nSD = E.nSD;
    mbSD = E.mbSD;
    %resolution in data points - 5 milliseconds
    rs = 50; %
    %select fast and slow rising onsets
    select_fast = 1; %peak reached within 40 ms of onsets
    fast_pk = 400; %data points at 10000 Hz
    select_slow = 1;
    %method:
    %1: detection of spikes based on thresholds
    %2: detection of spikes based on fitting gamma functions -- not coded up
    %method = 1;
    %load raw electrophysiology vector
    try
        if E.el_choice == 1
            load(IOI.res.el{s1});
        else
            %             [dir0 fil0 ext0] = fileparts(IOI.res.el{s1});
            load(IOI.res.el2{s1});
            el = el2;
        end
    catch
        [dir0 fil0 ext0] = fileparts(IOI.res.el{s1});
        fil = fullfile(elDir,[fil0 ext0]);
        load(fil);
    end
    if mean(abs(el))>10
        el = el/E.sf; %normalization
    end
    %remove time stamps for actual or spurious stimulations
    %ind = el>2; ind2 = [ind(7:end) false false false false false false];
    %el(ind)=el(ind2);
    el0 = el;
    if E.use_epilepsy_convention, sgn = -1; else sgn = 1; end
    if E.write_pictures
        fdir = fullfile(cdir,'fig');
        if ~exist(fdir,'dir'), mkdir(fdir); end
        ns = length(el);
        ds = 10;
        lp = linspace(0,ns/sf,ns/10);
    end
    %apply high pass filter to remove drifts
    if E.electro_hpf_butter_On
        cutoff = E.hpf_butter_freq; %in Hz
        order = E.hpf_butter_order;
        el = ButterHPF(sf,cutoff,order,el0);
        if E.write_pictures
            h = figure; plot(lp,sgn*el0(1:ds:end),'k'); hold on; plot(lp,sgn*el(1:ds:end),'r');
            xlabel('time (s)')
            ylabel('LFP (a.u.)')
            title('Detected onsets')
            filen1 = fullfile(fdir,['el_hpf_S' gen_num_str(s1,2) '.fig']);
            saveas(h,filen1,'fig'); %save as .fig
            filen2 = fullfile(fdir,['el_hpf_S' gen_num_str(s1,2) '.tiff']);
            print(h, '-dtiffn', filen2);
            close(h);
        end
    end
    %LPF
    if E.electro_lpf_butter_On
        cutoff = E.lpf_butter_freq; %in Hz
        order = E.lpf_butter_order;
        %forward
        el = ButterLPF(sf,cutoff,order,el);
        %backward
        el = el(end:-1:1);
        el = ButterLPF(sf,cutoff,order,el);
        el = el(end:-1:1);
        if E.write_pictures
            h = figure; plot(lp,sgn*el0(1:ds:end),'k'); hold on; plot(lp,sgn*el(1:ds:end),'r'); hold off
            xlabel('time (s)')
            ylabel('LFP (a.u.)')
            title('Detected onsets')
            filen1 = fullfile(fdir,['el_hpf_lpf_S' gen_num_str(s1,2) '.fig']);
            saveas(h,filen1,'fig'); %save as .fig
            filen2 = fullfile(fdir,['el_hpf_lpf_S' gen_num_str(s1,2) '.tiff']);
            print(h, '-dtiffn', filen2);
            close(h);
        end
    end
    
    SD0 = std(el);
    MN = mean(el);
    dP = floor(sf*nP);
    eSD = max(SD0,mbSD);
    %find onsets peaks: pkh: peak height; pk: onset time at sf sampling frequency
    [pkh pk] = findpeaks(el,'MINPEAKHEIGHT',MN+nSD*max(SD0,mbSD),'MINPEAKDISTANCE',rs);
    if isempty(pk)
        disp('No onsets detected initially. Chosen minimum SD too high, setting it to zero');
        [pkh pk] = findpeaks(el,'MINPEAKHEIGHT',MN+nSD*SD0,'MINPEAKDISTANCE',rs);
    end
    %[pkh2 pk2] = findpeaks(-el,'MINPEAKHEIGHT',MN+nSD*max(SD0,mbSD),'MINPEAKDISTANCE',rs);
    if E.write_pictures
        h = figure; plot(lp,sgn*el(1:ds:end),'k'); hold on; stem(pk/sf,sgn*(MN+nSD*eSD)*ones(1,length(pk)),'r');
        xlabel('time (s)')
        ylabel('LFP (a.u.)')
        title('Detected onsets')
        stem(25,sgn*MN,'b'); stem(50,sgn*SD0,'g'); %stem(75,sgn*(MN+nSD*SD),'m'); hold off
        filen1 = fullfile(fdir,['el_stem_tmp_S' gen_num_str(s1,2) '.fig']);
        saveas(h,filen1,'fig'); %save as .fig
        filen2 = fullfile(fdir,['el_stem_tmp_S' gen_num_str(s1,2) '.tiff']);
        print(h, '-dtiffn', filen2);
        close(h);
    end
    [npkh npk dur] = ioi_find_good_peaks(pk,pkh,dP,el,fast_pk,E,MN,eSD);
    pk = npk/sf;
    pkh = npkh;
    if E.write_pictures
        h = figure; plot(lp,sgn*el(1:ds:end),'k'); hold on; stem(pk,sgn*(MN+nSD*eSD)*ones(1,length(pk)),'r');
        xlabel('time (s)')
        ylabel('LFP (a.u.)')
        title('Detected onsets')
        stem(25,sgn*MN,'b'); stem(50,sgn*SD0,'g'); %stem(75,sgn*(MN+nSD*SD),'m'); hold off
        filen1 = fullfile(fdir,['el_stem_S' gen_num_str(s1,2) '.fig']);
        saveas(h,filen1,'fig'); %save as .fig
        filen2 = fullfile(fdir,['el_stem_S' gen_num_str(s1,2) '.tiff']);
        print(h, '-dtiffn', filen2);
        close(h);
    end
    %histogram of spike intervals
    if E.write_pictures
        df = diff(pk);
        h = figure; hist(df,120);
        xlabel('time (s)')
        ylabel('Number of onset intervals')
        title('Histogram of onset intervals')
        legend(sprintf( '%s\n%s\n%s', ['Total spikes: ' int2str(length(pk))],...
            ['Mean filt el: ' num2str(MN,'%.3f')],...
            ['SD filt el: ' num2str(SD0,'%.3f')]));
        filen1 = fullfile(fdir,['el_hist_S' gen_num_str(s1,2) '.fig']);
        saveas(h,filen1,'fig'); %save as .fig
        filen2 = fullfile(fdir,['el_hist_S' gen_num_str(s1,2) '.tiff']);
        print(h, '-dtiffn', filen2);
        close(h);
        df0 = df(df<1); %spike intervals less than 1 s
        h = figure; hist(df0,120); legend(['Spike intervals < 1 s: ' int2str(length(df0))]);
        xlabel('time (s)')
        ylabel('Number of onset intervals')
        title('Zoomed-in histogram of onset intervals')
        filen1 = fullfile(fdir,['el_hist_zoom_S' gen_num_str(s1,2) '.fig']);
        saveas(h,filen1,'fig'); %save as .fig
        filen2 = fullfile(fdir,['el_hist_zoom_S' gen_num_str(s1,2) '.tiff']);
        print(h, '-dtiffn', filen2);
        close(h);
    end
    A = ioi_analyze_LFP(E,el,pk,IOI,s1,fdir);
catch exception
    disp(exception.identifier)
    disp(exception.stack(1))
end