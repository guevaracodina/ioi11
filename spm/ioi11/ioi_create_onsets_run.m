function out = ioi_create_onsets_run(job)
%select onsets, default is stimulation based
stim_choice=0;
if isfield(job.stim_choice,'electro_stims')
    stim_choice = 1;
    %sampling frequency
    E.sf = job.stim_choice.electro_stims.sf;
    %number of standard deviations
    E.nSD = job.stim_choice.electro_stims.nSD;
    %minimum standard deviation imposed
    E.mbSD = job.stim_choice.electro_stims.mbSD;
    %minimal distance between peaks in seconds
    E.dP = job.stim_choice.electro_stims.dP/1000;
    %HPF
    if isfield(job.stim_choice.electro_stims.electro_hpf_butter,'electro_hpf_butter_On')
        E.electro_hpf_butter_On = 1;
        E.hpf_butter_freq = job.stim_choice.electro_stims.electro_hpf_butter.electro_hpf_butter_On.electro_hpf_butter_freq;
        E.hpf_butter_order = job.stim_choice.electro_stims.electro_hpf_butter.electro_hpf_butter_On.electro_hpf_butter_order;
    else
        E.electro_hpf_butter_On = 0;
    end
    %LPF
    if isfield(job.stim_choice.electro_stims.electro_lpf_butter,'electro_lpf_butter_On')
        E.electro_lpf_butter_On = 1;
        E.lpf_butter_freq = job.stim_choice.electro_stims.electro_lpf_butter.electro_lpf_butter_On.electro_lpf_butter_freq;
        E.lpf_butter_order = job.stim_choice.electro_stims.electro_lpf_butter.electro_lpf_butter_On.electro_lpf_butter_order;
    else
        E.electro_lpf_butter_On = 0;
    end
    E.write_pictures = job.stim_choice.electro_stims.write_pictures;
    E.use_epilepsy_convention = job.stim_choice.electro_stims.use_epilepsy_convention;
    E.electrophysiology_onset_name = job.stim_choice.electro_stims.electrophysiology_onset_name;
    E.tb = job.stim_choice.electro_stims.tb;
    E.ta = job.stim_choice.electro_stims.ta;
end
if isfield(job.stim_choice,'manual_stims')
    stim_choice = 2;
    if ~isempty(job.stim_choice.manual_stims.onset_type_info)
        oti = job.stim_choice.manual_stims.onset_type_info;
    else
        oti = [];
    end
end

%select a subset of sessions
if isfield(job.session_choice,'select_sessions')
    all_sessions = 0;
    selected_sessions = job.session_choice.select_sessions.selected_sessions;
else
    all_sessions = 1;
end

elDir = job.elDir;
%Big loop over subjects
for SubjIdx=1:length(job.IOImat)
    try
        tic
        clear IOI
        %Load IOI.mat information
        IOImat = job.IOImat{SubjIdx};               
        [dir_ioimat dummy] = fileparts(job.IOImat{SubjIdx});
        if isfield(job.IOImatCopyChoice,'IOImatCopy')
            newDir = job.IOImatCopyChoice.IOImatCopy.NewIOIdir;
            newDir = fullfile(dir_ioimat,newDir);
            if ~exist(newDir,'dir'),mkdir(newDir); end
            IOImat = fullfile(newDir,'IOI.mat');
        else
            newDir = dir_ioimat;
        end
        try
            load(IOImat);
        catch
            load(job.IOImat{SubjIdx});
        end
        if ~isfield(IOI.res,'OnsetsOK') || job.force_redo
            %loop over sessions
            for s1=1:length(IOI.sess_res)
                if all_sessions || sum(s1==selected_sessions)
                    if all_sessions
                        selected_sessions = 1:length(IOI.sess_res);
                    end
                    switch stim_choice
                        case 0
                            %Just check that default stimulations are present
                            if ~isfield(IOI.sess_res{s1},'onsets')
                                disp(['No onsets for subject ' int2str(SubjIdx) ', session ' int2str(s1)]);
                            end
                        case 1
                            %Electrophysiology, for each subject and session
                            %try to specify a valid folder for electrophysiology
                            elDir0 = [];
                            if ~isempty(elDir)
                                if length(job.IOImat) == length(elDir)
                                    elDir0 = elDir{SubjIdx};
                                end
                            end
                            [pkh ons] = ioi_get_onsets(IOI,s1,E,newDir,elDir0); %pk in seconds; pkh in arbitrary units
                            ot = 1;
                            IOI.sess_res{s1}.E = E; %Electrophysiology structure used for detection
                            IOI.sess_res{s1}.names{ot} = E.electrophysiology_onset_name;
                            IOI.sess_res{s1}.onsets{ot} = ons;
                            IOI.sess_res{s1}.durations{ot} = IOI.dev.TR;
                            IOI.sess_res{s1}.parameters{ot} = pkh;
                        case 2
                            if isempty(oti)
                                spm_input(['Subject ' int2str(SubjIdx) ', Session ' int2str(s1)],'-1','d');
                                p = 1;
                                ot = 0;
                                while p
                                    p = spm_input('Add an onset type?',2,'y/n');
                                    if p == 'y', p = 1; else p = 0; end
                                    if p
                                        ot = ot + 1;
                                        IOI.sess_res{s1}.names{ot} = spm_input('Enter onset name',3,'s');
                                        IOI.sess_res{s1}.onsets{ot} = spm_input('Enter onsets in seconds',4,'e','',NaN);
                                        IOI.sess_res{s1}.durations{ot} = spm_input('Enter durations in seconds',5,'e','',NaN);
                                        IOI.sess_res{s1}.parameters{ot} = spm_input('Enter amplitude in seconds',6,'e','',NaN);
                                        try close(h2); end
                                    end
                                end
                            else
                                %find whether specified for each subject or to apply to all subjects
                                if length(oti) == length(job.IOImat) || length(oti) == 1
                                    %OK
                                    if length(oti) == 1
                                        subj = 1;
                                    else
                                        subj = SubjIdx;
                                    end
                                    if length(oti{subj}) == length(selected_sessions) || length(oti{subj}) == 1
                                        if length(oti{subj}) == 1
                                            sess = 1;
                                        else
                                            sess = find(s1==selected_sessions);
                                        end
                                        %loop over onset types
                                        for ot=1:length(oti{subj}{sess})
                                            IOI.sess_res{s1}.names{ot} = oti{subj}{sess}(ot).onset_name;
                                            IOI.sess_res{s1}.onsets{ot} = oti{subj}{sess}(ot).onset_times;
                                            IOI.sess_res{s1}.durations{ot} = oti{subj}{sess}(ot).onset_durations;
                                            IOI.sess_res{s1}.parameters{ot} = oti{subj}{sess}(ot).onset_amplitude;
                                        end
                                    else
                                        disp('Problem with number of sessions for which onsets were specified');
                                    end
                                else
                                    disp('Problem with number of subjects for which onsets were specified');
                                end
                            end
                    end
                end
            end
            IOI.res.OnsetsOK = 1;
            save(IOImat,'IOI');
        end
        
        toc
        disp(['Subject ' int2str(SubjIdx) ' complete']);
        out.IOImat{SubjIdx} = IOImat;
    catch exception
        disp(exception.identifier)
        disp(exception.stack(1))
        out.IOImat{SubjIdx} = job.IOImat{SubjIdx};
    end
end
end
function [pkh pk] = ioi_get_onsets(IOI,s1,E,cdir,elDir)
try
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
        load(IOI.res.el{s1});
    catch
        [dir0 fil0 ext0] = fileparts(IOI.res.el{s1});
        fil = fullfile(elDir,[fil0 ext0]);
        load(fil);
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
            h = figure; plot(lp,sgn*el0(1:ds:end),'k'); hold on; plot(lp,-el(1:ds:end),'r');
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
            h = figure; plot(lp,sgn*el0(1:ds:end),'k'); hold on; plot(lp,-el(1:ds:end),'r'); hold off
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
        disp(['No onsets detected initially. Chosen minimum SD too high, setting it to zero']);
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
    npk = [];
    npkh = [];
    cf = 0.5; %?
    %Narrow down to good events, by making sure it goes above then below baseline
    for i=1:length(pk)
        good = 0;
        %check that next peak is at least dP points later
        if (i>1 && pk(i)-pk(i-1) > dP) || (i==1)
            if i== 1
                pkm = -Inf;
            else
                %minimum over interval of 2*dP points
                pkm = min(el(pk(i-1):pk(i)));
            end
            if pkm < MN-cf*eSD
                %went below baseline, found new peak
                
                %calculate max
                [dummy pkM]= max(el(pk(i):min(length(el),(pk(i)+3*fast_pk))));
                
                if select_fast
                    if pkM <= fast_pk
                        %if pkM > 1
                        good =1;
                        %end
                    end
                end
                if select_slow
                    if  pkM> fast_pk
                        good = 1;
                    end
                end
                %exclude if previous spike was very close and current spike
                %is not deep -- to do
                
                if good
                    npk =[npk pk(i)];
                    npkh = [npkh pkh(i)];
                end
            end
        end
    end
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
    A = private_analyze_LFP(E,el,pk,IOI,s1,fdir);
catch exception
    disp(exception.identifier)
    disp(exception.stack(1))
end
end

function A = private_analyze_LFP(E,el,ons,IOI,s,fdir)
n = length(ons);
sf = E.sf; %Hz, sampling frequency
tb = E.tb; %0.2; %time before, in seconds
ta = E.ta; %0.5; %time after, in seconds
Pf = 10; %rescaling factor for power, to put it on a similar scale in plots
lb = ceil(sf*tb);
la = ceil(sf*ta);
m = zeros(n,lb+la);
lp = linspace(-tb,ta,(ta+tb)*sf);
if E.use_epilepsy_convention
    el = -el;
end
for i=1:n
    st = ceil(ons(i)*sf);
    try
        %store electrophysiology data
        m(i,:) = el(st-lb+1:st+la);
    catch
        %data after last onset or before first onset might be incomplete
    end
end
h = figure; plot(lp,m','k')
xlabel('seconds')
filen1 = fullfile(fdir,['Onset_profile_S' gen_num_str(s,2) '.fig']);
saveas(h,filen1,'fig'); %save as .fig
filen2 = fullfile(fdir,['Onset_profile_S' gen_num_str(s,2) '.tiff']);
print(h, '-dtiffn', filen2);
close(h);

h = figure; plot(lp,m') %with colors
xlabel('seconds')
filen1 = fullfile(fdir,['Onset_profile_colors_S' gen_num_str(s,2) '.fig']);
saveas(h,filen1,'fig'); %save as .fig
filen2 = fullfile(fdir,['Onset_profile_colors_S' gen_num_str(s,2) '.tiff']);
print(h, '-dtiffn', filen2);
close(h);

%Calculate some properties of the spike
%max amplitude
aM = max(m(:,lb:lb+la),[],2);
%min amplitude
am = min(m(:,lb:lb+la),[],2);
%difference between max and min amplitude
ad = aM-am;
%area
A0 = mean(abs(m(:,lb:lb+la)),2);
%signed area
sA0 = mean(m(:,lb:lb+la),2);
%power
P = Pf*mean(m(:,lb:lb+la).*m(:,lb:lb+la),2);
%baseline
b = mean(m(:,1:lb),2);
%previous quantities normalized to baseline:
mb = m-repmat(b,1,lb+la);
%max amplitude wrt baseline
aMb = max(mb(:,lb:lb+la),[],2);
%min amplitude wrt baseline
amb = min(mb(:,lb:lb+la),[],2);
%area wrt baseline
Ab = mean(abs(mb(:,lb:lb+la)),2);
%signed area wrt baseline
sAb = mean(mb(:,lb:lb+la),2);
%power wrt baseline
Pb = Pf*mean(mb(:,lb:lb+la).*mb(:,lb:lb+la),2);

%Store
A{s}.n = n;
A{s}.m = m;
A{s}.b = b;
A{s}.aM = aM;
A{s}.am = am;
A{s}.ad = ad;
A{s}.A = A0;
A{s}.P = P;
A{s}.mb = mb;
A{s}.aMb = aMb;
A{s}.amb = amb;
A{s}.Ab = Ab;
A{s}.sAb = sAb;
A{s}.Pb = Pb;
A{s}.Y = [aM am ad A0 sA0 P b aMb amb Ab sAb Pb];
%Protocole
try
    A{s}.Pt = IOI.sess_res{s}.scan_info.protocoleStimulation;
end
A{s}.ons = ons;
%clear aM am ad A sA P b mb aMb amb Ab sAb Pb
end