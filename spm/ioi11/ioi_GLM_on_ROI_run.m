function out = ioi_GLM_on_ROI_run(job)
%Volterra
volt = job.volt;
%bases
bases = job.bases;
%HPF filter
if isfield(job.hpf_butter,'hpf_butter_On')
    hpf_butter_On = 1;
    hpf_butter_freq = job.hpf_butter.hpf_butter_On.hpf_butter_freq;
    hpf_butter_order = job.hpf_butter.hpf_butter_On.hpf_butter_order;
else
    hpf_butter_On = 0;
end
%LPF filter
fwhm = job.lpf_gauss.fwhm1;

%select a subset of sessions
if isfield(job.session_choice,'select_sessions')
    all_sessions = 0;
    selected_sessions = job.session_choice.select_sessions.selected_sessions;
else
    all_sessions = 1;
end
%select a subset of ROIs
if isfield(job.ROI_choice,'select_ROIs')
    all_ROIs = 0;
    selected_ROIs = job.ROI_choice.select_ROIs.selected_ROIs;
else
    all_ROIs = 1;
end
%select onsets
if isfield(job.stim_choice,'electro_stims')
    electro_stims = 1;
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
else
    %default stims
    electro_stims = 0;
end

%Big loop over subjects
for SubjIdx=1:length(job.IOImat)
    try
        tic
        clear IOI ROI SPM
        %Load IOI.mat information
        IOImat = job.IOImat{SubjIdx};
        load(IOImat);
        if ~isfield(IOI.res,'ROIOK')
            disp(['No ROI available for subject ' int2str(SubjIdx) ' ... skipping series extraction']);
        else
            if ~isfield(IOI.res,'seriesOK')
                disp(['Extracted series not available for subject ' int2str(SubjIdx) ' ... skipping GLM']);
            else
                if ~isfield(IOI.res,'GLMOK') || job.force_redo
                    if isfield(IOI,'X'), IOI = rmfield(IOI,'X'); end
                    [dir_ioimat dummy] = fileparts(job.IOImat{SubjIdx});
                    if isfield(job.IOImatCopyChoice,'IOImatCopy')
                        newDir = job.IOImatCopyChoice.IOImatCopy.NewIOIdir;
                        newDir = fullfile(dir_ioimat,newDir);
                        if ~exist(newDir,'dir'),mkdir(newDir); end
                        IOImat = fullfile(newDir,'IOI.mat');
                        cdir = newDir; %current directory
                    else
                        cdir = dir_ioimat;
                    end
                    %loop over sessions
                    for s1=1:length(IOI.sess_res)
                        if all_sessions || sum(s1==selected_sessions)
                            %Electrophysiology, for each subject and session
                            if electro_stims
                                [pkh ons] = ioi_get_onsets(IOI,s1,E,cdir); %pk in seconds; pkh in arbitrary units
                                dur = 1;
                                name = '';
                                %                                 separate_slow_fast = 1;
                                %                                 if separate_slow_fast
                                %                                     name{1}
                                %                                 end
                            else
                                ons = IOI.sess_res{s1}.onsets{1}; %already in seconds *IOI.dev.TR;
                                dur = IOI.sess_res{s1}.durations{1}; %*IOI.dev.TR;
                                name = IOI.sess_res{s1}.names{1};
                            end
                            %convolve with hemodynamic response function
                            [X U] = ioi_get_X(IOI,name,ons,dur,s1,bases,volt);
                            IOI.Sess(s1).U = U; %store onsets for each session
                            IOI.X{s1}.X0 = X;
                            %filter X - HPF
                            if hpf_butter_On
                                X = ButterHPF(1/IOI.dev.TR,hpf_butter_freq,hpf_butter_order,X);
                            end
                            %add a constant
                            X = [X ones(size(X,1),1)];
                            %get K for low pass filtering:
                            K = get_K(X,fwhm,IOI.dev.TR);
                            %filter X - LPF
                            %calculate X inverse
                            %Xm = pinv(X);
                            KX = spm_sp('Set', spm_filter_HPF_LPF_WMDL(K,X));
                            KX.X = full(KX.X); %Filtered X
                            Xm = spm_sp('x-',KX); % projector
                            %covariance
                            bcov = Xm * K.KL;
                            bcov = bcov * bcov';
                            IOI.X{s1}.X = X;
                            IOI.X{s1}.Xm = Xm;
                            IOI.X{s1}.bcov = bcov;
                            %load ROI
                            load(IOI.ROI.ROIfname);
                            %approximate calculation of effective degrees of freedom
                            [trRV trRVRV] = approx_trRV(KX.X,Xm,K.KL);
                            IOI.X{s1}.trRV = trRV;
                            IOI.X{s1}.trRVRV = trRVRV;
                            IOI.X{s1}.erdf = (trRV)^2/trRVRV;
                            %loop over available colors
                            for c1=1:length(IOI.sess_res{s1}.fname)
                                if ~(IOI.color.eng(c1)==IOI.color.laser)
                                    %Loop over ROIs
                                    for r1=1:length(ROI)
                                        if all_ROIs || sum(r1==selected_ROIs)
                                            y = ROI{r1}{s1,c1};
                                            if ~isempty(y)
                                                %filtering of the data: HPF
                                                if hpf_butter_On
                                                    y =  ButterHPF(1/IOI.dev.TR,hpf_butter_freq,hpf_butter_order,y);
                                                end
                                                %filtering of the data: LPF (Gaussian)
                                                y = spm_filter_HPF_LPF_WMDL(K,y')';
                                                %GLM inversion: calculating beta and residuals - would be SPM.Vbeta,
                                                b = Xm * y'; % beta : least square estimate
                                                IOI.X{1,s1}.b{r1,c1} = b;
                                                %Compute t stat
                                                res = y'-X*b;
                                                res2 = sum(res.^2);
                                                IOI.X{1,s1}.r(r1,c1) = res2; % Residuals
                                                IOI.X{1,s1}.t(r1,c1) = b(1)/(res2*bcov(1,1)/trRV)^0.5;
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                    
                    IOI.res.GLMOK = 1;
                    save(IOImat,'IOI');
                end
            end
        end
        toc
        disp(['Subject ' int2str(SubjIdx) ' complete']);
        out.IOImat{SubjIdx} = IOImat;
        
    catch exception
        disp(exception.identifier)
        disp(exception.stack(1))
    end
end
end
function [pkh pk] = ioi_get_onsets(IOI,s1,E,cdir)
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
    load(IOI.res.el{s1});
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
        el = ButterLPF(sf,cutoff,order,el0);
        if E.write_pictures
            h = figure; plot(lp,sgn*el0(1:ds:end),'k'); hold on; plot(lp,-el(1:ds:end),'r'); hold off
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
    if E.write_pictures
        h = figure; plot(lp,sgn*el(1:ds:end),'k'); hold on; stem(pk/sf,sgn*(MN+nSD*eSD)*ones(1,length(pk)),'r'); 
        stem(25,sgn*MN,'b'); stem(50,sgn*SD0,'g'); %stem(75,sgn*(MN+nSD*SD),'m'); hold off
        filen1 = fullfile(fdir,['el_stem_tmp_S' gen_num_str(s1,2) '.fig']);
        saveas(h,filen1,'fig'); %save as .fig
        filen2 = fullfile(fdir,['el_stem_tmp_S' gen_num_str(s1,2) '.tiff']);
        print(h, '-dtiffn', filen2);
        close(h);
    end
    npk = [];
    npkh = [];
    cf = 0.5;
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
        filen1 = fullfile(fdir,['el_hist_zoom_S' gen_num_str(s1,2) '.fig']);
        saveas(h,filen1,'fig'); %save as .fig
        filen2 = fullfile(fdir,['el_hist_zoom_S' gen_num_str(s1,2) '.tiff']);
        print(h, '-dtiffn', filen2);
        close(h);       
    end
    A = private_analyze_LFP(E,el,pk,IOI,s1);
catch exception
    disp(exception.identifier)
    disp(exception.stack(1))
end
end

function K = get_K(X,fwhm,TR)
svec = 1:size(X,1);
K.HParam.type = 'none';
K.LParam.FWHM = fwhm;
K.LParam.type = 'Gaussian';
K = struct( 'HParam', K.HParam,...
    'row', svec ,...
    'RT', TR,...
    'LParam', K.LParam);
K = spm_filter_HPF_LPF_WMDL(K);
K.row = 1:length(svec);
end

function A = private_analyze_LFP(E,el,ons,IOI,s)
    n = length(ons);
    sf = E.sf; %Hz, sampling frequency
    tb = 0.2; %time before, in seconds
    ta = 0.5; %time after, in seconds
    Pf = 10; %rescaling factor for power, to put it on a similar scale in plots
    lb = ceil(sf*tb);
    la = ceil(sf*ta);
    m = zeros(n,lb+la);
    lp = linspace(-tb,ta,(ta+tb)*sf);
    for i=1:n
        st = ceil(ons(i)*sf);
        try
            %store electrophysiology data
            m(i,:) = el(st-lb+1:st+la);
        catch
            %data after last onset or before first onset might be incomplete
        end
    end
    figure; plot(lp,m')
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