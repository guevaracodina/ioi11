function out = ioi_GLM_on_ROI_run(job)
%spm_shift = 32; %annoying shift in spm, present on X and U, etc.
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
save_figures = job.save_figures;
generate_figures = job.generate_figures;
figure_show_stim = 0; %job.figure_show_stim;
figure_rebase_to_zero_at_stim = job.figure_rebase_to_zero_at_stim;
use_onset_amplitudes = job.use_onset_amplitudes;
include_flow = job.include_flow;
show_mse = job.show_mse;
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
                    else
                        newDir = dir_ioimat;
                    end
                    if save_figures
                        dir_fig = fullfile(newDir,'fig');
                        if ~exist(dir_fig,'dir'),mkdir(dir_fig);end
                    end
                    %loop over sessions
                    for s1=1:length(IOI.sess_res)
                        if all_sessions || sum(s1==selected_sessions)
                            %Electrophysiology, for each subject and session
                            
                            %TO-DO: generalize to more than 1 onset
                            %type
                            ot = 1;
                            ons = IOI.sess_res{s1}.onsets{ot}; %already in seconds *IOI.dev.TR;
                            dur = IOI.sess_res{s1}.durations{ot}; %*IOI.dev.TR;
                            name = IOI.sess_res{s1}.names{ot};
                            if use_onset_amplitudes
                                amp = IOI.sess_res{s1}.parameters{ot};
                            else
                                amp = [];
                            end
                            %convolve with hemodynamic response function
                            [Xtmp U] = ioi_get_X(IOI,name,ons,dur,amp,s1,bases,volt);
                            IOI.Sess(s1).U = U; %store onsets for each session
                            
                            %loop over available colors
                            for c1=1:length(IOI.sess_res{s1}.fname)
                                if ~iscell(Xtmp)
                                    X = Xtmp;
                                else
                                    X = Xtmp{c1};
                                end
                                if ~isempty(X)
                                    IOI.X{s1}.X0 = X;
                                    %filter X - HPF
                                    if hpf_butter_On
                                        X = ButterHPF(1/IOI.dev.TR,hpf_butter_freq,hpf_butter_order,X);
                                    end
                                    %add a constant
                                    X = [X ones(size(X,1),1)];
                                    %get K for low pass filtering:
                                    K = get_K(1:size(X,1),fwhm,IOI.dev.TR);
                                    %filter X - LPF
                                    %calculate X inverse
                                    %Xm = pinv(X);
                                    %Xu = X(:,1);
                                    %filter forward
                                    X1 = spm_filter_HPF_LPF_WMDL(K,X);
                                    %filter backward
                                    X1 = X1(end:-1:1,:);
                                    X1 = spm_filter_HPF_LPF_WMDL(K,X1);
                                    X1 = X1(end:-1:1,:);
                                    KX = spm_sp('Set', X1);
                                    KX.X = full(KX.X); %Filtered X
                                    Xm = spm_sp('x-',KX); % projector
                                    %covariance
                                    bcov = Xm * K.KL;
                                    bcov = bcov * bcov';
                                    %approximate calculation of effective degrees of freedom
                                    [trRV trRVRV] = approx_trRV(KX.X,Xm,K.KL);
                                    if ~iscell(Xtmp)
                                        cX = c1;
                                    else
                                        cX = 1;
                                    end
                                    IOI.X{s1}.X{cX} = X;
                                    IOI.X{s1}.Xm{cX} = Xm;
                                    IOI.X{s1}.bcov{cX} = bcov;
                                    IOI.X{s1}.trRV{cX} = trRV;
                                    IOI.X{s1}.trRVRV{cX} = trRVRV;
                                    IOI.X{s1}.erdf{cX} = (trRV)^2/trRVRV;
                                    %load ROI
                                    if ~isempty(job.ROImat)
                                        load(job.ROImat{SubjIdx});
                                    else
                                        try
                                            load(IOI.ROI.ROIfname);
                                        catch
                                            load(fullfile(dir_ioimat,'ROI.mat'));
                                        end
                                    end
                                    
                                    %Loop over ROIs
                                    for r1=1:length(ROI)
                                        if all_ROIs || sum(r1==selected_ROIs)
                                            y = ROI{r1}{s1,c1};
                                            if ~isempty(y)
                                                %yu = y; %unfiltered data
                                                %filtering of the data: HPF
                                                if hpf_butter_On
                                                    y =  ButterHPF(1/IOI.dev.TR,hpf_butter_freq,hpf_butter_order,y);
                                                end
                                                %yu = y;
                                                %filtering of the data: LPF (Gaussian), forward
                                                y = spm_filter_HPF_LPF_WMDL(K,y')';
                                                %filter backward
                                                y = y(end:-1:1);
                                                y = spm_filter_HPF_LPF_WMDL(K,y')';
                                                y = y(end:-1:1);
                                                %GLM inversion: calculating beta and residuals - would be SPM.Vbeta,
                                                b = Xm * y'; % beta : least square estimate
                                                %Compute t stat
                                                res = y'-X*b;
                                                res2 = sum(res.^2);
                                                IOI.X{s1}.b{r1,c1} = b;
                                                IOI.X{s1}.r(r1,c1) = res2; % Residuals
                                                IOI.X{s1}.mse(r1,c1) = res2/length(y);
                                                IOI.X{s1}.t(r1,c1) = b(1)/(res2*bcov(1,1)/trRV)^0.5;
                                                if volt
                                                    IOI.X{s1}.t2(r1,c1) = b(2)/(res2*bcov(2,2)/trRV)^0.5;
                                                end
                                                IOI.X{s1}.yf{r1,c1} = y; %filtered data
                                                %IOI.X{s1}.yu{r1,c1} = yu;
                                                IOI.X{s1}.yp{r1,c1} = X * b; %IOI.X{s1}.b{r1,c1}; %predicted
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                    IOI.res.GLMOK = 1;
                    save(IOImat,'IOI');
                    if ~exist('ROI','var')
                        %load ROI
                        if ~isempty(job.ROImat)
                            load(job.ROImat{SubjIdx});
                        else
                            try
                                load(IOI.ROI.ROIfname);
                            catch
                                load(fullfile(dir_ioimat,'ROI.mat'));
                            end
                        end
                    end
                    %Plot figures of fits
                    %line specification - (yp or yf) x color
                    lp1{1} = ':'; lp1{2} = '-'; %lp1{3} = '--'; lp1{4} = '-.';
                    ctotal = [];
                    if isfield(IOI.color,'HbO')
                        lp2{IOI.color.eng==IOI.color.HbO} = 'r'; %HbO
                        lp2{IOI.color.eng==IOI.color.HbR} = 'b'; %HbR
                        ctotal = [ctotal find(IOI.color.eng==IOI.color.HbO) ...
                            find(IOI.color.eng==IOI.color.HbR)];
                        lp3{IOI.color.eng==IOI.color.HbO} = 'HbO'; %HbO
                        lp3{IOI.color.eng==IOI.color.HbR} = 'HbR'; %HbR
                    end
                    if include_flow
                        if isfield(IOI.color,'flow')
                            lp2{IOI.color.eng==IOI.color.flow} = 'k'; %Flow
                            ctotal = [ctotal find(IOI.color.eng==IOI.color.flow)];
                            lp3{IOI.color.eng==IOI.color.flow} = 'Flow'; %Flow
                        end
                    end
                    %loop over sessions
                    h1 = 0;
                    for s1=1:length(IOI.sess_res)
                        if all_sessions || sum(s1==selected_sessions)
                            %Loop over ROIs
                            for r1=1:length(ROI)
                                if all_ROIs || sum(r1==selected_ROIs)
                                    h1 = h1 + 1;
                                    h(h1) = figure;
                                    legstr = {};
                                    %loop over available colors
                                    for c1=1:length(IOI.sess_res{s1}.fname) %or ctotal
                                        try
                                            if include_flow || ~(IOI.color.eng(c1)==IOI.color.flow)
                                                if ~isempty(IOI.X{s1}.yf{r1,c1})
                                                    if length(IOI.X{s1}.X)>1
                                                        X = IOI.X{s1}.X{c1};
                                                    else
                                                        X = IOI.X{s1}.X{1};
                                                    end
                                                    %lp = linspace(0,(size(X,1)-spm_shift)*IOI.dev.TR,size(X,1)-spm_shift);
                                                    lp = linspace(0,size(X,1)*IOI.dev.TR,size(X,1));
                                                    %filtered data
                                                    yf = IOI.X{s1}.yf{r1,c1};
                                                    %reconstructed data
                                                    yp = IOI.X{s1}.yp{r1,c1};
                                                    if figure_rebase_to_zero_at_stim
                                                        u0 = full(IOI.Sess(s1).U.u(33:end)');
                                                        u0(u0==0) = NaN;
                                                        for k0=1:length(u0)
                                                            if u0(k0) > 0
                                                                yf(k0:end) = yf(k0:end)-yf(k0);
                                                                yp(k0:end) = yp(k0:end)-yp(k0);
                                                            end
                                                        end
                                                    end
                                                    
                                                    t = IOI.X{s1}.t(r1,c1);
                                                    if volt==2, t2 = IOI.X{s1}.t2(r1,c1); end
                                                    plot(lp,yf,[lp1{1} lp2{c1}]); hold on
                                                    plot(lp,yp,[lp1{2} lp2{c1}]); hold on
                                                    legstr = [legstr; [lp3{c1} 'f']];
                                                    if show_mse
                                                        mse_str = [',' sprintf('%2.3f',IOI.X{s1}.mse(r1,c1))];
                                                    else
                                                        mse_str = '';
                                                    end
                                                        
                                                    if volt == 1
                                                        legstr = [legstr; [lp3{c1} 'p(' sprintf('%2.1f',t) mse_str ')']];
                                                    else
                                                        legstr = [legstr; [lp3{c1} 'p(' sprintf('%2.1f',t) ',' sprintf('%2.1f',t2) mse_str ')']];
                                                    end
                                                end
                                            end
                                        end
                                    end
                                    %plot stims
                                    if figure_show_stim
                                        u0 = full(IOI.Sess(s1).U.u(33:end)');
                                        u0(u0==0) = NaN;
                                        stem(lp,u0,'k');
                                    end
                                    tit = ['Session ' int2str(s1) ', ROI ' int2str(r1)];
                                    title(tit);
                                    legend(gca,legstr);
                                    ioi_save_figures(save_figures,generate_figures,h(h1),tit,dir_fig);
                                end
                            end
                        end
                    end
                end
            end
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
