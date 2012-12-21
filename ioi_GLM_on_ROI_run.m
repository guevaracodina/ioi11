function out = ioi_GLM_on_ROI_run(job)
%spm_shift = 32; %annoying shift in spm, present on X and U, etc.
%Volterra
volt = job.volt;
%bases
bases = job.bases;

[all_sessions selected_sessions] = ioi_get_sessions(job);
[all_ROIs selected_ROIs] = ioi_get_ROIs(job);
%filters
HPF = ioi_get_HPF(job);
LPF = ioi_get_LPF(job);

save_figures = job.save_figures;
generate_figures = job.generate_figures;
figure_show_stim = 0; %job.figure_show_stim;
figure_rebase_to_zero_at_stim = job.figure_rebase_to_zero_at_stim;
use_onset_amplitudes = job.use_onset_amplitudes;
include_flow = job.IC.include_flow; %other colors not yet supported
show_mse = job.show_mse;
onset_choice = job.onset_choice;
%Big loop over subjects
for SubjIdx=1:length(job.IOImat)
    try
        tic
        clear ROI SPM
        %Load IOI.mat information
        [IOI IOImat dir_ioimat]= ioi_get_IOI(job,SubjIdx);
        
        if ~isfield(IOI.res,'ROIOK')
            disp(['No ROI available for subject ' int2str(SubjIdx) ' ... skipping series extraction']);
        else
            if ~isfield(IOI.res,'seriesOK')
                disp(['Extracted series not available for subject ' int2str(SubjIdx) ' ... skipping GLM']);
            else
                if ~isfield(IOI.res,'GLMOK') || job.force_redo
                    if isfield(IOI,'X'), IOI = rmfield(IOI,'X'); end
                    if save_figures
                        dir_fig = fullfile(dir_ioimat,'fig');
                        if ~exist(dir_fig,'dir'),mkdir(dir_fig);end
                    end
                    %loop over sessions
                    for s1=1:length(IOI.sess_res)
                        if all_sessions || sum(s1==selected_sessions)
                            %Electrophysiology, for each subject and session
                            
                            %TO-DO: generalize to more than 1 onset
                            %type
                            %******by Cong on 12/1/2/19                                                       
                            if isfield(job.do_F_test,'F_test_enabled')
                                model_choice = job.do_F_test.F_test_enabled.model_choice;
                                onset_choice = 0;
%                                 if model_choice                                
%                                 onset_choice = [0 1]; %************reduced model is the detection
%                                 else                            
%                                 onset_choice = [0 2];   %**************reduced model is the stim
%                                 end
                            end
                            
                             
                       
                            switch onset_choice
                                case 0 %**onsets from stim and spontaneous activity
                                    %onsets from detection
                                    ot = 1;
                                    ons{ot} = IOI.sess_res{s1}.onsets{ot}; %already in seconds *IOI.dev.TR;
                                    dur{ot} = IOI.sess_res{s1}.durations{ot}; %*IOI.dev.TR;
                                    name{ot} = IOI.sess_res{s1}.names{ot};
                                    if use_onset_amplitudes
                                        amp{ot} = IOI.sess_res{s1}.parameters{ot}';
                                    else
                                        amp{ot} = [];
                                    end
                                    
                                  
                                    %onsets from stimulation.
                                    ot = 2;
                                    ons{ot} = IOI.sess_res{s1}.onsets{ot}; %already in seconds *IOI.dev.TR;
                                    dur{ot} = IOI.sess_res{s1}.durations{ot}; %*IOI.dev.TR;
                                    name{ot} = IOI.sess_res{s1}.names{ot};
                                    if use_onset_amplitudes
                                        amp{ot} = IOI.sess_res{s1}.parameters{ot}';
                                    else
                                        amp{ot} = [];
                                    end
                                case 1 %***************onsets from detection
                                    ot = 1;
                                    ons = IOI.sess_res{s1}.onsets{ot}; %already in seconds *IOI.dev.TR;
                                    dur = IOI.sess_res{s1}.durations{ot}; %*IOI.dev.TR;
                                    name = IOI.sess_res{s1}.names{ot};
                                    if use_onset_amplitudes
                                        amp = IOI.sess_res{s1}.parameters{ot}';
                                    else
                                        amp = [];
                                    end                                    
                                case 2   %*************onsets from stim
                                    ot = 2;
                                    ons= IOI.sess_res{s1}.onsets{ot}; %already in seconds *IOI.dev.TR;
                                    dur = IOI.sess_res{s1}.durations{ot}; %*IOI.dev.TR;
                                    name = IOI.sess_res{s1}.names{ot};
                                    if use_onset_amplitudes
                                        amp = IOI.sess_res{s1}.parameters{ot}';
                                    else
                                        amp = [];
                                    end
                            end
                            
                            %***************end
                            %convolve with hemodynamic response function
                            [Xtmp U] = ioi_get_X(IOI,name,ons,dur,amp,s1,bases,volt);
                            IOI.Sess(s1).U = U; %store onsets for each session
                            
                            %loop over available colors
                            for c1=1:length(IOI.sess_res{s1}.fname)
                                if ~iscell(Xtmp)
                                    X = Xtmp;
                                else
%                                     X = Xtmp{c1};
                                   %**********by Cong
                                    X_spikes = Xtmp{1};
                                    X_stim = Xtmp{2};
                                    X=[X_spikes X_stim];
                                    %**end
                                end
                                if ~isempty(X)
                                    IOI.X{s1}.X0 = X;
                                    %filter X - HPF
                                    if HPF.hpf_butter_On
                                        X = ButterHPF(1/IOI.dev.TR,HPF.hpf_butter_freq,HPF.hpf_butter_order,X);
                                    end
                                    
                                    if isfield(job.do_F_test,'F_test_enabled')                                    
                                        X = [X ones(size(X,1),1)];
                                        if model_choice                                        
                                            X_spikes = [X(1) ones(size(X(1),1),1)];
                                        else
                                            X_stim = [X(2) ones(size(X(2),1),1)];
                                        end                                      
                                    else                                    
                                    %add a constant
                                    X = [X ones(size(X,1),1)];
                                    end
                                    %get K for low pass filtering:
                                    K = get_K(1:size(X,1),LPF.fwhm1,IOI.dev.TR);
                                    %filter X - LPF
                                    %calculate X inverse
                                    %Xm = pinv(X);
                                    %Xu = X(:,1);
                                    %filter forward
                                    X1 = ioi_filter_HPF_LPF_WMDL(K,X);
                                    %filter backward
%                                     X1 = X1(end:-1:1,:);
%                                     X1 = ioi_filter_HPF_LPF_WMDL(K,X1);
%                                     X1 = X1(end:-1:1,:);
                                    KX = spm_sp('Set', X1);
                                    KX.X = full(KX.X); %Filtered X
                                    Xm = spm_sp('x-',KX); % projector
                                    %covariance
                                    bcov = Xm * K.KL;
                                    bcov = bcov * bcov';
                                    %approximate calculation of effective degrees of freedom
                                    [trRV trRVRV] = approx_trRV(KX.X,Xm,K.KL);
%                                     if ~iscell(Xtmp)
%                                         cX = c1;
%                                     else
%                                         cX = 1;
%                                     end
                                    cX = c1;
                                    IOI.X{s1}.X{cX} = X;
                                    IOI.X{s1}.Xm{cX} = Xm;
                                    IOI.X{s1}.bcov{cX} = bcov;
                                    IOI.X{s1}.trRV{cX} = trRV;
                                    IOI.X{s1}.trRVRV{cX} = trRVRV;
                                    IOI.X{s1}.erdf{cX} = (trRV)^2/trRVRV;
                                    %************by Cong 
                                    if isfield(job.do_F_test,'F_test_enabled')
                                        if model_choice                                        
                                            X_spikes = [X(:,1) ones(size(X(:,1),1),1)];
                                            F = X_spikes;
                                        else
                                            X_stim = [X(:,2) ones(size(X(:,2),1),1)];
                                            F = X_stim;
                                        end  
                                     K = get_K(1:size(F,1),LPF.fwhm1,IOI.dev.TR);
                                    %filter X - LPF
                                    %calculate X inverse
                                    %Xm = pinv(X);
                                    %Xu = X(:,1);
                                    %filter forward
                                    F1 = ioi_filter_HPF_LPF_WMDL(K,F);
                                    %filter backward
%                                     X1 = X1(end:-1:1,:);
%                                     X1 = ioi_filter_HPF_LPF_WMDL(K,X1);
%                                     X1 = X1(end:-1:1,:);
                                    KF = spm_sp('Set', F1);
                                    KF.F = full(KF.X); %Filtered X
                                    Fm = spm_sp('x-',KF); % projector
                                    %covariance
                                    bcov = Fm * K.KL;
                                    bcov = bcov * bcov';
                                    %approximate calculation of effective degrees of freedom
                                    [trRV trRVRV] = approx_trRV(KF.F,Fm,K.KL);
%                                     if ~iscell(Xtmp)
%                                         cX = c1;
%                                     else
%                                         cX = 1;
%                                     end
                                    cX = c1;
                                    IOI.F{s1}.F{cX} = F;
                                    IOI.F{s1}.Fm{cX} = Fm;
                                    IOI.F{s1}.bcov{cX} = bcov;
                                    IOI.F{s1}.trRV{cX} = trRV;
                                    IOI.F{s1}.trRVRV{cX} = trRVRV;
                                    IOI.F{s1}.erdf{cX} = (trRV)^2/trRVRV;
                                    end
                                    %****************end
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
                                                if HPF.hpf_butter_On
                                                    y =  ButterHPF(1/IOI.dev.TR,HPF.hpf_butter_freq,HPF.hpf_butter_order,y);
                                                end
                                                %yu = y;
                                                %filtering of the data: LPF (Gaussian), forward
                                                y = ioi_filter_HPF_LPF_WMDL(K,y')';
                                                %filter backward
%                                                 y = y(end:-1:1);
%                                                 y = ioi_filter_HPF_LPF_WMDL(K,y')';
%                                                 y = y(end:-1:1);
                                                %GLM inversion: calculating beta and residuals - would be SPM.Vbeta,
                                                b = Xm * y'; % beta : least square estimate                                           
                                                %Compute t stat
                                                res = y'-X*b;
                                                res2 = sum(res.^2);
                                                IOI.X{s1}.b{r1,c1} = b;
                                                IOI.X{s1}.r(r1,c1) = res2; % Residuals
                                                IOI.X{s1}.mse(r1,c1) = res2/length(y);
                                                IOI.X{s1}.t(r1,c1) = b(1)/(res2*bcov(1,1)/trRV)^0.5;
                                                %%%%%%%%%%%%%%%du f-test by
                                                %%%%%%%%%%%%%%%cong on
                                                %%%%%%%%%%%%%%%12/12/19
                                                if isfield(job.do_F_test,'F_test_enabled')
                                                b_re = Fm * y';
                                                res_re = y'- F*b_re;
                                                res_re2 = sum(res_re.^2);
                                                IOI.F{s1}.b_re{r1,c1} = b_re;
                                                IOI.F{s1}.r_re(r1,c1) = res_re2; % Residuals
%                                                 IOI.F{s1}.mse_re(r1,c1) = res_re2/length(y);                                               
                                                degree_X = length(y)-size(X,2);
                                                F_value_numerator = (res_re2 - res2)/(size(X,2)-size(F,2));
                                                F_value_denominator = res2/degree_X;
                                                F_value = F_value_numerator/F_value_denominator;
                                                IOI.F{s1}.f(r1,c1) = F_value;
                                                
                                                end
                                                %****************end
                                                if volt==2
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
                    lp1{1} = ':'; lp1{2} = '-'; 
                    %%%%%%%%%%%%%%%%%by Cong on 2012_07_06
                    lp1{3} = '--'; lp1{4} = '-.';
                    %*********************************
                    ctotal = [];
                    if isfield(IOI.color,'HbO')
                        lp2{IOI.color.eng==IOI.color.HbO} = 'r'; %HbO
                        lp2{IOI.color.eng==IOI.color.HbR} = 'b'; %HbR                       
                        ctotal = [ctotal find(IOI.color.eng==IOI.color.HbO) ...
                            find(IOI.color.eng==IOI.color.HbR) ];
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
                                
                            %**************************************
                                c0 = length(ctotal)-0;
                            %**************************************************
                                    if all_ROIs || sum(r1==selected_ROIs)
                                        h1 = h1 + 1;
                                        h(h1) = figure;
                                        legstr = {};
                                        %loop over available colors
                                        %**************************************
                                            for c1=1:ctotal(c0)

                                        %**************************************
                                        %for c1=1:ctotal(); %length(IOI.sess_res{s1}.fname) %or ctotal
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
                                                            if isfield(job.do_F_test,'F_test_enabled')
                                                            f = IOI.F{s1}.f(r1,c1);
                                                            end
                                                            if volt==2, t2 = IOI.X{s1}.t2(r1,c1); end
                                                            %***********************
                                                            %plot(lp,yf,[lp1{1} lp2{c1}]); hold on
                                                            %*************************
                                                            plot(lp,yp,[lp1{2} lp2{c1}]); hold on
%                                                             if c1==6
%                                                             plot(lp,yp,[lp1{2} lp2{c1}]); hold on  
                                                            
                                                            %*********************************
%                                                             legstr = [legstr; [lp3{c1} 'f']];
                                                            %************************************
                                                            if show_mse
                                                                mse_str = [',' sprintf('%2.3f',IOI.X{s1}.mse(r1,c1))];
                                                            else
                                                                mse_str = '';
                                                            end

                                                            if volt == 1
                                                                if isfield(job.do_F_test,'F_test_enabled')
                                                                legstr = [legstr; [lp3{c1} 'p(' sprintf('%2.1f',t) mse_str ',' sprintf('%2.1f',f) ')']];
                                                                else
                                                                legstr = [legstr; [lp3{c1} 'p(' sprintf('%2.1f',t) mse_str ')']];
                                                                end
                                                            else
                                                                legstr = [legstr; [lp3{c1} 'p(' sprintf('%2.1f',t) ',' sprintf('%2.1f',t2) mse_str ')']];
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        %**********************************
                                    end
                                    %**************************************
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
                                %end
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
