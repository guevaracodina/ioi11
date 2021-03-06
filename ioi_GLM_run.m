function out = ioi_GLM_run(job)
%spm_shift = 32; %annoying shift in spm, present on X and U, etc.
%Volterra
volt = job.volt;
%bases
bases = job.bases;
[all_sessions selected_sessions] = ioi_get_sessions(job);
%filters
HPF = ioi_get_HPF(job);
%*********by Cong on 12 /11/08
% LPF = ioi_get_LPF(job);
%*****end
onset_choice=job.onset_choice;
if isfield(job,'remove_stims')
    rmi = job.remove_stims;
else
    rmi = '';
end
if isfield(job,'use_stims')
    ust = job.use_stims;
else
    ust = '';
end
%find data selection mode
if isfield(job.data_selection_choice,'ROI_mode')
    image_mode = 0;
    ROImat = job.data_selection_choice.ROI_mode.ROImat;
   job.ROI_choice = job.data_selection_choice.ROI_mode.ROI_choice;
    [all_ROIs selected_ROIs] = ioi_get_ROIs(job);
    show_mse = job.data_selection_choice.ROI_mode.show_mse;
    figure_show_stim = job.data_selection_choice.ROI_mode.figure_show_stim;
    figure_rebase_to_zero_at_stim = job.data_selection_choice.ROI_mode.figure_rebase_to_zero_at_stim;
else
    image_mode = 1;
    %Gaussian spatial LPF
    if isfield(job.data_selection_choice.image_mode.spatial_LPF,'spatial_LPF_On')
        radius = job.data_selection_choice.image_mode.spatial_LPF.spatial_LPF_On.spatial_LPF_radius;
        spatial_LPF = 1;
    else
        spatial_LPF = 0;
    end
    job.shrinkage_choice = job.data_selection_choice.image_mode.shrinkage_choice;
    [shrinkage_choice SH] = ioi_get_shrinkage_choice(job);
end


if isfield(job.vasomotion_choice,'vasomotion_on')
    vasomotion_on = 1;
else
    vasomotion_on = 0;
end
%Colorbar - to allow specifying common min and max on all charts
if isfield(job.override_colorbar,'colorbar_override')
    Im.cbar.c_min = job.override_colorbar.colorbar_override.colorbar_min;
    Im.cbar.c_max = job.override_colorbar.colorbar_override.colorbar_max;
    Im.cbar.colorbar_override = 1;
else
    Im.cbar.colorbar_override = 0;
end
save_figures = job.save_figures;
generate_figures = job.generate_figures;
save_beta_mse = job.save_beta_mse;
use_onset_amplitudes = job.use_onset_amplitudes;
IC = job.IC; %colors to include

%Big loop over subjects
for SubjIdx=1:length(job.IOImat)
    try
        tic
        clear IOI ROI SPM
        %Load IOI.mat information
        IOImat = job.IOImat{SubjIdx};
        load(IOImat);
        if ~isfield(IOI.res,'ROIOK') && ~image_mode
            disp(['No ROI available for subject ' int2str(SubjIdx) ' ... skipping series extraction']);
        else
            if ~isfield(IOI.res,'seriesOK') && ~image_mode
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
                    if IC.include_HbT
                        if ~isfield(IOI.color,'HbT')
                            IOI.color.HbT = 'T';
                            IOI.color.eng = [IOI.color.eng IOI.color.HbT];
                        end
                    end
                    %save shrunk images
                    if image_mode 
                        if shrinkage_choice
                            IOI = ioi_save_shrunk_images(IOI,job,SH,dir_ioimat);
                        end
                    end
                    %Overwrite old IOI
                    save(job.IOImat{SubjIdx},'IOI');
                    %loop over sessions
                    for s1=1:length(IOI.sess_res)
                        if all_sessions || sum(s1==selected_sessions)
                            %Electrophysiology, for each subject and session
                            
                            %TO-DO: generalize to more than 1 onset type
                            %************************by Cong on 12/11/08
                                                       switch onset_choice
                                case 0 %**onsets from stim and detection
                                    %onsets from detection
                                    ot = 1;
                                    ons{ot} = IOI.sess_res{s1}.onsets{ot}; %already in seconds *IOI.dev.TR;
                                    dur{ot} = IOI.sess_res{s1}.durations{ot}; %*IOI.dev.TR;
                                    name{ot} = IOI.sess_res{s1}.names{ot};
                                    if use_onset_amplitudes
                                        amp{ot} = IOI.sess_res{s1}.parameters{ot};
                                    else
                                        amp{ot} = [];
                                    end
                                    
                                    %*******************by cong on 12/11/05
                                    %onsets from stimulation.
                                    ot = 2;
                                    ons{ot} = IOI.sess_res{s1}.onsets{ot}; %already in seconds *IOI.dev.TR;
                                    dur{ot} = IOI.sess_res{s1}.durations{ot}; %*IOI.dev.TR;
                                    name{ot} = IOI.sess_res{s1}.names{ot};
                                    if use_onset_amplitudes
                                        amp{ot} = IOI.sess_res{s1}.parameters{ot};
                                    else
                                        amp{ot} = [];
                                    end
                                case 1 %***************onsets from detection
                                    ot = 1;
                                    ons = IOI.sess_res{s1}.onsets{ot}; %already in seconds *IOI.dev.TR;
                                    dur = IOI.sess_res{s1}.durations{ot}; %*IOI.dev.TR;
                                    name = IOI.sess_res{s1}.names{ot};
                                    if use_onset_amplitudes
                                        amp = IOI.sess_res{s1}.parameters{ot};
                                    else
                                        amp = [];
                                    end                                    
                                case 2   %*************onsets from stim
                                    ot = 2;
                                    ons= IOI.sess_res{s1}.onsets{ot}; %already in seconds *IOI.dev.TR;
                                    dur = IOI.sess_res{s1}.durations{ot}; %*IOI.dev.TR;
                                    name = IOI.sess_res{s1}.names{ot};
                                    if use_onset_amplitudes
                                        amp = IOI.sess_res{s1}.parameters{ot};
                                    else
                                        amp = [];
                                    end
                            end
                            %********end
% %                             ot = job.which_onset_type;
% %                             ons = IOI.sess_res{s1}.onsets{ot}; %already in seconds *IOI.dev.TR;
% %                             dur = IOI.sess_res{s1}.durations{ot}; %*IOI.dev.TR;
% %                             name = IOI.sess_res{s1}.names{ot};
% %                             if use_onset_amplitudes
% %                                 amp = IOI.sess_res{s1}.parameters{ot};
% %                             else
% %                                 amp = [];
% %                             end
                            %remove onsets
                            %Alternative to this ioi_remove_onsets: create
                            %a function from similar code in
                            %ioi_stim_mean_run and call this function
                            %before line just above that defines ons, dur,
                            %**********i donot know if it is removed. I
                            %think it should be removed because it has been
                            %removed in create onsets. 
% %                             [ons amp IOI] = ioi_remove_onsets(ons, amp, rmi, ust,IOI,s1,ot);
                            %convolve with hemodynamic response function
                            [Xtmp U] = ioi_get_X(IOI,name,ons,dur,amp,s1,bases,volt);
                            IOI.Sess(s1).U = U; %store onsets for each session
                            
                            %loop over available colors
                            for c1=1:length(IOI.color.eng)-1 %(IOI.sess_res{s1}.fname)
                                doColor = ioi_doColor(IOI,c1,IC);
                                if doColor
                                    %select design matrix
                                    if ~iscell(Xtmp)
                                        X = Xtmp;
                                    else
                                        X = Xtmp{c1};
                                    end
                                    
                                    if ~isempty(X)
%                                         if image_mode || vasomotion_on
%*************************by Cong on 12/11/08
                                        if image_mode 
                                            %********end 
                                            %put all the data for this color, and session, into memory
                                            %Note that this takes several minutes to load per session
                                            %Y is typically 3 GB or larger
                                            y = ioi_get_images(IOI,1:IOI.sess_res{s1}.n_frames,c1,s1,dir_ioimat,shrinkage_choice);
                                            %calculate vasomotion regressor
                                            [nx ny nt] = size(y);
                                            if vasomotion_on
                                                vaso = squeeze(mean(mean(y,2),1));
                                                X = [X vaso];
                                            end
                                        end
                                        IOI.X{s1}.X0 = X;
                                        %filter X - HPF
                                        if HPF.hpf_butter_On
                                            X = ButterHPF(1/IOI.dev.TR,HPF.hpf_butter_freq,HPF.hpf_butter_order,X);
                                        end
                                        %add a constant
                                        X = [X ones(size(X,1),1)];
                                        %get K for low pass filtering:
                                        LPF.fwhm1 = job.lpf_gauss.fwhm1;
                                        K = get_K(1:size(X,1),LPF.fwhm1,IOI.dev.TR);
                                        %filter X - LPF
                                        %calculate X inverse
                                        %Xm = pinv(X);
                                        %Xu = X(:,1);
                                        %filter forward
                                        X1 = ioi_filter_HPF_LPF_WMDL(K,X);
                                        %filter backward
                                        %                                         X1 = X1(end:-1:1,:);
                                        %                                         X1 = ioi_filter_HPF_LPF_WMDL(K,X1);
                                        %                                         X1 = X1(end:-1:1,:);
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
                                        
                                        
                                        if image_mode
                                            %set possible Inf values of Y to max of non Inf values of Y
                                            indInf = isinf(y(:));
                                            indNaN = isnan(y(:));
                                            y(indInf) = 0;
                                            maxY = max(y(:));
                                            y(indInf) = maxY;
                                            y(indNaN) = maxY;
                                            clear indNaN indInf
                                            if spatial_LPF
                                                Ks.k1 = nx;
                                                Ks.k2 = ny;
                                                Ks.radius = radius;
                                                Ks = ioi_spatial_LPF('set',Ks);
                                                %Gaussian spatial low pass filter
                                                for i1=1:nt
                                                    y(:,:,i1) = ioi_spatial_LPF('lpf',Ks,squeeze(y(:,:,i1)));
                                                end
                                            end
                                            %reshape
                                            y = reshape(y,nx*ny,nt);
                                            %GLM on whole images
                                            if HPF.hpf_butter_On
                                                y =  ButterHPF(1/IOI.dev.TR,HPF.hpf_butter_freq,HPF.hpf_butter_order,y);
                                            end
                                            
                                            %filtering of the data: LPF (Gaussian), forward
                                            y = ioi_filter_HPF_LPF_WMDL(K,y')'; %takes about 2 minutes
                                            %filter backward
                                            %                                             y = y(:,end:-1:1);
                                            %                                             y = ioi_filter_HPF_LPF_WMDL(K,y')';
                                            %                                             y = y(:,end:-1:1);
                                            %                                             if any(isnan(y(:))) || any(isinf(y(:)))
                                            %                                                 a0 =1;
                                            %                                             end
                                            %GLM inversion: calculating beta and residuals - would be SPM.Vbeta,
                                            b = Xm * y'; % beta : least square estimate
                                            %Compute t stat
                                            %res = ; %a few seconds
                                            %ytilde = y'-X*b;
                                            res2 = sum((y'-X*b).^2);
                                            clear y
                                            mse = res2/length(res2);
                                            t = b(1,:)./(res2*bcov(1,1)/trRV).^0.5;
                                            %                                             if any(isnan(t(:))) || any(isinf(t(:)))
                                            %                                                 a0 =1;
                                            %                                             end
                                            if volt == 2
                                                t2 = b(2,:)./(res2*bcov(2,2)/trRV).^0.5;
                                            end
                                            %reshape
                                            b = reshape(b,[size(b,1) nx ny]);
                                            res2 = reshape(res2,[nx ny]);
                                            mse = reshape(mse,[nx ny]);
                                            t = reshape(t,[nx ny]);
                                            if volt == 2
                                                t2 = reshape(t2,[nx ny]);
                                            end
                                            Im.save_figures = save_figures;
                                            Im.generate_figures = generate_figures;
                                            Im0 = Im; %for images other than t-stats
                                            Im0.cbar.colorbar_override = 0;
                                            %save as nifti and store path to data
                                            name_str = ['S' gen_num_str(s1,2) '_' IOI.color.eng(c1)];
                                            name_title = ['S' int2str(s1) ' ' IOI.color.eng(c1)];
                                            for i0=1:size(b,1)
                                                b_names{i0} = fullfile(newDir,['beta' int2str(i0) '_' name_str]);
                                                b_title{i0} = ['beta' int2str(i0) ' ' name_title];
                                            end
                                            mse_name = fullfile(newDir,['mse_' name_str]);
                                            mse_title = ['mse ' name_title];
                                            t_name =  fullfile(newDir,['t_' name_str]);
                                            t_title = ['t ' name_title];
                                            if volt == 2
                                                t2_name =  fullfile(newDir,['t2_' name_str]);
                                                t2_title = ['t2 ' name_title];
                                            end
                                            vx = [1 1 1];
                                            %sign for mask -- no longer used
                                            %                                             try
                                            %                                                 if IOI.color.eng(c1) == IOI.color.red || IOI.color.eng(c1) == IOI.color.HbR
                                            %                                                     sgn = -1;
                                            %                                                 else
                                            %                                                     sgn = 1;
                                            %                                                 end
                                            %                                             catch
                                            %                                                 sgn = 1;
                                            %                                             end
                                            sgn =1;
                                            thold = 1.95; %threshold on t-stats
                                            if save_beta_mse
                                                for i0=1:size(b,1)
                                                    ioi_save_images(squeeze(b(i0,:,:)),b_names{i0},vx,Im0,b_title{i0});
                                                end
                                                ioi_save_images(mse,mse_name,vx,Im0,mse_title);
                                            end
                                            ioi_save_images(t,t_name,vx,Im,t_title);
                                            if save_beta_mse
                                                ioi_save_masked_images(t,t_name,vx,Im0,t_title,sgn,thold);
                                            end
                                            if volt == 2
                                                ioi_save_images(t2,t2_name,vx,Im,t2_title);
                                                if save_beta_mse
                                                    ioi_save_masked_images(t2,t2_name,vx,Im0,t2_title,sgn,thold);
                                                end
                                            end
                                            
                                            IOI.X{s1}.b{c1} = b_names;
                                            %IOI.X{s1}.r(c1) = res_name; % Residuals
                                            IOI.X{s1}.mse{c1} = mse_name; %res2/length(y);
                                            IOI.X{s1}.t{c1} = t_name;
                                            if volt == 2
                                                IOI.X{s1}.t2{c1} = t2_name;
                                            end
                                            %IOI.X{s1}.yf{c1} = y; %filtered data
                                            %IOI.X{s1}.yu{r1,c1} = yu;
                                            %IOI.X{s1}.yp{c1} = X * b;
                                        else
                                            %load ROI
                                            if ~isempty(ROImat)
                                                load(ROImat{SubjIdx});
                                            else
                                                try
                                                    load(IOI.ROI.ROIfname);
                                                catch
                                                    load(fullfile(dir_ioimat,'ROI.mat'));
                                                end
                                            end
                                            %GLM on ROIs
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
                                                        %                                                         y = y(end:-1:1);
                                                        %                                                         y = ioi_filter_HPF_LPF_WMDL(K,y')';
                                                        %                                                         y = y(end:-1:1);
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
                        end
                    end
                    IOI.job = job;
                    IOI.res.GLMOK = 1;
                    save(IOImat,'IOI');
                    if image_mode
                        %What to plot and save
                        
                    else
                        if ~exist('ROI','var')
                            %load ROI
                            if ~isempty(ROImat)
                                load(ROImat{SubjIdx});
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
                        if IC.include_flow
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
                                        for c1=1:ctotal; %length(IOI.sess_res{s1}.fname) %or ctotal
                                            try
                                                if IC.include_flow || ~(IOI.color.eng(c1)==IOI.color.flow)
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