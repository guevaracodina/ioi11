function out = ioi_GLM_run(job)
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
%find data selection mode
if isfield(job.data_selection_choice,'ROI_mode')
    image_mode = 0;
    ROImat = job.data_selection_choice.ROI_mode.ROImat;
    %select a subset of ROIs
    if isfield(job.data_selection_choice.ROI_mode.ROI_choice,'select_ROIs')
        all_ROIs = 0;
        selected_ROIs = job.data_selection_choice.ROI_mode.ROI_choice.select_ROIs.selected_ROIs;
    else
        all_ROIs = 1;
    end
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
end
if isfield(job.shrinkage_choice,'configuration_shrink')
    shrink_x = job.shrinkage_choice.configuration_shrink.shrink_x;
    shrink_y = job.shrinkage_choice.configuration_shrink.shrink_y;
    shrinkage_choice = 1;
else
    shrinkage_choice = 0;
end
save_figures = job.save_figures;
generate_figures = job.generate_figures;
use_onset_amplitudes = job.use_onset_amplitudes;
include_flow = job.include_flow;
include_OD = job.include_OD;
include_HbT = job.include_HbT;
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
                    if include_HbT
                        if ~isfield(IOI.color,'HbT')
                            IOI.color.HbT = 'T';
                            IOI.color.eng = [IOI.color.eng IOI.color.HbT];
                        end
                    end
                    %save shrunk images
                    if shrinkage_choice
                        for s1=1:length(IOI.sess_res)
                            if all_sessions || sum(s1==selected_sessions)
                                for c1=1:length(IOI.color.eng)
                                    doColor = ioi_doColor(IOI,c1,include_OD,include_flow,include_HbT);
                                    fname0 = {};
                                    if doColor                                        
                                        if IOI.color.eng(c1) == 'T'
                                            doHbT = 1;
                                            [cHbR cHbO] = ioi_find_HbRHbO(IOI,s1);
                                            fname = IOI.sess_res{s1}.fname{cHbR};
                                            fname2 = IOI.sess_res{s1}.fname{cHbO};
                                        else
                                            doHbT = 0;
                                            fname = IOI.sess_res{s1}.fname{c1};                                            
                                        end
                                        for f1=1:length(fname)
                                            V = spm_vol(fname{f1});
                                            Y = spm_read_vols(V);
                                            if doHbT
                                                V2 = spm_vol(fname2{f1});
                                                Y2 = spm_read_vols(V2);
                                                Y = Y+Y2;
                                            end
                                            %shrink by averaging
                                            Y0 = zeros(size(Y(1:shrink_x:(end-shrink_x+1),1:shrink_y:(end-shrink_y+1),:)));
                                            for i1=1:shrink_x
                                                for i2=1:shrink_y
                                                    Y0 = Y0 + Y(i1:shrink_x:(end-shrink_x+i1),i2:shrink_y:(end-shrink_y+i2),:);
                                                end
                                            end
                                            %save images
                                            [dir0 fil0 ext0] = fileparts(fname{f1});
                                            if doHbT
                                                fil0 = regexprep(fil0, ['_' IOI.color.HbR '_'],  ['_' IOI.color.HbT '_']);
                                            end
                                            tn = fullfile(dir0,[fil0 '_shrunk_' int2str(shrink_x) 'x' int2str(shrink_y) ext0]);
                                            fname0 = [fname0; tn]; 
                                            ioi_save_nifti(Y0,tn,[1 1 1]);
                                        end
                                    end
                                    IOI.sess_shrunk{s1}.fname{c1} = fname0;
                                end
                            end
                        end
                    end
                    %Overwrite old IOI
                    save(job.IOImat{SubjIdx},'IOI');
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
                            for c1=1:length(IOI.color.eng) %(IOI.sess_res{s1}.fname)
                                doColor = ioi_doColor(IOI,c1,include_OD,include_flow,include_HbT);                                       
                                if doColor
                                    %select design matrix
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
                                        X1 = ioi_filter_HPF_LPF_WMDL(K,X);
                                        %filter backward
                                        X1 = X1(end:-1:1,:);
                                        X1 = ioi_filter_HPF_LPF_WMDL(K,X1);
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
                                        
                                        if image_mode
                                            %put all the data for this color, and session, into memory
                                            %Note that this takes several minutes to load per session
                                            %Y is typically 3 GB or larger
                                            y = ioi_get_images(IOI,1:IOI.sess_res{s1}.n_frames,c1,s1,dir_ioimat);
                                            %set possible Inf values of Y to max of non Inf values of Y
                                            
                                            [nx ny nt] = size(y);
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
                                            if hpf_butter_On
                                                y =  ButterHPF(1/IOI.dev.TR,hpf_butter_freq,hpf_butter_order,y);
                                            end
                                            
                                            %filtering of the data: LPF (Gaussian), forward
                                            y = ioi_filter_HPF_LPF_WMDL(K,y')'; %takes about 2 minutes
                                            %filter backward
                                            y = y(:,end:-1:1);
                                            y = ioi_filter_HPF_LPF_WMDL(K,y')';
                                            y = y(:,end:-1:1);
%                                             if any(isnan(y(:))) || any(isinf(y(:)))
%                                                 a0 =1;
%                                             end
                                            %GLM inversion: calculating beta and residuals - would be SPM.Vbeta,
                                            b = Xm * y'; % beta : least square estimate
                                            %Compute t stat
                                            %res = ; %a few seconds
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
                                            for i0=1:size(b,1)
                                                ioi_save_images(squeeze(b(i0,:,:)),b_names{i0},vx,Im,b_title{i0});
                                            end
                                            ioi_save_images(mse,mse_name,vx,Im,mse_title);
                                            ioi_save_images(t,t_name,vx,Im,t_title);
                                            ioi_save_masked_images(t,t_name,vx,Im,t_title,sgn,thold);
                                            if volt == 2
                                                ioi_save_images(t2,t2_name,vx,Im,t2_title);
                                                ioi_save_masked_images(t2,t2_name,vx,Im,t2_title,sgn,thold);
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
                                                        if hpf_butter_On
                                                            y =  ButterHPF(1/IOI.dev.TR,hpf_butter_freq,hpf_butter_order,y);
                                                        end
                                                        %yu = y;
                                                        %filtering of the data: LPF (Gaussian), forward
                                                        y = ioi_filter_HPF_LPF_WMDL(K,y')';
                                                        %filter backward
                                                        y = y(end:-1:1);
                                                        y = ioi_filter_HPF_LPF_WMDL(K,y')';
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
                        end
                    end
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
                                        for c1=1:ctotal; %length(IOI.sess_res{s1}.fname) %or ctotal
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

