function ioi_stim_figures(job,IOI,onsets_list,ROI,dir_fig)
[all_sessions selected_sessions] = ioi_get_sessions(job);
[all_ROIs selected_ROIs] = ioi_get_ROIs(job);
EBS = 3; %error bar step
add_error_bars = job.add_error_bars;
window_after = round(job.window_after/IOI.dev.TR);
Ma = IOI.res.Ma; %mean of "after segments", by region, stimulus type, and by color and region
%GMa = IOI.res.GMa;
Da = IOI.res.Da; %standard error
%Dma = IOI.res.Dma; %mean of standard error - single number
%Sa = IOI.res.Sa;
IC = job.IC; %colors to include
save_figures = job.save_figures;
generate_figures = job.generate_figures;

%Global Figures
ctotal = [];
h1 = 0;
%line specification - ROI x color
lp1{1} = '-'; lp1{2} = ':'; lp1{3} = '--'; lp1{4} = '-.';
if isfield(IOI.color,'HbO') && IC.include_HbO
    lp2{IOI.color.eng==IOI.color.HbO} = 'r'; %HbO
    lp3{IOI.color.eng==IOI.color.HbO} = 'HbO'; %HbO
    ctotal = [ctotal find(IOI.color.eng==IOI.color.HbO)];
end
if isfield(IOI.color,'HbR') && IC.include_HbR
    lp2{IOI.color.eng==IOI.color.HbR} = 'b'; %HbR
    lp3{IOI.color.eng==IOI.color.HbR} = 'HbR'; %HbR
    ctotal = [ctotal find(IOI.color.eng==IOI.color.HbR)];
end
if IC.include_HbT
    lp2{IOI.color.eng==IOI.color.HbT} = 'g';
    lp3{IOI.color.eng==IOI.color.HbT} = 'HbT'; %HbT
    ctotal = [ctotal find(IOI.color.eng==IOI.color.HbT)];
end
if isfield(IOI.color,'flow') && IC.include_flow
    lp2{IOI.color.eng==IOI.color.flow} = 'k'; %Flow
    lp3{IOI.color.eng==IOI.color.flow} = 'Flow';
    ctotal = [ctotal find(IOI.color.eng==IOI.color.flow)];
end
if IC.include_OD
    %to do...
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%global figures, only for 1st onset
if exist('GMa','var')
    ls = linspace(-job.window_offset,job.window_after-job.window_offset,window_after);
    if length(GMa) <= length(lp1) %plot identifiable series
        h1 = h1 + 1;
        h(h1) = figure;
        for r1 = 1:size(GMa,1)
            for c1 = ctotal
                if ~add_error_bars
                    plot(ls,GMa{r1,1}{c1},[lp1{r1} lp2{c1}]); hold on
                else
                    if r1 == 1
                        errorbar(ls(2:end-1),GMa{r1,1}{c1}(2:end-1),GDa{r1,1}{c1}(2:end-1),[lp1{r1} lp2{c1}]); hold on
                    else
                        plot(ls,GMa{r1,1}{c1},[lp1{r1} lp2{c1}]); hold on
                    end
                end
            end
        end
        tit = [IOI.subj_name ' Mean ROI'];
        title(tit);
        ioi_save_figures(save_figures,generate_figures,h(h1),tit,dir_fig);
    else %plot all series with random colors, skip error bars
        ColorSet = varycolor(size(GMa,1));
        for c1 = ctotal
            h1 = h1 + 1;
            h(h1) = figure;
            set(gca, 'ColorOrder', ColorSet);
            hold all
            for r1=1:size(GMa,1)
                plot(ls,GMa{r1,1,c1}); %hold on
            end
            hc1 = ioi_set_colorbar(gcf,size(GMa,1));
            legend off
            set(gcf, 'Colormap', ColorSet);
            
            if save_figures
                tit = [IOI.subj_name ' Mean ' IOI.color.eng(c1) ' allROI'];
                title(tit);
                ioi_save_figures(save_figures,generate_figures,h(h1),tit,dir_fig);
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figures by session and by stimulus type
ls = linspace(-job.window_offset,job.window_after-job.window_offset,window_after);
for s1=1:length(IOI.sess_res)
    if all_sessions || sum(s1==selected_sessions)
        %loop over onset type
        for m1=1:length(onsets_list{s1})
            if any(m1==job.which_onset_type)
                %Figures with colors combined on same figure - one figure per ROI
                %loop over ROIs
                r2 = 0;
                for r1=1:length(ROI)
                    if all_ROIs || sum(r1==selected_ROIs)
                        r2 = r2+1;
                        h1 = h1 + 1;
                        h(h1) = figure;
                        leg_str = {};
                        for c1 = ctotal
                            try %in case we do not have this color
                                if ~add_error_bars
                                    plot(ls,Ma{r2,m1}{c1,s1},[lp1{1} lp2{c1}]); hold on
                                else
                                    errorbar(ls(2:EBS:end-1),Ma{r2,m1}{c1,s1}(2:EBS:end-1),Da{r2,m1}{c1,s1}(2:EBS:end-1),[lp1{1} lp2{c1}]); hold on
                                end
                                leg_str = [leg_str; lp3{c1}];
                            end
                            if job.extract_HRF
                                if job.include_nlinfit
                                    plot(ls,F{r2,m1}{c1,s1}.yp,[lp1{2} lp2{c1}]); hold on
                                    leg_str = [leg_str; [lp3{c1} '-NL']];
                                end
                                try
                                    plot(ls,H{r2,m1}{c1,s1}.yp,[lp1{3} lp2{c1}]); hold on
                                    leg_str = [leg_str; [lp3{c1} '-EM']];
                                catch
                                end
                            end
                        end
                        legend(leg_str);
                        if save_figures
                            if isfield(IOI.res.ROI{r1},'name')
                                tit = [IOI.subj_name ' ROI ' int2str(r1) ' ' IOI.res.ROI{r1}.name ', Session ' int2str(s1) ', Stimulus ' int2str(m1)];
                            else
                                tit = [IOI.subj_name ' ROI ' int2str(r1) ' ' int2str(r1) ', Session ' int2str(s1) ', Stimulus ' int2str(m1)];
                            end
                            title(tit);
                            ioi_save_figures(save_figures,generate_figures,h(h1),tit,dir_fig);
                        end
                        
                        %Generate a figure for each stim
                        lsp{1} = '-r'; lsp{2} = '-b'; lsp{3} = '-g';
                        if job.make_eachStim_figures
                            d = IOI.res.Sa{r1,m1}{ctotal(1),1};
                            [nst n] = size(d);
                            for sn1=1:nst
                                h = figure; hold on
                                lp0 = linspace(0,n*IOI.dev.TR,n);
                                cc1 = 0;
                                for c1 = ctotal %[6 7 9]
                                    cc1 = cc1+1;
                                    d = IOI.res.Sa{r1,m1}{c1,1};
                                    try
                                        plot(lp0,d(sn1,:),lsp{cc1});
                                    catch
                                        disp(['Problem: ' IOI.subj_name ', Session ' int2str(s1) ', ROI ' int2str(r1) ', Stim type ' int2str(m1) ', Color: ' IOI.color.eng(c1) ', Stim number:' gen_num_str(sn1,2)])
                                    end
                                end
                                tit = [IOI.subj_name ', Session ' int2str(s1) ', ROI ' int2str(r1) ', Stim type ' int2str(m1) ', Stim number ' gen_num_str(sn1,2)];
                                title(tit)
                                xlabel('Time (s)');
                                ylabel('Changes (%)');
                                pathGF4 = fullfile(dir_fig,'StimEach');
                                if ~exist(pathGF4,'dir'), mkdir(pathGF4); end
                                ioi_save_figures(save_figures,generate_figures,h,tit,pathGF4);
                            end
                        end
                    end
                end
                
                %Generate a combined figure for all stims
                if job.make_stim_figures
                    for c1 = ctotal
                        for r1=1:length(ROI)
                            if all_ROIs || sum(r1==selected_ROIs)
                                d = IOI.res.Sa{r1,m1}{c1,1};
                                [nst n] = size(d);
                                lp0 = linspace(0,n*IOI.dev.TR,n);
                                h = figure; plot(lp0,d);
                                tit = [IOI.subj_name ', Session ' int2str(s1) ', Stim type ' int2str(m1) ', Color ' IOI.color.eng(c1) ', ROI ' int2str(r1)];
                                title(tit)
                                xlabel('Time (s)');
                                ylabel('Changes (%)');
                                pathGF3 = fullfile(dir_fig,'StimGroup');
                                if ~exist(pathGF3,'dir'), mkdir(pathGF3); end
                                ioi_save_figures(save_figures,generate_figures,h,tit,pathGF3);
                            end
                        end
                    end
                end
                
                %Figures with ROIs combined on same figure - one figure per color
                if size(Ma,1) <= length(lp1) %plot identifiable ROIs
                    for c1 = ctotal
                        try %if we do not have this color
                            h1 = h1 + 1;
                            h(h1) = figure;
                            r2 = 0;
                            leg = {};
                            for r1=1:length(ROI)
                                if all_ROIs || sum(r1==selected_ROIs)                                    
                                    r2 = r2+1;
                                    if ~add_error_bars
                                        plot(ls,Ma{r2,m1}{c1,s1},[lp1{r2} lp2{c1}]); hold on
                                    else
                                        errorbar(ls(2:end-1),Ma{r2,m1}{c1,s1}(2:end-1),Da{r2,m1}{c1,s1}(2:end-1),[lp1{r2} lp2{c1}]); hold on
                                    end
                                    try
                                        leg = [leg; IOI.res.ROI{r1}.name];
                                    catch
                                        leg = [leg; ['ROI ' int2str(r1)]];
                                    end
                                end
                            end
                            legend(gca,leg);
                            if save_figures
                                tit = [IOI.subj_name ' Color ' IOI.color.eng(c1) ', Session ' int2str(s1) ', Stimulus ' int2str(m1)];
                                title(tit);
                                ioi_save_figures(save_figures,generate_figures,h(h1),tit,dir_fig);
                            end
                        end
                    end
                else
                    %plot all series with random colors, skip error bars
                    ColorSet = varycolor(size(Ma,1));
                    for c1 = ctotal
                        h1 = h1 + 1;
                        h(h1) = figure;
                        set(gca, 'ColorOrder', ColorSet);
                        hold all
                        for r1=1:size(Ma,1)
                            plot(ls,Ma{r1,m1}{c1,s1});
                        end
                        legend off
                        set(gcf, 'Colormap', ColorSet);
                        hc1 = ioi_set_colorbar(gcf,size(Ma,1));
                        if save_figures
                            tit = [IOI.subj_name ' Color ' IOI.color.eng(c1) ', Session ' int2str(s1) ', Stimulus ' int2str(m1)];
                            title(tit);
                            ioi_save_figures(save_figures,generate_figures,h(h1),tit,dir_fig);
                        end
                    end
                end
            end
        end
        lsp{1} = '.c'; %lsp{2} = '.m'; lsp{3} = '.g'; 
        %remove jumps options:
        TR = IOI.dev.TR;
        OP.Sb = 4; %number of standard deviations
        OP.Nr = 1/TR; %number of points removed before and after
        OP.Mp = 10/TR; %size of gaps to be filled
        OP.sf = 1/TR; %sampling frequency
        OP.ubf = 1; %use Butterworth filter
        OP.bf = 0.01; %Butterworth HPF cutoff
        OP.bo = 2; %Butterworth order
        %onsets
        armonsets = [];
        %TT = IOI.sess_res{s1}.n_frames*TR; %minimum: 14.33 min
        Stm = job.which_onset_type;
        onsets = IOI.sess_res{s1}.onsets{Stm};
        %durations = IOI.sess_res{s1}.durations{Stm}; %small variability, but why?
        %nonsets = length(onsets); %OK, au moins 10 stims de chaque type
        %tmonsets = onsets(end);
        try
            rmonsets = IOI.onsets_removed{s1}{Stm};
            armonsets = [armonsets rmonsets];
        end
        armonsets = sort(armonsets);
        %lonsets = max(tmonsets(:));
        %dlonsets = TT-lonsets; %varies from 0 to 22 seconds
        
        stm_len1 = 10;
        stm_len2 = 2*stm_len1;
        stm_len3 = stm_len1/2;
        if job.make_timeCourse_figures
            for c1 = ctotal
                for r1=1:length(ROI)
                    if all_ROIs || sum(r1==selected_ROIs)
                        if IOI.color.eng(c1) == IOI.color.HbT
                            tmp_d = ROI{r1}{s1,IOI.color.eng == IOI.color.HbO}+...
                                        ROI{r1}{s1,IOI.color.eng == IOI.color.HbR};
                        else
                            tmp_d = ROI{r1}{s1,c1};
                        end
                        Ns = length(tmp_d);
                        lp = linspace(0,Ns*TR,Ns);
                        %add a HPF
                        rmi = job.remove_stims;
                        if ~isempty(rmi)
                            if IOI.color.eng(c1) == IOI.color.flow;
                                OP.ubf = 0;
                            else
                                OP.ubf = 1;
                            end
                            tmp_d2 = ioi_remove_jumps(tmp_d,OP);
                        else
                            tmp_d2 = tmp_d;
                        end
                        d = ButterHPF(1/TR,OP.bf,OP.bo,tmp_d2);
                        
                        h = figure; plot(lp,d); hold on; plot(lp,tmp_d,'y'); plot(lp,tmp_d2,'g')
                        try
                            stem(rmi,stm_len1*ones(1,length(rmi)),'xk');
                            stem(armonsets,stm_len2*ones(1,length(armonsets)),'+r');
                        end
                        stem(onsets,stm_len3*ones(1,length(onsets)),lsp{1});
                        
                        tit = [IOI.subj_name ', Session ' int2str(s1) ', Color ' IOI.color.eng(c1) ', ROI ' int2str(r1)];
                        title(tit);
                        pathGF1 = fullfile(dir_fig,'TimeCourses');
                        if ~exist(pathGF1,'dir'), mkdir(pathGF1); end
                        ioi_save_figures(save_figures,generate_figures,h,tit,pathGF1);
                    end
                end
            end
        end
    end
end