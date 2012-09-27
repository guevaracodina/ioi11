function out = ioi_stim_mean_run(job)
[all_sessions selected_sessions] = ioi_get_sessions(job);
[all_ROIs selected_ROIs] = ioi_get_ROIs(job);
EBS = 3; %error bar step
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

extract_HRF = job.extract_HRF;
%Other options
generate_global = job.generate_global;
IC = job.IC; %colors to include
PGM = [];

%save_figures
save_figures = job.save_figures;
generate_figures = job.generate_figures;
include_nlinfit = job.include_nlinfit;
add_error_bars = job.add_error_bars;

for SubjIdx=1:length(job.IOImat)
    try
        tic
        clear ROI onsets_list M
        %Load IOI.mat information
        [IOI IOImat dir_ioimat]= ioi_get_IOI(job,SubjIdx);
        
        if ~isfield(IOI,'dev')
            IOI.dev.TR = 0.2;
        end
        
        window_after = round(job.window_after/IOI.dev.TR);
        %window_before = round(job.window_before/IOI.dev.TR);
        %window_offset = round(job.window_offset/IOI.dev.TR);
        if ~isfield(IOI.res,'seriesOK')
            disp(['No extracted time series available for subject ' int2str(SubjIdx) ' ... skipping series extraction']);
        else
            if ~isfield(IOI.res,'meanOK') || job.force_redo
                
                %load ROI
                if ~isempty(job.ROImat)
                    load(job.ROImat{SubjIdx});
                else
                    try
                        load(IOI.ROI.ROIfname);
                    catch
                        ts1 = strfind(dir_ioimat,'Res');
                        ts2 = strfind(dir_ioimat,filesep);
                        ts3 = ts2(ts2>ts1);
                        tdir = dir_ioimat(1:ts3(2));
                        load(fullfile(tdir,'ROI.mat'));
                    end
                end
                if save_figures
                    dir_fig = fullfile(dir_ioimat,'fig');
                    if ~exist(dir_fig,'dir'),mkdir(dir_fig);end
                end
                
                %get stimulation information - Careful, here onset duration is ignored!
                Ns = length(IOI.sess_res);
                if IC.include_HbT
                    if ~isfield(IOI.color,'HbT')
                        IOI.color.HbT = 'T';
                        IOI.color.eng = [IOI.color.eng IOI.color.HbT];
                    end
                end
                %restric onsets
                [IOI onsets_list pars_list] = ioi_restrict_onsets(IOI,job,rmi,ust);
                %                 onsets_list=onsets_list{2};
                %Check whether there is the same number of onset types in
                %each session; this is a
                %check if averaging over all selected sessions is sensible
                %in which case, global_M = 1 (M is the number of onset types)
                if generate_global
                    global_M = 1; found_first_M = 0;
                    for s1=1:Ns
                        if all_sessions || sum(s1==selected_sessions)
                            M{s1} = length(onsets_list{s1});
                            if ~found_first_M
                                PGM = M{s1};
                                found_first_M = 1;
                            else
                                if ~(M{s1}==PGM) %PGM = possible_global_M;
                                    %no global M possible
                                    global_M = 0;
                                end
                            end
                        end
                    end
                else
                    global_M = 0;
                end
                
                %find max M to loop over
                maxM = length(onsets_list{1});
                if Ns > 1
                    for s1=2:Ns
                        if all_sessions || sum(s1==selected_sessions)
                            if length(onsets_list{s1}) > maxM
                                maxM = length(onsets_list{s1});
                            end
                        end
                    end
                end
                                               
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Compute mean over stimulations
                [IOI U0] = ioi_stim_mean_call(job,IOI,ROI,maxM,global_M,...
                    PGM,onsets_list,pars_list);
                IOI.res.meanOK = 1;
                save(IOImat,'IOI');
                
                Ma = IOI.res.Ma; %mean of "after segments", by region, stimulus type, and by color and region
                Da = IOI.res.Da; %standard error
                Dma = IOI.res.Dma; %mean of standard error - single number
                Sa = IOI.res.Sa;
                %Define M structure for Expectation-Maximization algorithm
                %Fit difference of two gamma functions
                if extract_HRF
                    IOI = ioi_extract_HRF(job,IOI,ROI,maxM,global_M,Ma,GMa);
                    save(IOImat,'IOI');
                end
                %Global Figures
                ctotal = [];
                h1 = 0;
                if generate_figures || save_figures
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
                                                if extract_HRF
                                                    if include_nlinfit
                                                        plot(ls,F{r2,m1}{c1,s1}.yp,[lp1{2} lp2{c1}]); hold on
                                                        leg_str = [leg_str; [lp3{c1} '-NL']];
                                                    end
                                                    %                                                     plot(ls,H{r2,m1}{c1,s1}.yp,[lp1{3} lp2{c1}]); hold on
                                                    %                                                     leg_str = [leg_str; [lp3{c1} '-EM']];
                                                    %***************************************change by cong on 06/28
                                                    try
                                                        plot(ls,H{r2,m1}{c1,s1}.yp,[lp1{3} lp2{c1}]); hold on
                                                        leg_str = [leg_str; [lp3{c1} '-EM']];
                                                    catch
                                                    end
                                                    %*****************************    end
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
                                                plot(ls,Ma{r1,m1}{c1,s1}); %hold on
                                                %****************************************** changed by Cong on 06/28
                                                %                                                 try
                                                %                                                 plot(ls,Ma{r1,m1}{c1,s1}); %hold on
                                                %                                                 catch
                                                %                                                 end
                                                %**********************************************fini
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