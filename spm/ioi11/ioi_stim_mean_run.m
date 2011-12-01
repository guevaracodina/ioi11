function out = ioi_stim_mean_run(job)
%select a subset of sessions
if isfield(job.session_choice,'select_sessions')
    all_sessions = 0;
    selected_sessions = job.session_choice.select_sessions.selected_sessions;
else
    all_sessions = 1;
end
%select onsets
if isfield(job.stim_choice,'electro_stims')
    electro_stims = 1;
else
    %default stims
    electro_stims = 0;
end
%select nature of stimulation data to be used for averaging
if isfield(job.ROI_choice,'select_ROIs')
    all_ROIs = 0;
    selected_ROIs = job.ROI_choice.select_ROIs.selected_ROIs;
else
    all_ROIs = 1;
end
%HPF
if isfield(job.hpf_butter,'hpf_butter_On')
    HPF.hpf_butter_On = 1;
    HPF.hpf_butter_freq = job.hpf_butter.hpf_butter_On.hpf_butter_freq;
    HPF.hpf_butter_order = job.hpf_butter.hpf_butter_On.hpf_butter_order;
else
    HPF.hpf_butter_On = 0;
end
try 
    remove_segment_drift = job.remove_segment_drift;
catch
    remove_segment_drift = 1;
end
try
    fit_3_gamma = job.fit_3_gamma;
catch
    fit_3_gamma = 0;
end

%Other options
generate_global = job.generate_global;
include_flow = job.include_flow;
extract_HRF = job.extract_HRF;
%save_figures
save_figures = job.save_figures;
normalize_choice = job.normalize_choice;
%get size of windows before and after stimulation to average on, in data points

generate_figures = job.generate_figures;
add_error_bars = job.add_error_bars;
for SubjIdx=1:length(job.IOImat)
    try
        tic
        clear IOI ROI onsets_list M
        %Load IOI.mat information
        IOImat = job.IOImat{SubjIdx};
        load(IOImat);
        if ~isfield(IOI,'dev')
            IOI.dev.TR = 0.2;
        end
        window_after = round(job.window_after/IOI.dev.TR);
        window_before = round(job.window_before/IOI.dev.TR);
        
        if ~isfield(IOI.res,'seriesOK')
            disp(['No extracted time series available for subject ' int2str(SubjIdx) ' ... skipping series extraction']);
        else
            if ~isfield(IOI.res,'meanOK') || job.force_redo
                [dir_ioimat dummy] = fileparts(job.IOImat{SubjIdx});
                if isfield(job.IOImatCopyChoice,'IOImatCopy')
                    newDir = job.IOImatCopyChoice.IOImatCopy.NewIOIdir;
                    newDir = fullfile(dir_ioimat,newDir);
                    if ~exist(newDir,'dir'),mkdir(newDir); end
                    IOImat = fullfile(newDir,'IOI.mat');
                else
                    newDir = dir_ioimat;
                end
                if isfield(IOI.ROI,'ROIfname')
                    try
                        load(IOI.ROI.ROIfname);
                    catch
                        load(fullfile(dir_ioimat,'ROI.mat'));
                    end
                end
                if save_figures
                    dir_fig = fullfile(newDir,'fig');
                    if ~exist(dir_fig,'dir'),mkdir(dir_fig);end
                end
                
                %get stimulation information - Careful, here onset duration is ignored!
                Ns = length(IOI.sess_res);
                %loop over sessions
                for s1=1:Ns
                    if all_sessions || sum(s1==selected_sessions)
                        if ~electro_stims %default stims
                            onsets_list{s1} = IOI.sess_res{s1}.onsets;
                        else
                            onsets_list{s1} = {};
                            for i0=1:length(IOI.Sess(s1).U)
                                onsets_list{s1} = [onsets_list{s1}; IOI.Sess(s1).U(i0).ons];
                            end
                        end
                    end
                end
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
                %Preallocate the arrays
                if all_ROIs
                    selected_ROIs = 1:length(ROI);
                end
                %Keep track later (for GLM) of which ROIs were selected
                IOI.res.mean.selected_ROIs = selected_ROIs;
                Nroi = length(selected_ROIs);
                Nc = length(IOI.color.eng);
                
                if global_M
                    GSa = cell(Nroi,PGM);
                    GSb = cell(Nroi,PGM);
                    GMa = cell(Nroi,PGM);
                    GMb = cell(Nroi,PGM);
                    GDa = cell(Nroi,PGM);
                    GDb = cell(Nroi,PGM);
                    %for c1=1:Nc
                    
                end
                Sa = cell(Nroi,maxM);
                Sb = cell(Nroi,maxM);
                Ma = cell(Nroi,maxM);
                Mb = cell(Nroi,maxM);
                Da = cell(Nroi,maxM);
                Db = cell(Nroi,maxM);
                
                %loop over onset types
                for m1=1:maxM 
                    %loop over colors
                    for c1=1:length(IOI.color.eng)
                        %loop over ROIs
                        r2 = 0;
                        for r1=1:length(ROI)
                            if all_ROIs || sum(r1==selected_ROIs)
                                r2 = r2 + 1;
                                Gtmp_array_before = zeros(1,window_before);
                                Gtmp_array_after = zeros(1,window_after);
                                GSb{r2,m1}{c1} = [];
                                GSa{r2,m1}{c1} = [];
                                kb = 0; %counter of segments before onsets
                                ka = 0; %counter of segments after onsets
                                kb2 = 0; %counter of skipped segments before onsets
                                ka2 = 0; %counter of skipped segments after onsets
                                %loop over sessions
                                for s1=1:length(IOI.sess_res)
                                    if all_sessions || sum(s1==selected_sessions)
                                        tmp_array_before = zeros(1,window_before);
                                        tmp_array_after = zeros(1,window_after);
                                        try
                                            tmp_d = ROI{r1}{s1,c1};
                                            %normalize flow
                                            if IOI.color.eng(c1)==IOI.color.flow
                                                tmp_d = tmp_d/mean(tmp_d); %or median
                                            end
                                        catch
                                            tmp_d = [];
                                        end
                                        if ~isempty(tmp_d)
                                            if HPF.hpf_butter_On
                                                tmp_d = ButterHPF(1/IOI.dev.TR,HPF.hpf_butter_freq,HPF.hpf_butter_order,tmp_d);
                                            end
                                            Sb{r2,m1}{c1,s1} = [];
                                            Sa{r2,m1}{c1,s1} = [];
                                            %loop over onsets for that session
                                            U = round(onsets_list{s1}{m1}/IOI.dev.TR); %in data points
                                            for u1=1:length(U)
                                                clear tmp_median;
                                                try
                                                    tmp1 = tmp_d(U(u1)-window_before:U(u1)-1);
                                                    switch normalize_choice
                                                        case 1
                                                            tmp_median = median(tmp1);
                                                        case 2
                                                            tmp_median = tmp1(end);
                                                    end
                                                    tmp_array_before = tmp_array_before + tmp1-tmp_median;
                                                    Gtmp_array_before = Gtmp_array_before + tmp1-tmp_median;
                                                    kb = kb+1;
                                                    Sb{r2,m1}{c1,s1} = [Sb{r2,m1}{c1,s1};tmp1-tmp_median];
                                                    if global_M && m1 <= possible_global_M
                                                        GSb{r2,m1}{c1} = [GSb{r2,m1}{c1};tmp1-tmp_median];
                                                    end
                                                catch
                                                    kb2 = kb2+1;
                                                    if kb2 < 3 && r2 == 1
                                                        disp(['Could not include segment before onset ' int2str(u1) ' at time ' num2str(U(u1)*IOI.dev.TR) ...
                                                            ' for session ' int2str(s1) ' for ROI ' int2str(r1) ...
                                                            ' for color ' int2str(c1) ' for onset type ' int2str(m1) ...
                                                            ' in global average over all sessions... skipping ' int2str(kb2) ' so far']);
                                                    end
                                                end
                                                try
                                                    tmp1 = tmp_d(U(u1):U(u1)+window_after-1);
                                                    if ~exist('tmp_median','var') || normalize_choice == 2
                                                        tmp_median = tmp1(1);
                                                    end
                                                    if remove_segment_drift
                                                        if normalize_choice == 1
                                                            %use the same length as window before to estimate end value of segment
                                                            tmp_end = mean(tmp1(end-window_before:end));
                                                        else
                                                            tmp_end = tmp1(end);
                                                        end
                                                        slope = tmp_median + linspace(0,1,length(tmp1))*(tmp_end-tmp_median);
                                                        tmp1 = tmp1 - slope;
                                                    else
                                                        tmp1 = tmp1-tmp_median;
                                                    end
                                                    tmp_array_after = tmp_array_after + tmp1;
                                                    Gtmp_array_after = Gtmp_array_after + tmp1;
                                                    ka = ka+1;
                                                    Sa{r2,m1}{c1} = [Sa{r2,m1}{c1};tmp1];
                                                    if global_M && m1 <= possible_global_M
                                                        GSa{r2,m1}{c1} = [GSa{r2,m1}{c1};tmp1];
                                                    end
                                                catch
                                                    ka2 = ka2+1;
                                                    if ka2<3 && r2 == 1
                                                        disp(['Could not include segment after onset ' int2str(u1) ' at time ' num2str(U(u1)*IOI.dev.TR) ...
                                                            ' for session ' int2str(s1) ' for ROI ' int2str(r1) ...
                                                            ' for color ' int2str(c1) ' for onset type ' int2str(m1) ...
                                                            ' in global average over all sessions... skipping ' int2str(ka2) ' so far']);
                                                        
                                                    end
                                                end
                                                %                                                     if any(isnan(tmp_array_after))
                                                %                                                         a=1;
                                                %                                                     end
                                            end
                                            %                                             else
                                            %                                                 if ~isfield(IOI.color,'contrast') || (isfield(IOI.color,'contrast') && ~(IOI.color.eng(c1)==IOI.color.contrast))
                                            %                                                     disp(['Skipped session ' int2str(s1) ' for color ' IOI.color.eng(c1) ' for ROI ' int2str(r1)]);
                                            %                                                 end
                                            Mb{r2,m1}{c1,s1} = tmp_array_before/kb; %global mean before
                                            Ma{r2,m1}{c1,s1} = tmp_array_after/ka; %global mean after
                                            Da{r2,m1}{c1,s1} = std(Sa{r2,m1}{c1,s1},0,1)/sqrt(ka); %SEM
                                            Db{r2,m1}{c1,s1} = std(Sb{r2,m1}{c1,s1},0,1)/sqrt(kb); %SEM
                                            Dma{r2,m1}{c1,s1} = mean(Da{r2,m1}{c1,s1});
                                            Dmb{r2,m1}{c1,s1} = mean(Db{r2,m1}{c1,s1});
                                        end
                                    end
                                end
                                
                                if global_M && m1 <= possible_global_M
                                    GMb{r2,m1}{c1} = tmp_array_before/kb; %global mean before
                                    GMa{r2,m1}{c1} = tmp_array_after/ka; %global mean after
                                    GDa{r2,m1}{c1} = std(GSa{r2,m1}{c1},0,1)/sqrt(ka); %SEM
                                    GDb{r2,m1}{c1} = std(GSb{r2,m1}{c1},0,1)/sqrt(kb); %SEM
                                    GDma{r2,m1}{c1} = mean(GDa{r2,m1}{c1});
                                    GDmb{r2,m1}{c1} = mean(GDb{r2,m1}{c1});
                                end
                                if (ka2>0 || kb2 > 0) && r2 == 1
                                    disp(['Skipped ' int2str(ka2) ' segments after onsets and ' int2str(kb2) ' segments before onsets']);
                                end
                            end
                        end
                    end
                end
                
                IOI.res.meanOK = 1;
                %Global results
                if global_M
                    try
                        IOI.res.GMa = GMa; %mean of all segments
                        IOI.res.GMb = GMb;
                        IOI.res.GDa = GDa; %standard error
                        IOI.res.GDb = GDb;
                        IOI.res.GDma = GDma; %mean of standard error - single number
                        IOI.res.GDmb = GDmb;
                        IOI.res.GSa = GSa; %each segment
                        IOI.res.GSb = GSb;
                    end
                end
                %Session results
                IOI.res.Ma = Ma; %mean of "after segments", by region, stimulus type, and by color and region
                IOI.res.Mb = Mb; %mean of "before segments"
                IOI.res.Da = Da; %standard error
                IOI.res.Db = Db;
                IOI.res.Dma = Dma; %mean of standard error - single number
                IOI.res.Dmb = Dmb;
                IOI.res.Sa = Sa; %each segment
                IOI.res.Sb = Sb;
                save(IOImat,'IOI');
                %Define M structure for Expectation-Maximization algorithm                 
                %Fit difference of two gamma functions
                if extract_HRF                    
                    x = linspace(0,job.window_after-IOI.dev.TR,window_after);
                    %loop over onset types
                    for m1=1:maxM
                        %loop over colors
                        for c1=1:length(IOI.color.eng)
                            %loop over ROIs
                            r2 = 0;
                            for r1=1:length(ROI)
                                if all_ROIs || sum(r1==selected_ROIs)
                                    r2 = r2 + 1;
                                    %loop over sessions
                                    for s1=1:length(IOI.sess_res)
                                        if all_sessions || sum(s1==selected_sessions)
                                            d = Ma{r2,m1}{c1,s1};
                                            if ~isempty(d)
                                                F{r2,m1}{c1,s1} = ioi_nlinfit(x,d,IOI.color,c1,include_flow);
                                                H{r2,m1}{c1,s1} = ioi_HDM_hrf(IOI.dev.TR,x,d,IOI.color,c1,include_flow,fit_3_gamma);
                                            end
                                        end
                                    end
                                    if global_M
                                        d = GMa{r2,m1}{c1};
                                        if ~isempty(d)
                                            GF{r2,m1}{c1} = ioi_nlinfit(x,d,IOI.color,c1,include_flow);
                                            GH{r2,m1}{c1,s1} = ioi_HDM_hrf(IOI.dev.TR,x,d,IOI.color,c1,include_flow,fit_3_gamma);
                                        end
                                    end
                                end
                            end
                        end
                    end
                    if global_M
                        IOI.res.GF = GF; %mean of all segments
                        IOI.res.GH = GH;
                    end
                    %Session results
                    IOI.res.F = F;
                    IOI.res.H = H;
                    save(IOImat,'IOI');
                end
                %Global Figures
                ctotal = [];
                h1 = 0;
                if generate_figures || save_figures
                    %line specification - ROI x color
                    lp1{1} = '-'; lp1{2} = ':'; lp1{3} = '--'; lp1{4} = '-.';
                    if isfield(IOI.color,'HbO')
                        lp2{IOI.color.eng==IOI.color.HbO} = 'r'; %HbO
                        lp2{IOI.color.eng==IOI.color.HbR} = 'b'; %HbR
                        ctotal = [ctotal find(IOI.color.eng==IOI.color.HbO) ...
                            find(IOI.color.eng==IOI.color.HbR)];
                    end
                    if include_flow
                        if isfield(IOI.color,'flow')
                            lp2{IOI.color.eng==IOI.color.flow} = 'k'; %Flow
                            ctotal = [ctotal find(IOI.color.eng==IOI.color.flow)];
                        end
                    end
                    %global figures, only for 1st onset
                    if exist('GMa','var')
                        ls = linspace(0,job.window_after,window_after);
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
                            tit = 'Mean_ROI';
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
                                hc1 = set_colorbar(gcf,size(GMa,1));
                                legend off
                                set(gcf, 'Colormap', ColorSet);
                                tit = ['Mean_' IOI.color.eng(c1) '_allROI'];
                                title(tit);
                                ioi_save_figures(save_figures,generate_figures,h(h1),tit,dir_fig);
                            end
                        end
                    end
                    
                    %Figures by session and by stimulus type
                    ls = linspace(0,job.window_after,window_after);
                    for s1=1:length(IOI.sess_res)
                        if all_sessions || sum(s1==selected_sessions)
                            %loop over onset type
                            for m1=1:length(onsets_list{s1})
                                %Figures with colors combined on same figure - one figure per ROI
                                %loop over ROIs
                                r2 = 0;
                                for r1=1:length(ROI)
                                    if all_ROIs || sum(r1==selected_ROIs)
                                        r2 = r2+1;
                                        h1 = h1 + 1;
                                        h(h1) = figure;
                                        for c1 = ctotal
                                            if ~add_error_bars
                                                plot(ls,Ma{r2,m1}{c1,s1},[lp1{1} lp2{c1}]); hold on
                                            else
                                                errorbar(ls(2:end-1),Ma{r2,m1}{c1,s1}(2:end-1),Da{r2,m1}{c1,s1}(2:end-1),[lp1{1} lp2{c1}]); hold on
                                            end
                                            if extract_HRF
                                                plot(ls,F{r2,m1}{c1,s1}.yp,[lp1{2} lp2{c1}]); hold on
                                                plot(ls,H{r2,m1}{c1,s1}.yp,[lp1{3} lp2{c1}]); hold on
                                            end
                                        end
                                        if ~extract_HRF
                                        if include_flow
                                            legend(gca,'HbO','HbR','Flow');
                                        else
                                            legend(gca,'HbO','HbR');
                                        end
                                        else
                                        if include_flow
                                            legend(gca,'HbO','HbO-NL','HbO-EM','HbR','HbR-NL','HbR-EM','Flow','Flow-NL','Flow-EM');
                                        else
                                            legend(gca,'HbO','HbO-NL','HbO-EM','HbR','HbR-NL','HbR-EM');
                                        end   
                                        end
                                        tit = ['ROI ' int2str(r1) ', Session ' int2str(s1) ', Stimulus ' int2str(m1)];
                                        title(tit);
                                        ioi_save_figures(save_figures,generate_figures,h(h1),tit,dir_fig);
                                    end
                                end
                                
                                %Figures with ROIs combined on same figure - one figure per color
                                if size(Ma,1) <= length(lp1) %plot identifiable ROIs
                                    for c1 = ctotal
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
                                                leg = [leg; ['ROI ' int2str(r1)]];
                                            end
                                        end
                                        legend(gca,leg);
                                        tit = ['Color ' IOI.color.eng(c1) ', Session ' int2str(s1) ', Stimulus ' int2str(m1)];
                                        title(tit);
                                        ioi_save_figures(save_figures,generate_figures,h(h1),tit,dir_fig);
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
                                        end
                                        legend off
                                        set(gcf, 'Colormap', ColorSet);
                                        hc1 = set_colorbar(gcf,size(Ma,1));
                                        tit = ['Color ' IOI.color.eng(c1) ', Session ' int2str(s1) ', Stimulus ' int2str(m1)];
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
        toc
        disp(['Subject ' int2str(SubjIdx) ' complete']);
        out.IOImat{SubjIdx} = IOImat;
    catch exception
        disp(exception.identifier)
        disp(exception.stack(1))
    end
end
end

function hc1 = set_colorbar(gcf,len)
figure(gcf)
hc1 = colorbar;
set(hc1, 'YLim', [1 len+1]);
y_tick = linspace(1, len, len)'+0.49;
set(hc1, 'YTick', y_tick);
%set(hc1, 'YTickMode', 'Manual');
set(hc1, 'FontSize', 12);
%Customize here number of decimals
set(hc1,'YTickLabel',sprintf('%.0f |',get(hc1,'YTick')'));
end