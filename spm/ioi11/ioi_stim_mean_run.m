function out = ioi_stim_mean_run(job)
[all_sessions selected_sessions] = ioi_get_sessions(job);
[all_ROIs selected_ROIs] = ioi_get_ROIs(job);
%filters
HPF = ioi_get_HPF(job);
LPF = ioi_get_LPF(job);

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
if isfield(job,'remove_stims_SD')
    remove_stims_SD = job.remove_stims_SD;
else
    remove_stims_SD = 1;
end
remove_segment_drift = job.remove_segment_drift;
fit_3_gamma = job.fit_3_gamma;
extract_HRF = job.extract_HRF;
%Other options
generate_global = job.generate_global;
IC = job.IC; %colors to include

%save_figures
save_figures = job.save_figures;
generate_figures = job.generate_figures;
normalize_choice = job.normalize_choice;
include_nlinfit = job.include_nlinfit;
window_offset = job.window_offset;
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
        window_offset = round(window_offset/IOI.dev.TR);
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
                    dir_fig = fullfile(newDir,'fig');
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
                [IOI onsets_list] = ioi_restrict_onsets(IOI,job,rmi,ust);
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
                %Preallocate the arrays
                if all_ROIs
                    selected_ROIs = 1:length(ROI);
                end
                %Keep track later (for GLM) of which ROIs were selected
                IOI.res.mean.selected_ROIs = selected_ROIs;
                Nroi = length(selected_ROIs);
                %Nc = length(IOI.color.eng);
                
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
                Dma = cell(Nroi,maxM);
                Dmb = cell(Nroi,maxM);
                Rs = cell(Nroi,maxM);
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
                                Gkb = 0; %counter of segments before onsets
                                Gka = 0; %counter of segments after onsets
                                Gkb2 = 0; %counter of skipped segments before onsets
                                Gka2 = 0; %counter of skipped segments after onsets
                                %loop over sessions
                                for s1=1:length(IOI.sess_res)
                                    if all_sessions || sum(s1==selected_sessions)
                                        tmp_array_before = zeros(1,window_before);
                                        tmp_array_after = zeros(1,window_after);
                                        kb = 0; %counter of segments before onsets
                                        ka = 0; %counter of segments after onsets
                                        kb2 = 0; %counter of skipped segments before onsets
                                        ka2 = 0; %counter of skipped segments after onsets
                                        %loop over sessions
                                        try
                                            if IC.include_HbT
                                                if IOI.color.eng(c1) == IOI.color.HbT
                                                    tmp_d = ROI{r1}{s1,IOI.color.eng == IOI.color.HbO}+...
                                                        ROI{r1}{s1,IOI.color.eng == IOI.color.HbR};
                                                else
                                                    tmp_d = ROI{r1}{s1,c1};
                                                end
                                            else
                                                tmp_d = ROI{r1}{s1,c1};
                                            end
                                            normalize_flow = 0;
                                            if normalize_flow
                                            %normalize flow
                                            if isfield(IOI.color,'flow')
                                                if IOI.color.eng(c1)==IOI.color.flow
                                                    tmp_d = tmp_d/mean(tmp_d); %or median
                                                end
                                            end
                                            end
                                            remove_jumps = 0;
                                            if remove_jumps
                                            %remove jumps options:
                                            OP.Sb = 4; %number of standard deviations
                                            OP.Nr = 1/IOI.dev.TR; %number of points removed before and after
                                            OP.Mp = 10/IOI.dev.TR; %size of gaps to be filled
                                            OP.sf = 1/IOI.dev.TR; %sampling frequency
                                            OP.ubf = 1; %use Butterworth filter
                                            OP.bf = 0.01; %Butterworth HPF cutoff
                                            OP.bo = 2; %Butterworth order
                                            
                                            %if ~isempty(rmi{i0})
                                            if IOI.color.eng(c1) == IOI.color.flow
                                                OP.ubf = 0;
                                            else
                                                OP.ubf = 1;
                                            end
                                            tmp_d = ioi_remove_jumps(tmp_d,OP);
                                            %end
                                            end
                                        catch
                                            tmp_d = [];
                                        end
                                        
                                        if ~isempty(tmp_d)
                                            if HPF.hpf_butter_On
                                                tmp_d = ButterHPF(1/IOI.dev.TR,HPF.hpf_butter_freq,HPF.hpf_butter_order,tmp_d);
                                            end
                                            if LPF.lpf_gauss_On
                                                K = get_K(1:length(tmp_d),LPF.fwhm1,IOI.dev.TR);
                                                y = tmp_d;
                                                %forward
                                                y = ioi_filter_HPF_LPF_WMDL(K,y')';
                                                %backward
                                                %                                                 y = y(end:-1:1);
                                                %                                                 y = ioi_filter_HPF_LPF_WMDL(K,y')';
                                                %                                                 y = y(end:-1:1);
                                                tmp_d = y;
                                            end
                                            
                                            first_pass = 1;
                                            second_pass = 0;
                                            removeSeg = [];
                                            while first_pass
                                                
                                                Sb{r2,m1}{c1,s1} = [];
                                                Sa{r2,m1}{c1,s1} = [];
                                                %loop over onsets for that session
                                                if ~isempty(onsets_list{s1}{m1})
                                                    U = round(onsets_list{s1}{m1}/IOI.dev.TR)-window_offset; %in data points
                                                    U0{s1} = ioi_get_U(IOI,[],U,0,s1); %only for plotting stims
                                                    
                                                    for u1=1:length(U)
                                                        if ~any(u1==removeSeg)
                                                            clear tmp_median;
                                                            try
                                                                tmp1 = tmp_d(U(u1)-window_before:U(u1)-1);
                                                                switch normalize_choice
                                                                    case 1
                                                                        tmp_median = median(tmp1);
                                                                    case 2
                                                                        tmp_median = tmp1(end);
                                                                    case 3
                                                                        tmp_median = mean(tmp1);
                                                                end
                                                                tmp_array_before = tmp_array_before + tmp1-tmp_median;
                                                                Gtmp_array_before = Gtmp_array_before + tmp1-tmp_median;
                                                                kb = kb+1;
                                                                Gkb = Gkb+1;
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
                                                                    switch normalize_choice
                                                                        case {1,3}
                                                                            %use the same length as window before to estimate end value of segment
                                                                            tmp_end = mean(tmp1(end-window_before:end));
                                                                        case 2
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
                                                                Gka = Gka+1;
                                                                Sa{r2,m1}{c1,s1} = [Sa{r2,m1}{c1,s1};tmp1];
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
                                                        %
                                                    end
                                                    
                                                    Mb{r2,m1}{c1,s1} = tmp_array_before/kb; %global mean before
                                                    Ma{r2,m1}{c1,s1} = tmp_array_after/ka; %global mean after
                                                    Da{r2,m1}{c1,s1} = std(Sa{r2,m1}{c1,s1},0,1)/sqrt(ka); %SEM
                                                    Db{r2,m1}{c1,s1} = std(Sb{r2,m1}{c1,s1},0,1)/sqrt(kb); %SEM
                                                    Dma{r2,m1}{c1,s1} = mean(Da{r2,m1}{c1,s1});
                                                    Dmb{r2,m1}{c1,s1} = mean(Db{r2,m1}{c1,s1});
                                                    Rs{r2,m1}{c1,s1} = removeSeg;
                                                    
                                                    if remove_stims_SD
                                                        if second_pass
                                                            first_pass = 0;
                                                        else
                                                            %removeSeg = [];
                                                            if ~isempty(tmp_array_after)
                                                                meanA = mean(Ma{r2,m1}{c1,s1});
                                                                tmpSeg = Sa{r2,m1}{c1,s1};
                                                                meanSeg = mean(tmpSeg,2);
                                                                meanSd = std(tmpSeg(:));
                                                                for a0=1:length(meanSeg)
                                                                    if abs(meanSeg(a0)-meanA) > 1.5*meanSd %very strong criterion perhaps better to keep it at 2*meanSd
                                                                        removeSeg = [removeSeg a0];
                                                                    end
                                                                end
                                                            end
                                                            second_pass = 1;
                                                        end
                                                    else
                                                        first_pass = 0;
                                                    end
                                                else
                                                    first_pass = 0;
                                                end
                                            end
                                        end
                                    end
                                end
                                
                                if global_M && m1 <= possible_global_M
                                    GMb{r2,m1}{c1} = tmp_array_before/Gkb; %global mean before
                                    GMa{r2,m1}{c1} = tmp_array_after/Gka; %global mean after
                                    GDa{r2,m1}{c1} = std(GSa{r2,m1}{c1},0,1)/sqrt(Gka); %SEM
                                    GDb{r2,m1}{c1} = std(GSb{r2,m1}{c1},0,1)/sqrt(Gkb); %SEM
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
                IOI.res.Rs = Rs;
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
                                            try
                                                d = Ma{r2,m1}{c1,s1};
                                                if ~isempty(d)
                                                    F{r2,m1}{c1,s1} = ioi_nlinfit(x,d,IOI.color,c1,IC.include_flow);
                                                    H{r2,m1}{c1,s1} = ioi_HDM_hrf(IOI.dev.TR,x,d,IOI.color,c1,IC.include_flow,fit_3_gamma);
                                                end
                                            catch
                                                F{r2,m1}{c1,s1} = [];
                                                H{r2,m1}{c1,s1} = [];
                                            end
                                        end
                                    end
                                    if global_M
                                        try
                                            d = GMa{r2,m1}{c1};
                                            if ~isempty(d)
                                                GF{r2,m1}{c1} = ioi_nlinfit(x,d,IOI.color,c1,IC.include_flow);
                                                GH{r2,m1}{c1,s1} = ioi_HDM_hrf(IOI.dev.TR,x,d,IOI.color,c1,IC.include_flow,fit_3_gamma);
                                            end
                                        catch
                                            GF{r2,m1}{c1} = [];
                                            GH{r2,m1}{c1,s1} = [];
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
                                                        errorbar(ls(2:end-1),Ma{r2,m1}{c1,s1}(2:end-1),Da{r2,m1}{c1,s1}(2:end-1),[lp1{1} lp2{c1}]); hold on
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