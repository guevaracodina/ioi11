function out = ioi_stim_mean_image_run(job)
%select a subset of sessions
if isfield(job.session_choice,'select_sessions')
    all_sessions = 0;
    selected_sessions = job.session_choice.select_sessions.selected_sessions;
else
    all_sessions = 1;
end

%HPF
if isfield(job.hpf_butter,'hpf_butter_On')
    HPF.hpf_butter_On = 1;
    HPF.hpf_butter_freq = job.hpf_butter.hpf_butter_On.hpf_butter_freq;
    HPF.hpf_butter_order = job.hpf_butter.hpf_butter_On.hpf_butter_order;
else
    HPF.hpf_butter_On = 0;
end
%LPF
if isfield(job.lpf_choice,'lpf_gauss_On')
    LPF.lpf_gauss_On = 1;
    LPF.fwhm1 = job.lpf_choice.lpf_gauss_On.fwhm1;
else
    LPF.lpf_gauss_On = 0;
end
%Spatial filter
if isfield(job.spatial_LPF,'spatial_LPF_On')
    radius = job.spatial_LPF.spatial_LPF_On.spatial_LPF_radius;
    spatial_LPF = 1;
else
    spatial_LPF = 0;
end
%shrinking of images
if isfield(job.shrinkage_choice,'configuration_shrink')
    shrink_x = job.shrinkage_choice.configuration_shrink.shrink_x;
    shrink_y = job.shrinkage_choice.configuration_shrink.shrink_y;
    try
        force_shrink_recompute = job.shrinkage_choice.configuration_shrink.force_shrink_recompute;
    catch
        force_shrink_recompute = 0;
    end
    shrinkage_choice = 1;
else
    shrinkage_choice = 0;
end

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
%Other options
include_OD = job.include_OD;
include_flow = job.include_flow;
include_HbT = job.include_HbT;
%save_figures
save_figures = job.save_figures;
generate_figures = job.generate_figures;
normalize_choice = job.normalize_choice;
try
    window_offset = job.window_offset;
catch
    window_offset = 0;
    job.window_offset = 0;
end
for SubjIdx=1:length(job.IOImat)
    try
        
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
        
        if ~isfield(IOI.res,'meanimageOK') || job.force_redo
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
            
            %get stimulation information - Careful, here onset duration is ignored!
            Ns = length(IOI.sess_res);
            if include_HbT
                if ~isfield(IOI.color,'HbT')
                    IOI.color.HbT = 'T';
                    IOI.color.eng = [IOI.color.eng IOI.color.HbT];
                end
            end
            Nc = length(IOI.color.eng);
            %save shrunk images
            if shrinkage_choice
                IOI = ioi_save_shrunk_images(IOI,job,SH,dir_ioimat);
            end
            
            %loop over sessions
            for s1=1:Ns
                if all_sessions || sum(s1==selected_sessions)
                    tmp_onsets = IOI.sess_res{s1}.onsets;
                    %remove onsets
                    if ~isempty(rmi)
                        IOI.rmi = rmi;
                        for k0=1:length(tmp_onsets)
                            r_onsets = round(tmp_onsets{k0});
                            k_onsets = [];
                            rem_onsets = [];
                            for u0=1:length(r_onsets)
                                if any(r_onsets(u0) == rmi) || any(r_onsets(u0) == rmi + 1) ||  any(r_onsets(u0) == rmi - 1)
                                    %skip this onset
                                    rem_onsets = [rem_onsets tmp_onsets{k0}(u0)];
                                else
                                    k_onsets = [k_onsets tmp_onsets{k0}(u0)];
                                end
                            end
                            IOI.onsets_kept{s1}{k0} = k_onsets;
                            IOI.onsets_removed{s1}{k0} = rem_onsets;
                            tmp_onsets{k0} = k_onsets;
                        end
                    else
                        IOI.onsets_kept{s1} = tmp_onsets;
                        IOI.onsets_removed{s1} = '';
                    end
                    if ~isempty(ust)
                        for k0=1:length(tmp_onsets)
                            try
                                tmp_onsets{k0} = tmp_onsets{k0}(ust);
                            catch
                                disp('Problem restricting number of onsets to user-specified -- probably not enough onsets');
                            end
                            try
                                if ~any(k0 == job.which_onset_type)
                                    tmp_onsets{k0} = ''; %remove these onset types
                                end
                            end
                        end
                    end
                    onsets_list{s1} = tmp_onsets;
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
            
            
            %Keep track later (for GLM) of which ROIs were selected
            Nc = length(IOI.color.eng);
            
            
            %retravailler
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
                                        if include_HbT
                                            if IOI.color.eng(c1) == IOI.color.HbT
                                                tmp_d = ROI{r1}{s1,IOI.color.eng == IOI.color.HbO}+...
                                                    ROI{r1}{s1,IOI.color.eng == IOI.color.HbR};
                                            else
                                                tmp_d = ROI{r1}{s1,c1};
                                            end
                                        else
                                            tmp_d = ROI{r1}{s1,c1};
                                        end
                                        %normalize flow
                                        if isfield(IOI.color,'flow')
                                            if IOI.color.eng(c1)==IOI.color.flow
                                                tmp_d = tmp_d/mean(tmp_d); %or median
                                            end
                                        end
                                        
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
            
            IOI.res.meanimageOK = 1;
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
        end
                
        disp(['Subject ' int2str(SubjIdx) ' complete']);
        out.IOImat{SubjIdx} = IOImat;
    catch exception
        disp(exception.identifier)
        disp(exception.stack(1))
        out.IOImat{SubjIdx} = job.IOImat{SubjIdx};
    end
end