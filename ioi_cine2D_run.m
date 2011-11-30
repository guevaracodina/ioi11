function out = ioi_cine2D_run(job)
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
%HPF
if isfield(job.hpf_butter,'hpf_butter_On')
    HPF.hpf_butter_On = 1;
    HPF.hpf_butter_freq = job.hpf_butter.hpf_butter_On.hpf_butter_freq;
    HPF.hpf_butter_order = job.hpf_butter.hpf_butter_On.hpf_butter_order;
else
    HPF.hpf_butter_On = 0;
end
include_flow = job.include_flow;
normalize_choice = job.normalize_choice;
%get size of windows before and after stimulation to average on, in data points
for SubjIdx=1:length(job.IOImat)
    try
        tic
        clear IOI onsets_list M
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
            if ~isfield(IOI.res,'cineOK') || job.force_redo
                [dir_ioimat dummy] = fileparts(job.IOImat{SubjIdx});
                if isfield(job.IOImatCopyChoice,'IOImatCopy')
                    newDir = job.IOImatCopyChoice.IOImatCopy.NewIOIdir;
                    newDir = fullfile(dir_ioimat,newDir);
                    if ~exist(newDir,'dir'),mkdir(newDir); end
                    IOImat = fullfile(newDir,'IOI.mat');
                else
                    newDir = dir_ioimat;
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
                Nc = length(IOI.color.eng);
                %loop over sessions
                for s1=1:length(IOI.sess_res)
                    if all_sessions || sum(s1==selected_sessions)
                        
                        %loop over colors
                        for c1=5 %1:length(IOI.color.eng)
                            %loop over onset types
                            for m1=1:length(onsets_list{s1})
                                %Image arrays 
                                
                                %loop over onsets for that session
                                U = round(onsets_list{s1}{m1}/IOI.dev.TR); %in data points
                                for u1=1:length(U)
                                    %Frames to get:
                                    fr_before = U(u1)-window_before:U(u1)-1;
                                    fr_after = U(u1):U(u1)+window_after-1;
                                    %Load images before and after stimulus
                                    Yb = get_images(IOI,fr_before,c1,s1,dir_ioimat);
                                    Ya = get_images(IOI,fr_after,c1,s1,dir_ioimat);
                                    
                                        
                                        
                                        
                                
                                kb = 0; %counter of segments before onsets
                                ka = 0; %counter of segments after onsets
                                kb2 = 0; %counter of skipped segments before onsets
                                ka2 = 0; %counter of skipped segments after onsets
                                
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
                                    end
                                    Mb{r2,m1}{c1,s1} = tmp_array_before/kb; %global mean before
                                end
                            end
                        end
                        
                        if (ka2>0 || kb2 > 0) && r2 == 1
                            disp(['Skipped ' int2str(ka2) ' segments after onsets and ' int2str(kb2) ' segments before onsets']);
                        end
                    end
                end
                
            end
            IOI.res.cineOK = 1;
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

function Ya = get_images(IOI,frames,c1,s1,dir_ioimat)
Ya = [];
first_pass = 1;
fmod_previous = 1;
for i=1:length(frames)
    frames_loaded = 0;
    if i==1
        %find number of frames per block
        if length(IOI.sess_res{s1}.si) > 1
            fpb = IOI.sess_res{s1}.si{2}-IOI.sess_res{s1}.si{1};
        else
            fpb = IOI.sess_res{s1}.ei{1}-IOI.sess_res{s1}.si{1}+1;
        end
    end
    fmod = mod(frames(i),fpb);
    if fmod < fmod_previous 
        frames_loaded = 0;
    end
    fct = ceil(frames(i)/fpb);
    if ~frames_loaded
    try        
        V = spm_vol(IOI.sess_res{s1}.fname{c1}{fct});
        Y = spm_read_vols(V);
        frames_loaded = 1;
    catch
        [dir0 fil0 ext0] = fileparts(IOI.sess_res{s1}.fname{c1}{fct});
        fsep = strfind(dir0,filesep);
        res = strfind(dir0,'Res');
        fgr = fsep(fsep > res);
        tdir = dir0(1:fgr(2)); 
        tfname = fullfile(dir_ioimat,['S' gen_num_str(s1,2)],[fil0 ext0]);
        try
            V = spm_vol(tfname);
            Y = spm_read_vols(V);
            frames_loaded = 1;            
        catch
            %file does not exist
            return
        end
    end
    if first_pass
        %initialize Ya
        Ya = zeros(size(Y,1),size(Y,2),length(frames));
        first_pass = 0;
    end
    end
    
    Ya(:,:,i) = squeeze(Y(:,:,1,fmod));
end
end