function out = ioi_stim_mean_image_run(job)
%select a subset of sessions
[all_sessions selected_sessions] = ioi_get_sessions(job);
%filters
HPF = ioi_get_HPF(job);
LPF = ioi_get_LPF(job);

%Spatial filter
if isfield(job.spatial_LPF,'spatial_LPF_On')
    radius = job.spatial_LPF.spatial_LPF_On.spatial_LPF_radius;
    spatial_LPF = 1;
else
    spatial_LPF = 0;
end
%shrinking of images
[shrinkage_choice SH] = ioi_get_shrinkage_choice(job);

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
IC = job.IC; %colors to include

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
interactive_mode = job.interactive_mode;
do_another_stim = 0; %Boolean, only turned on in interactive mode
for SubjIdx=1:length(job.IOImat)
    try
        clear IOI ROI onsets_list M
        %Load IOI.mat information
        IOImat = job.IOImat{SubjIdx};
        load(IOImat);
        %check whether to skip any calculations, and create new IOI directory
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
            %create directory for figures
            if save_figures
                dir_fig = fullfile(newDir,'fig');
                if ~exist(dir_fig,'dir'),mkdir(dir_fig);end
            end
            if ~isfield(IOI,'dev')
                IOI.dev.TR = 0.2;
            end
            %rescale windows to data points (round off to integers)
            window_after = round(job.window_after/IOI.dev.TR);
            window_before = round(job.window_before/IOI.dev.TR);
            window_offset = round(window_offset/IOI.dev.TR);
            %Include HbT
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
            %restric onsets
            [IOI onsets_list] = ioi_restrict_onsets(IOI,job,rmi,ust);
            
            for s1=1:length(IOI.sess_res)
                if all_sessions || sum(s1==selected_sessions)
                    %loop over colors
                    for c1=1:Nc
                        %check whether this color should be done or not
                        doColor = ioi_doColor(IOI,c1,include_OD,include_flow,include_HbT,include_HbR,include_HbO);
                        if doColor
                            %load image data
                            y_first_pass = 1;
                            y = ioi_get_images(IOI,1:IOI.sess_res{s1}.n_frames,c1,s1,dir_ioimat,shrinkage_choice);
                            [nx ny nt] = size(y);
                            %loop over onset types
                            for m1=1:length(onsets_list{s1})
                                if isempty(job.which_onset_type) || any(m1==job.which_onset_type)
                                    onset_first_pass = 1;
                                    %do at least once, or as many times as requested by the user
                                    while onset_first_pass || interactive_mode || do_another_stim
                                        onset_first_pass = 0;
                                        if ~onset_first_pass && interactive_mode && do_another_stim
                                            %display interactive figure
                                        end
                                        
                                        if y_first_pass
                                            Sa = cell(Nroi,maxM);
                                            Sb = cell(Nroi,maxM);
                                            Ma = cell(Nroi,maxM);
                                            Mb = cell(Nroi,maxM);
                                            Da = cell(Nroi,maxM);
                                            Db = cell(Nroi,maxM);
                                            y_first_pass = 0;
                                        end
                                        Gtmp_array_before = zeros(1,window_before);
                                        Gtmp_array_after = zeros(1,window_after);
                                        GSb{r2,m1}{c1} = [];
                                        GSa{r2,m1}{c1} = [];
                                        Gkb = 0; %counter of segments before onsets
                                        Gka = 0; %counter of segments after onsets
                                        Gkb2 = 0; %counter of skipped segments before onsets
                                        Gka2 = 0; %counter of skipped segments after onsets
                                        %loop over sessions
                                        
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
                                                y = ioi_filter_HPF_LPF_WMDL(K,y')';
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
                                                        end
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
                                
                                
                                if (ka2>0 || kb2 > 0) && r2 == 1
                                    disp(['Skipped ' int2str(ka2) ' segments after onsets and ' int2str(kb2) ' segments before onsets']);
                                end
                            end
                        end
                    end
                end
            end
            
            IOI.res.meanimageOK = 1;
            
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