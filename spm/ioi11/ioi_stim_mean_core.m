function [Rs Sb Sa Ma Da Dma U0 GK] = ioi_stim_mean_core(job,IOI,d,onsets_list,pars_list,...
    Rs,Sb,Sa,Ma,Da,Dma,GSb,GSa,global_M,PGM,r2,r1,m1,c1,s1,GK,norm1,norm2)
if isfield(job,'remove_stims_SD')
    remove_stims_SD = job.remove_stims_SD;
else
    remove_stims_SD = 1;
end
U0 = [];
Weigh_amplitude_of_threshold = job.Weigh_amplitude_of_threshold;
remove_segment_drift = job.remove_segment_drift;
normalize_choice = job.normalize_choice;

window_after = round(job.window_after/IOI.dev.TR);
window_before = round(job.window_before/IOI.dev.TR);

pass_twice = 1;
removeSeg = [];
for pass_cnt = 1:2
    %Two passes: the first pass is used to calculate the standard deviation, as
    %a criterion to remove stims; the average is then recalculated in the second pass
    if pass_cnt == 1 || (pass_cnt == 2 && pass_twice)
        kb = 0; %counter of segments before onsets
        ka = 0; %counter of segments after onsets
        kb2 = 0; %counter of skipped segments before onsets
        ka2 = 0; %counter of skipped segments after onsets
        Sb{r2,m1}{c1,s1} = []; %list of segments before
        Sa{r2,m1}{c1,s1} = []; %list of segments after
        %loop over onsets for that session
        if ~isempty(onsets_list{s1}{m1})
            U = round((onsets_list{s1}{m1}-job.window_offset)/IOI.dev.TR); %in data points
            U0{s1} = ioi_get_U(IOI,[],U,0,s1); %only for plotting stims
            
            for u1=1:length(U)
                if ~any(u1==removeSeg)
                    clear tmp_median;
                    try
                        tmp_median = d(U(u1)); %in case the following fails
                        tmp1 = d(U(u1)-window_before:U(u1)-1);
                        switch normalize_choice
                            case 1
                                tmp_median = median(tmp1);
                            case 2
                                tmp_median = tmp1(end);
                            case 3
                                tmp_median = mean(tmp1);
                        end
                        kb = kb+1;
                        GK.Gkb = GK.Gkb+1;
                        Sb{r2,m1}{c1,s1} = [Sb{r2,m1}{c1,s1};tmp1-tmp_median];
                        if global_M && m1 <= PGM
                            GSb{r2,m1}{c1} = [GSb{r2,m1}{c1};tmp1-tmp_median];
                        end
                    catch
                        kb2 = kb2+1;
                        GK.Gkb2 = GK.Gkb2+1;
                        if kb2 < 3 && r2 == 1
                            disp(['Could not include segment before onset ' int2str(u1) ' at time ' num2str(U(u1)*IOI.dev.TR) ...
                                ' for session ' int2str(s1) ' for ROI ' int2str(r1) ...
                                ' for color ' int2str(c1) ' for onset type ' int2str(m1) ...
                                ' in global average over all sessions... skipping ' int2str(kb2) ' so far']);
                        end
                    end
                    try
                        tmp1 = d(U(u1):U(u1)+window_after-1);
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
                        
                        if Weigh_amplitude_of_threshold
                            tmp1 = tmp1/pars_list{s1}{m1}(u1);
                        end
                        
                        %normalize to get percent change
                        if job.mult_normalize_choice
                            tmp1 = tmp1/norm2;
                        else
                            tmp1 = tmp1/(norm1+tmp_median);
                        end
                        ka = ka+1;
                        GK.Gka = GK.Gka+1;
                        Sa{r2,m1}{c1,s1} = [Sa{r2,m1}{c1,s1};tmp1];
                        if global_M && m1 <= PGM
                            GSa{r2,m1}{c1} = [GSa{r2,m1}{c1};tmp1];
                        end
                    catch
                        ka2 = ka2+1;
                        GK.Gka2 = GK.Gka2+1;
                        if ka2<3 && r2 == 1
                            disp(['Could not include segment after onset ' int2str(u1) ' at time ' num2str(U(u1)*IOI.dev.TR) ...
                                ' for session ' int2str(s1) ' for ROI ' int2str(r1) ...
                                ' for color ' int2str(c1) ' for onset type ' int2str(m1) ...
                                ' in global average over all sessions... skipping ' int2str(ka2) ' so far']);
                        end
                    end
                end
            end
            
            Ma{r2,m1}{c1,s1} = mean(Sa{r2,m1}{c1,s1},1); %global mean after
            Da{r2,m1}{c1,s1} = std(Sa{r2,m1}{c1,s1},0,1)/sqrt(ka); %SEM
            Dma{r2,m1}{c1,s1} = mean(Da{r2,m1}{c1,s1}); %mean SEM
            Rs{r2,m1}{c1,s1} = removeSeg; %Segments that were removed based on standard deviation criterion
            
            if remove_stims_SD
                if pass_cnt == 1
                    meanA = mean(Ma{r2,m1}{c1,s1});
                    tmpSeg = Sa{r2,m1}{c1,s1};
                    meanSeg = mean(tmpSeg,2);
                    meanSd = std(tmpSeg(:));
                    for a0=1:length(meanSeg)
                        if abs(meanSeg(a0)-meanA) > job.std_choice*meanSd %very strong criterion perhaps better to keep it at 2*meanSd
                            removeSeg = [removeSeg a0];
                        end
                    end
                end
            else
                pass_twice = 0;
            end
        else
            pass_twice = 0;
        end
        GK.ka = ka;
        GK.kb = kb;
        GK.ka2 = ka2;
        GK.kb2 = kb2;
    end
end