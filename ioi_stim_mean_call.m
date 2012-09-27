function [IOI U0] = ioi_stim_mean_call(job,IOI,ROI,maxM,global_M,...
    PGM,onsets_list,pars_list)
[all_sessions selected_sessions] = ioi_get_sessions(job);
[all_ROIs selected_ROIs] = ioi_get_ROIs(job);
%Preallocate the arrays
if all_ROIs
    selected_ROIs = 1:length(ROI);
end
%Keep track later (for GLM) of which ROIs were selected
IOI.res.mean.selected_ROIs = selected_ROIs;
Nroi = length(selected_ROIs);
IC = job.IC; %colors to include
%filters
HPF = ioi_get_HPF(job);
LPF = ioi_get_LPF(job);
if ~isfield(IOI,'conc')
    baseline_hbt = 100;
    baseline_hbo = 60;
    baseline_hbr = 40;
    IOI.conc.baseline_hbt = baseline_hbt;
    IOI.conc.baseline_hbo = baseline_hbo;
    IOI.conc.baseline_hbr = baseline_hbr;
else
    baseline_hbt = IOI.conc.baseline_hbt;
    baseline_hbo = IOI.conc.baseline_hbo;
    baseline_hbr = IOI.conc.baseline_hbr;
end

if global_M
    GSa = cell(Nroi,PGM);
    GSb = cell(Nroi,PGM);
    GMa = cell(Nroi,PGM);
    GDa = cell(Nroi,PGM);
    GDma = cell(Nroi,PGM);
    %for c1=1:Nc
else
    GMa = [];
end
Sa = cell(Nroi,maxM);
Sb = cell(Nroi,maxM);
Ma = cell(Nroi,maxM);
Da = cell(Nroi,maxM);
Dma = cell(Nroi,maxM);
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
                %Gtmp_array_before = zeros(1,window_before);
                %Gtmp_array_after = zeros(1,window_after);
                GSb{r2,m1}{c1} = [];
                GSa{r2,m1}{c1} = [];
                GK.Gkb = 0; %counter of segments before onsets
                GK.Gka = 0; %counter of segments after onsets
                GK.Gkb2 = 0; %counter of skipped segments before onsets
                GK.Gka2 = 0; %counter of skipped segments after onsets
                %loop over sessions
                for s1=1:length(IOI.sess_res)
                    if all_sessions || sum(s1==selected_sessions)
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
                            switch IOI.color.eng(c1)
                                case 'T'
                                    norm1 = baseline_hbt;
                                    norm2 = norm1;
                                case 'O'
                                    norm1 = baseline_hbo;
                                    norm2 = norm1;
                                case 'D'
                                    norm1 = baseline_hbr;
                                    norm2 = norm1;
                                otherwise
                                    norm1 = 0;
                                    norm2 = mean(tmp_d);
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
                                tmp_DC = median(tmp_d);
                                tmp_d = ButterHPF(1/IOI.dev.TR,HPF.hpf_butter_freq,HPF.hpf_butter_order,tmp_d);
                                tmp_d = tmp_d + tmp_DC; %add back the DC component
                            end
                            if LPF.lpf_gauss_On
                                if ~LPF.apply_lpf_on_flow_only || (LPF.apply_lpf_on_flow_only && IOI.color.eng(c1) == IOI.color.flow)
                                    K = get_K(1:length(tmp_d),LPF.fwhm1,IOI.dev.TR);
                                    y = tmp_d;
                                    y = ioi_filter_HPF_LPF_WMDL(K,y')';
                                    tmp_d = y;
                                end
                            end
                            %averaging
                            [Rs Sb Sa Ma Da Dma U0 GK] = ioi_stim_mean_core(job,IOI,tmp_d,onsets_list,pars_list,...
                                Rs,Sb,Sa,Ma,Da,Dma,GSb,GSa,global_M,PGM,r2,r1,m1,c1,s1,GK,norm1,norm2);
                        end
                    end
                end
                
                if global_M && m1 <= PGM
                    GMa{r2,m1}{c1} = mean(GSa{r2,m1}{c1},1); %tmp_array_after/Gka; %global mean after
                    GDa{r2,m1}{c1} = std(GSa{r2,m1}{c1},0,1)/sqrt(GK.Gka); %SEM
                    GDma{r2,m1}{c1} = mean(GDa{r2,m1}{c1});
                end
                if (GK.ka2>0 || GK.kb2 > 0) && r2 == 1
                    disp(['Skipped ' int2str(GK.ka2) ' segments after onsets and ' int2str(GK.kb2) ' segments before onsets']);
                end
            end
        end
    end
end
%Global results
if global_M
    try
        IOI.res.GMa = GMa; %mean of all segments
        IOI.res.GDa = GDa; %standard error
        IOI.res.GDma = GDma; %mean of standard error - single number
        IOI.res.GSa = GSa; %each segment
        IOI.res.GSb = GSb;
    end
end
%Session results
IOI.res.Ma = Ma; %mean of "after segments", by region, stimulus type, and by color and region
IOI.res.Da = Da; %standard error
IOI.res.Dma = Dma; %mean of standard error - single number
IOI.res.Sa = Sa; %each segment
IOI.res.Sb = Sb;
IOI.res.Rs = Rs;