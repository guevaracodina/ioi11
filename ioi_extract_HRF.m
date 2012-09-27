function IOI = ioi_extract_HRF(job,IOI,ROI,maxM,global_M,Ma,GMa)
[all_sessions selected_sessions] = ioi_get_sessions(job);
[all_ROIs selected_ROIs] = ioi_get_ROIs(job);
x = linspace(0,job.window_after-IOI.dev.TR,window_after);
fit_3_gamma = job.fit_3_gamma;
IC = job.IC; %colors to include
include_nlinfit = job.include_nlinfit;

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
                                if include_nlinfit
                                    F{r2,m1}{c1,s1} = ioi_nlinfit(x,d,IOI.color,c1,IC.include_flow);
                                end
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
                            if include_nlinfit
                                GF{r2,m1}{c1} = ioi_nlinfit(x,d,IOI.color,c1,IC.include_flow);
                            end
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