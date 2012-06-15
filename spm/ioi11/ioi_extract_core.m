function [ROI IOI] = ioi_extract_core(IOI,job,mask)
[all_sessions selected_sessions] = ioi_get_sessions(job);
[all_ROIs selected_ROIs] = ioi_get_ROIs(job);
IC = job.IC;
%loop over sessions
for s1=1:length(IOI.sess_res)
    if all_sessions || sum(s1==selected_sessions)
        %loop over available colors
        for c1=1:length(IOI.sess_res{s1}.fname)
            doColor = ioi_doColor(IOI,c1,IC);
            if doColor
                colorOK = 1;
                if ~(IOI.color.eng(c1)==IOI.color.laser)
                    %skip laser - only extract for flow
                    fname_list = IOI.sess_res{s1}.fname{c1};
                    %initialize
                    for r1=1:length(IOI.res.ROI)
                        if all_ROIs || sum(r1==selected_ROIs)
                            ROI{r1}{s1,c1} = [];
                        end
                    end
                    %loop over files
                    for f1=1:length(fname_list)
                        try
                            fname = fname_list{f1};
                            vols = spm_vol(fname);
                            d = spm_read_vols(vols);
                            [d1 d2 d3 d4] = size(d);
                            if d1 <= 1 || d2 <= 1
                                colorOK = 0;
                            end
                        catch
                            colorOK = 0;
                        end
                        %time dimension in 3rd dimension for colors
                        %R, G, Y, but in 4th dimension for O, D, F
                        %Loop over ROIs
                        [IOI ROI] = ioi_extract_main(IOI,ROI,job,d,d3,d4,c1,s1,colorOK,mask);
                    end
                    if colorOK
                        disp(['ROIs for session ' int2str(s1) ' and color ' IOI.color.eng(c1) ' completed']);
                    end
                end
            end
        end
    end
end