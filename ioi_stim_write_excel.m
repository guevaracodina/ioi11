function ioi_stim_write_excel(job,IOI,onsets_list,ROI,dir_fig)
%code based on ioi_stim_figures
[all_sessions selected_sessions] = ioi_get_sessions(job);
[all_ROIs selected_ROIs] = ioi_get_ROIs(job);
EBS = 3; %error bar step
add_error_bars = job.add_error_bars;
window_after = round(job.window_after/IOI.dev.TR);
Ma = IOI.res.Ma; %mean of "after segments", by region, stimulus type, and by color and region
%GMa = IOI.res.GMa;
Da = IOI.res.Da; %standard error
%Dma = IOI.res.Dma; %mean of standard error - single number
%Sa = IOI.res.Sa;
IC = job.IC; %colors to include

%Global Figures
ctotal = [];
h1 = 0;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
                        leg_str = {};
                        for c1 = ctotal 
                            if isfield(IOI.res.ROI{r1},'name')
                                tit = [IOI.subj_name ' ROI ' int2str(r1) ' ' IOI.res.ROI{r1}.name ', Session ' int2str(s1) ', Stimulus ' int2str(m1)];
                            else
                                tit = [IOI.subj_name ' ROI ' int2str(r1) ' ' int2str(r1) ', Session ' int2str(s1) ', Stimulus ' int2str(m1)];
                            end
                            v = Ma{r2,m1}{c1,s1}';
                            ioi_write_excel_core(job.write_excel_text,v,'mean',tit,dir_fig);
                            v = Da{r2,m1}{c1,s1}';
                            ioi_write_excel_core(job.write_excel_text,v,'std',tit,dir_fig);
                        end
                    end
                end
                
                %Generate a combined figure for all stims
                for c1 = ctotal
                    for r1=1:length(ROI)
                        if all_ROIs || sum(r1==selected_ROIs)
                            d = IOI.res.Sa{r1,m1}{c1,1};
                            if job.write_points_before
                                db = IOI.res.Sb{r1,m1}{c1,1};
                                d = [db d];
                            end
                            tit = [IOI.subj_name ', Session ' int2str(s1) ', Stim type ' int2str(m1) ', Color ' IOI.color.eng(c1) ', ROI ' int2str(r1)];
                            pathGF3 = fullfile(dir_fig,'StimGroup');
                            if ~exist(pathGF3,'dir'), mkdir(pathGF3); end
                            %Write with time in rows and stims in columns
                            ioi_write_excel_core(job.write_excel_text,d','allstims',tit,pathGF3);
                        end
                    end
                end
            end
        end        
    end
end