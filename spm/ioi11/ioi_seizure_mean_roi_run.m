
function out = ioi_seizure_mean_roi_run(job)
%_______________________________________________________________________
% Copyright (C) 2011 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
[all_sessions selected_sessions] = ioi_get_sessions(job);
[all_ROIs selected_ROIs] = ioi_get_ROIs(job);
%filters
HPF = ioi_get_HPF(job);
LPF = ioi_get_LPF(job);
IC = job.IC; %colors to include
save_figures = job.save_figures;
generate_figures = job.generate_figures;
normalize_choice = job.normalize_choice;
window_offset = job.window_offset;
add_error_bars = job.add_error_bars;
resampleOn=1;
HG=cell(1,4);
for i=1:4
    HG{i}=struct();
end
for SubjIdx=1:length(job.IOImat)
    try
        tic
        [IOI IOImat dir_ioimat]= ioi_get_IOI(job,SubjIdx);
        if ~isfield(IOI,'dev')
            IOI.dev.TR = 0.2;
        end
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
        window_after = round(job.window_after/IOI.dev.TR);
        window_before = round(job.window_before/IOI.dev.TR);
        window_offset = round(window_offset/IOI.dev.TR);
        if ~isfield(IOI.res,'meanOK') || job.force_redo
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
                dir_fig = fullfile(dir_ioimat,'fig');
                if ~exist(dir_fig,'dir'),mkdir(dir_fig);
                end
            end
            Ns = length(IOI.sess_res);
            if IC.include_HbT
                if ~isfield(IOI.color,'HbT')
                    IOI.color.HbT = 'T';
                    IOI.color.eng = [IOI.color.eng IOI.color.HbT];
                end
            end
            if all_ROIs
                selected_ROIs = 1:length(ROI);
            end
            %********************later onsets will be used
            %                 [IOI onsets_list pars_list] = ioi_restrict_onsets(IOI,job,rmi,ust);
            IOI.res.mean.selected_ROIs = selected_ROIs;
            keep_going = spm_input('Enter the number of seizures ',1,'y/n',[1 0]);
            if keep_going
                Nsz = spm_input('Enter the number of seizures ','+1');
                for i=1:Nsz  %****number of seizures
                    tb{i} = spm_input('Enter the start time of the seizure ','+1');
                    ta{i} = spm_input('Enter the end time of the seizure ','+1');
                end
            end
            HG{1}.name1='HbO';
            HG{2}.name1='HbR';
            HG{3}.name1='Flow';
            HG{4}.name1='HbT';
            HG{1}.name2='HbO_baseline';
            HG{2}.name2='HbR_baseline';
            HG{3}.name2='Flow_baseline';
            HG{4}.name2='HbT_baseline';
            HG{1}.name3='HbO_seizure';
            HG{2}.name3='HbR_seizure';
            HG{3}.name3='Flow_seizure';
            HG{4}.name3='HbT_seizure';
            HG{1}.name2_mean='HbO_baseline_mean';
            HG{2}.name2_mean='HbR_baseline_mean';
            HG{3}.name2_mean='Flow_baseline_mean';
            HG{4}.name2_mean='HbT_baseline_mean';
            for s1=1:length(IOI.sess_res)
                if all_sessions || sum(s1==selected_sessions)
                    r2=0;
                    for r1=1:length(ROI)
                        if all_ROIs || sum(r1==selected_ROIs)
                            h1=0; ctotal = [];                          
                                for i=5:7
                                    if LPF.lpf_gauss_On
                                        K = get_K(1:length(ROI{1,r1}{s1,i}),LPF.fwhm1,IOI.dev.TR);
                                        ROI{1,r1}{s1,i} = ioi_filter_HPF_LPF_WMDL(K,(ROI{1,r1}{s1,i})')';
                                    end
                                end
                                for i=1:Nsz
                                    index1{i}=((tb{i}/IOI.dev.TR)-window_before-window_offset):((tb{i}/IOI.dev.TR)-window_offset); % the suration of baseline
                                    index{i}=((tb{i}/IOI.dev.TR)-window_offset):((ta{i}/IOI.dev.TR)+window_after);  % the duration of each seizure
                                    sz_duration(i)=length(index{i});
                                    if IC.include_HbO
                                        HG{1}.value1=ROI{1,r1}{s1,5};
                                        HG{1}.value1=HG{1}.value1-mean(HG{1}.value1)+baseline_hbo; %normalization
                                        %                                             if show_percent_changes
                                        HG{1}.value1=100*(HG{1}.value1/baseline_hbo-1);
                                        %                                             end
                                        HG{1}.value2(i,:)= HG{1}.value1(index1{i});    %%%%%%%%%%%%%calculate the baseline
                                        HG{1}.value3{i}= HG{1}.value1(index{i});     %%%%%%%%%%%%%%%%%%%calculate the seizure
                                    end
                                    if IC.include_HbR
                                        HG{2}.value1=ROI{1,r1}{s1,6};
                                        HG{2}.value1=HG{2}.value1-mean(HG{2}.value1)+baseline_hbr;
                                        %                                             if show_percent_changes
                                        HG{2}.value1=100*(HG{2}.value1/baseline_hbr-1);
                                        %                                             end
                                        HG{2}.value2(i,:)= HG{2}.value1(index1{i});
                                        HG{2}.value3{i}= HG{2}.value1(index{i});
                                    end
                                    if IC.include_flow
                                        HG{3}.value1=ROI{1,r1}{s1,7};
                                        HG{3}.value2(i,:)= HG{3}.value1(index1{i});
                                        HG{3}.value3{i}= HG{3}.value1(index{i});
                                    end
                                    if IC.include_HbT
                                        HG{4}.value1=(ROI{1,r1}{s1,5}+ROI{1,r1}{s1,6})-mean((ROI{1,r1}{s1,5}+ROI{1,r1}{s1,6}))+baseline_hbt;
                                        %                                             if show_percent_changes
                                        HG{4}.value1=100*(HG{4}.value1/baseline_hbt-1);
                                        %                                             end
                                        HG{4}.value2(i,:)= HG{4}.value1(index1{i});
                                        HG{4}.value3{i}= HG{4}.value1(index{i});
                                    end
                                end
                                 duration_mean=fix(mean(sz_duration));
                                %************calculate the mean of
                                %baseline
                                if resampleOn
                                    for i=1:Nsz
                                        
                                        if IC.include_HbO
                                            HG{1}.value2_mean=mean(HG{1}.value2,2);
                                            HG{1}.value3{i}= resample(HG{1}.value3{i},duration_mean,length(HG{1}.value3{i}));
                                            HG{1}.value4(i,:)=HG{1}.value3{i}-HG{1}.value2_mean(i);
                                            HG{1}.MeanSz=mean(HG{1}.value4);
                                            HG{1}.StdSz=std(HG{1}.value4,0,1)/sqrt(Nsz);
                                            IOI.sess_res{s1}.meanSeizure{1}=HG{1}.MeanSz;
                                            IOI.sess_res{s1}.StdSeizure{1}=HG{1}.StdSz;
                                        end
                                        if IC.include_HbR
                                            HG{2}.value2_mean=mean(HG{2}.value2,2);
                                            HG{2}.value3{i}= resample(HG{2}.value3{i},duration_mean,length(HG{2}.value3{i}));
                                            HG{2}.value4(i,:)=HG{2}.value3{i}-HG{2}.value2_mean(i);
                                            HG{2}.MeanSz=mean(HG{2}.value4);
                                            HG{2}.StdSz=std(HG{2}.value4,0,1)/sqrt(Nsz);
                                            IOI.sess_res{s1}.meanSeizure{1}=HG{2}.MeanSz;
                                            IOI.sess_res{s1}.StdSeizure{1}=HG{2}.StdSz;
                                        end
                                        if IC.include_flow
                                            HG{3}.value2_mean=mean(HG{3}.value2,2);
                                            HG{3}.value3{i}= resample(HG{3}.value3{i},duration_mean,length(HG{3}.value3{i}));
                                            HG{3}.value4(i,:)=HG{3}.value3{i}-HG{3}.value2_mean(i);
                                            HG{3}.MeanSz=mean(HG{3}.value4);
                                            HG{3}.StdSz=std(HG{3}.value4,0,1)/sqrt(Nsz);
                                            IOI.sess_res{s1}.meanSeizure{1}=HG{3}.MeanSz;
                                            IOI.sess_res{s1}.StdSeizure{1}=HG{3}.StdSz;
                                        end
                                        if IC.include_HbT
                                            HG{4}.value2_mean=mean(HG{4}.value2,2);
                                            HG{4}.value3{i}= resample(HG{4}.value3{i},duration_mean,length(HG{4}.value3{i}));
                                            HG{4}.value4(i,:)=HG{4}.value3{i}-HG{4}.value2_mean(i);
                                            HG{4}.MeanSz=mean(HG{4}.value4);
                                            HG{4}.StdSz=std(HG{4}.value4,0,1)/sqrt(Nsz);
                                            IOI.sess_res{s1}.meanSeizure{1}=HG{4}.MeanSz;
                                            IOI.sess_res{s1}.StdSeizure{1}=HG{4}.StdSz;
                                        end
                                    end
                                end
                                save(IOImat,'IOI');
                                                     
                           

                            %**************caculate the mean and errorbar
%                             for j=1:4
%                                 HG{j}.MeanSz=mean(HG{j}.value4);
%                                 HG{j}.StdSz=std(HG{j}.value4,0,1)/sqrt(Nsz);
%                                 IOI.sess_res{s1}.meanSeizure{j}=HG{j}.MeanSz;
%                                 IOI.sess_res{s1}.StdSeizure{j}=HG{j}.StdSz;
%                             end
%                             IOI.res.meanSZ=HG;
%                             IOI.res.StdSz=HG.StdSz;
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
                               % lp=linspace(-time_offset,length(HG{1}.MeanSz)/fs-time_offset,length(HG{1}.MeanSz));
                                ls=linspace(-job.window_offset,(length(HG{1}.MeanSz)*IOI.dev.TR)-job.window_offset,length(HG{1}.MeanSz));
                                r2 = r2+1;
                                h1 = h1 + 1;
                                h(h1) = figure;
                                leg_str = {};
                                for c1 = ctotal
                                    try %in case we do not have this color
                                        if ~add_error_bars
                                            plot(ls,HG{c1-4}.MeanSz,[lp1{1} lp2{c1}]); hold on
                                        else
                                            errorbar(ls(1:end),HG{c1-4}.MeanSz,HG{c1-4}.StdSz,[lp1{1} lp2{c1}]); hold on
                                        end
                                        leg_str = [leg_str; lp3{c1}];
                                    end
                                end
                                legend(leg_str);
                                if save_figures
                                    if isfield(IOI.res.ROI{r1},'name')
                                        tit = [IOI.subj_name ' ROI ' int2str(r1) ' ' IOI.res.ROI{r1}.name ', Session ' int2str(s1)];
                                    else
                                        tit = [IOI.subj_name ' ROI ' int2str(r1) ' ' int2str(r1) ', Session ' int2str(s1)];
                                    end
                                    title(tit);
                                    ioi_save_figures(save_figures,generate_figures,h(h1),tit,dir_fig);
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

function out = ioi_seizure_mean_roi_run(job)
%_______________________________________________________________________
% Copyright (C) 2011 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
[all_sessions selected_sessions] = ioi_get_sessions(job);
[all_ROIs selected_ROIs] = ioi_get_ROIs(job);
%filters
% HPF = ioi_get_HPF(job);
LPF = ioi_get_LPF(job);
IC = job.IC; %colors to include
save_figures = job.save_figures;
generate_figures = job.generate_figures;
% normalize_choice = job.normalize_choice;
window_offset = job.window_offset;
add_error_bars = job.add_error_bars;
resampleOn=1;
HG=cell(1,4);
for i=1:4
    HG{i}=struct();
end
for SubjIdx=1:length(job.IOImat)
    try
        tic
        [IOI IOImat dir_ioimat]= ioi_get_IOI(job,SubjIdx);
        if ~isfield(IOI,'dev')
            IOI.dev.TR = 0.2;
        end
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
        window_after = round(job.window_after/IOI.dev.TR);
        window_before = round(job.window_before/IOI.dev.TR);
        window_offset = round(window_offset/IOI.dev.TR);
        if ~isfield(IOI.res,'meanOK') || job.force_redo
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
                dir_fig = fullfile(dir_ioimat,'fig');
                if ~exist(dir_fig,'dir'),mkdir(dir_fig);
                end
            end
            Ns = length(IOI.sess_res);
            if IC.include_HbT
                if ~isfield(IOI.color,'HbT')
                    IOI.color.HbT = 'T';
                    IOI.color.eng = [IOI.color.eng IOI.color.HbT];
                end
            end
            if all_ROIs
                selected_ROIs = 1:length(ROI);
            end
            %********************later onsets will be used
            %                 [IOI onsets_list pars_list] = ioi_restrict_onsets(IOI,job,rmi,ust);
            IOI.res.mean.selected_ROIs = selected_ROIs;
            keep_going = spm_input('Enter the number of seizures ',1,'y/n',[1 0]);
            if keep_going
                Nsz = spm_input('Enter the number of seizures ','+1');
                for i=1:Nsz  %****number of seizures
                    tb{i} = spm_input('Enter the start time of the seizure ','+1');
                    ta{i} = spm_input('Enter the end time of the seizure ','+1');
                end
            end
            HG{1}.name1='HbO';
            HG{2}.name1='HbR';
            HG{3}.name1='Flow';
            HG{4}.name1='HbT';
            HG{1}.name2='HbO_baseline';
            HG{2}.name2='HbR_baseline';
            HG{3}.name2='Flow_baseline';
            HG{4}.name2='HbT_baseline';
            HG{1}.name3='HbO_seizure';
            HG{2}.name3='HbR_seizure';
            HG{3}.name3='Flow_seizure';
            HG{4}.name3='HbT_seizure';
            HG{1}.name2_mean='HbO_baseline_mean';
            HG{2}.name2_mean='HbR_baseline_mean';
            HG{3}.name2_mean='Flow_baseline_mean';
            HG{4}.name2_mean='HbT_baseline_mean';
            for s1=1:length(IOI.sess_res)
                if all_sessions || sum(s1==selected_sessions)
                    r2=0;
                    for r1=1:length(ROI)
                        if all_ROIs || sum(r1==selected_ROIs)
                            h1=0; ctotal = [];                          
                                for i=5:7
                                    if LPF.lpf_gauss_On
                                        K = get_K(1:length(ROI{1,r1}{s1,i}),LPF.fwhm1,IOI.dev.TR);
                                        ROI{1,r1}{s1,i} = ioi_filter_HPF_LPF_WMDL(K,(ROI{1,r1}{s1,i})')';
                                    end
                                end
                                for i=1:Nsz
                                    index1{i}=((tb{i}/IOI.dev.TR)-window_before-window_offset):((tb{i}/IOI.dev.TR)-window_offset); % the suration of baseline
                                    index{i}=((tb{i}/IOI.dev.TR)-window_offset):((ta{i}/IOI.dev.TR)+window_after);  % the duration of each seizure
                                    sz_duration(i)=length(index{i});
                                    if IC.include_HbO
                                        HG{1}.value1=ROI{1,r1}{s1,5};
                                        HG{1}.value1=HG{1}.value1-mean(HG{1}.value1)+baseline_hbo; %normalization
                                        %                                             if show_percent_changes
                                        HG{1}.value1=100*(HG{1}.value1/baseline_hbo-1);
                                        %                                             end
                                        HG{1}.value2(i,:)= HG{1}.value1(index1{i});    %%%%%%%%%%%%%calculate the baseline
                                        HG{1}.value3{i}= HG{1}.value1(index{i});     %%%%%%%%%%%%%%%%%%%calculate the seizure
                                    end
                                    if IC.include_HbR
                                        HG{2}.value1=ROI{1,r1}{s1,6};
                                        HG{2}.value1=HG{2}.value1-mean(HG{2}.value1)+baseline_hbr;
                                        %                                             if show_percent_changes
                                        HG{2}.value1=100*(HG{2}.value1/baseline_hbr-1);
                                        %                                             end
                                        HG{2}.value2(i,:)= HG{2}.value1(index1{i});
                                        HG{2}.value3{i}= HG{2}.value1(index{i});
                                    end
                                    if IC.include_flow
                                        HG{3}.value1=ROI{1,r1}{s1,7};
                                        HG{3}.value2(i,:)= HG{3}.value1(index1{i});
                                        HG{3}.value3{i}= HG{3}.value1(index{i});
                                    end
                                    if IC.include_HbT
                                        HG{4}.value1=(ROI{1,r1}{s1,5}+ROI{1,r1}{s1,6})-mean((ROI{1,r1}{s1,5}+ROI{1,r1}{s1,6}))+baseline_hbt;
                                        %                                             if show_percent_changes
                                        HG{4}.value1=100*(HG{4}.value1/baseline_hbt-1);
                                        %                                             end
                                        HG{4}.value2(i,:)= HG{4}.value1(index1{i});
                                        HG{4}.value3{i}= HG{4}.value1(index{i});
                                    end
                                end
                                 duration_mean=fix(mean(sz_duration));
                                %************calculate the mean of
                                %baseline
                                if resampleOn
                                    for i=1:Nsz
                                        
                                        if IC.include_HbO
                                            HG{1}.value2_mean=mean(HG{1}.value2,2);
                                            HG{1}.value3{i}= resample(HG{1}.value3{i},duration_mean,length(HG{1}.value3{i}));
                                            HG{1}.value4(i,:)=HG{1}.value3{i}-HG{1}.value2_mean(i);
                                            HG{1}.MeanSz=mean(HG{1}.value4);
                                            HG{1}.StdSz=std(HG{1}.value4,0,1)/sqrt(Nsz);
                                            IOI.sess_res{s1}.meanSeizure{1}=HG{1}.MeanSz;
                                            IOI.sess_res{s1}.StdSeizure{1}=HG{1}.StdSz;
                                        end
                                        if IC.include_HbR
                                            HG{2}.value2_mean=mean(HG{2}.value2,2);
                                            HG{2}.value3{i}= resample(HG{2}.value3{i},duration_mean,length(HG{2}.value3{i}));
                                            HG{2}.value4(i,:)=HG{2}.value3{i}-HG{2}.value2_mean(i);
                                            HG{2}.MeanSz=mean(HG{2}.value4);
                                            HG{2}.StdSz=std(HG{2}.value4,0,1)/sqrt(Nsz);
                                            IOI.sess_res{s1}.meanSeizure{1}=HG{2}.MeanSz;
                                            IOI.sess_res{s1}.StdSeizure{1}=HG{2}.StdSz;
                                        end
                                        if IC.include_flow
                                            HG{3}.value2_mean=mean(HG{3}.value2,2);
                                            HG{3}.value3{i}= resample(HG{3}.value3{i},duration_mean,length(HG{3}.value3{i}));
                                            HG{3}.value4(i,:)=HG{3}.value3{i}-HG{3}.value2_mean(i);
                                            HG{3}.MeanSz=mean(HG{3}.value4);
                                            HG{3}.StdSz=std(HG{3}.value4,0,1)/sqrt(Nsz);
                                            IOI.sess_res{s1}.meanSeizure{1}=HG{3}.MeanSz;
                                            IOI.sess_res{s1}.StdSeizure{1}=HG{3}.StdSz;
                                        end
                                        if IC.include_HbT
                                            HG{4}.value2_mean=mean(HG{4}.value2,2);
                                            HG{4}.value3{i}= resample(HG{4}.value3{i},duration_mean,length(HG{4}.value3{i}));
                                            HG{4}.value4(i,:)=HG{4}.value3{i}-HG{4}.value2_mean(i);
                                            HG{4}.MeanSz=mean(HG{4}.value4);
                                            HG{4}.StdSz=std(HG{4}.value4,0,1)/sqrt(Nsz);
                                            IOI.sess_res{s1}.meanSeizure{1}=HG{4}.MeanSz;
                                            IOI.sess_res{s1}.StdSeizure{1}=HG{4}.StdSz;
                                        end
                                    end
                                end
                                save(IOImat,'IOI');
                                                     
                           

                            %**************caculate the mean and errorbar
%                             for j=1:4
%                                 HG{j}.MeanSz=mean(HG{j}.value4);
%                                 HG{j}.StdSz=std(HG{j}.value4,0,1)/sqrt(Nsz);
%                                 IOI.sess_res{s1}.meanSeizure{j}=HG{j}.MeanSz;
%                                 IOI.sess_res{s1}.StdSeizure{j}=HG{j}.StdSz;
%                             end
%                             IOI.res.meanSZ=HG;
%                             IOI.res.StdSz=HG.StdSz;
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
                               % lp=linspace(-time_offset,length(HG{1}.MeanSz)/fs-time_offset,length(HG{1}.MeanSz));
                                ls=linspace(-job.window_offset,(length(HG{1}.MeanSz)*IOI.dev.TR)-job.window_offset,length(HG{1}.MeanSz));
                                r2 = r2+1;
                                h1 = h1 + 1;
                                h(h1) = figure;
                                leg_str = {};
                                for c1 = ctotal
                                    try %in case we do not have this color
                                        if ~add_error_bars
                                            plot(ls,HG{c1-4}.MeanSz,[lp1{1} lp2{c1}]); hold on
                                        else
                                            errorbar(ls(1:end),HG{c1-4}.MeanSz,HG{c1-4}.StdSz,[lp1{1} lp2{c1}]); hold on
                                        end
                                        leg_str = [leg_str; lp3{c1}];
                                    end
                                end
                                legend(leg_str);
                                if save_figures
                                    if isfield(IOI.res.ROI{r1},'name')
                                        tit = [IOI.subj_name ' ROI ' int2str(r1) ' ' IOI.res.ROI{r1}.name ', Session ' int2str(s1)];
                                    else
                                        tit = [IOI.subj_name ' ROI ' int2str(r1) ' ' int2str(r1) ', Session ' int2str(s1)];
                                    end
                                    title(tit);
                                    ioi_save_figures(save_figures,generate_figures,h(h1),tit,dir_fig);
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

