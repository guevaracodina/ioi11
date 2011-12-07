function out = ioi_explore_stim_run(job)
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
% %select nature of stimulation data to be used for averaging
% if isfield(job.ROI_choice,'select_ROIs')
%     all_ROIs = 0;
%     selected_ROIs = job.ROI_choice.select_ROIs.selected_ROIs;
% else
%     all_ROIs = 1;
% end
%save_figures
save_figures = job.save_figures;
normalize_choice = job.normalize_choice;
%get size of windows before and after stimulation to average on, in data points

generate_figures = job.generate_figures;
add_error_bars = job.add_error_bars;
which_ons = 1; %which onset type to use
for SubjIdx=1:length(job.IOImat)
    try
        tic
        clear IOI ROI onsets_list M
        %Load IOI.mat information
        IOImat = job.IOImat{SubjIdx};
        load(IOImat);
        window_after = round(job.window_after/IOI.dev.TR);
        window_before = round(job.window_before/IOI.dev.TR);
        
        %         if ~isfield(IOI.res,'seriesOK')
        %             disp(['No extracted time series available for subject ' int2str(SubjIdx) ' ... skipping series extraction']);
        %         else
        %if ~isfield(IOI.res,'meanOK') || job.force_redo
            [dir_ioimat dummy] = fileparts(job.IOImat{SubjIdx});
%             if isfield(IOI.ROI,'ROIfname')
%                 load(IOI.ROI.ROIfname)
%             end
            if save_figures
                dir_fig = fullfile(dir_ioimat,'fig');
                if ~exist(dir_fig,'dir'),mkdir(dir_fig);end
            end
            %get stimulation information - Careful, here onset duration is ignored!
            %loop over sessions
            for s1=1:length(IOI.sess_res)
                if all_sessions || sum(s1==selected_sessions)
                    if ~electro_stims %default stims
                        onsets_list{s1} = IOI.sess_res{s1}.onsets;
                    else
                        onsets_list{s1} = IOI.Sess(s1).U(which_ons).ons; %take first type of onsets
                    end
                end
            end
            
            %check if averaging over all selected sessions is sensible
            %in which case, global_M = 1 (M is the number of onset types)
            global_M = 1; found_first_M = 0;
            for s1=1:length(IOI.sess_res)
                if all_sessions || sum(s1==selected_sessions)
                    onsets_list{s1} = IOI.sess_res{s1}.onsets;
                    M{s1} = length(onsets_list{s1});
                    if ~found_first_M
                        possible_global_M = M{s1};
                        found_first_M = 1;
                    else
                        if ~(M{s1}==possible_global_M)
                            %no global M possible
                            global_M = 0;
                        end
                    end
                end
            end
            %compute mean over all sessions
            if global_M
                M = possible_global_M;
                %loop over onset types
                for m1=1:M %loop not used; which_ons
                    %loop over colors
                    for c1=1:length(IOI.color.eng)
                        %loop over ROIs
                        r2 = 0;
                        for r1=1:length(ROI)
                            if all_ROIs || sum(r1==selected_ROIs)
                                r2 = r2 + 1;
                                tmp_array_before = zeros(1,window_before);
                                tmp_array_after = zeros(1,window_after);
                                kb = 0; %counter of segments before onsets
                                ka = 0; %counter of segments after onsets
                                kb2 = 0; %counter of skipped segments before onsets
                                ka2 = 0; %counter of skipped segments after onsets
                                GSa{r2,m1}{c1} = [];
                                GSb{r2,m1}{c1} = [];
                                %loop over sessions
                                for s1=1:length(IOI.sess_res)
                                    if all_sessions || sum(s1==selected_sessions)
                                        tmp_d = ROI{r1}{s1,c1};
                                        if ~isempty(tmp_d)
                                            %loop over onsets for that session
                                            U = onsets_list{s1}{m1};
                                            for u1=1:length(U)
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
                                                    kb = kb+1;
                                                    GSb{r2,m1}{c1} = [GSb{r2,m1}{c1};tmp1-tmp_median];
                                                catch
                                                    kb2 = kb2+1;
                                                    if kb2 < 3
                                                        disp(['Could not include segment before onset ' int2str(u1) ' at time ' num2str(U(u1)) ...
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
                                                    tmp_array_after = tmp_array_after + tmp1-tmp_median;
                                                    ka = ka+1;
                                                    GSa{r2,m1}{c1} = [GSa{r2,m1}{c1};tmp1-tmp_median];
                                                catch
                                                    ka2 = ka2+1;
                                                    if ka2<3
                                                        disp(['Could not include segment after onset ' int2str(u1) ' at time ' num2str(U(u1)) ...
                                                            ' for session ' int2str(s1) ' for ROI ' int2str(r1) ...
                                                            ' for color ' int2str(c1) ' for onset type ' int2str(m1) ...
                                                            ' in global average over all sessions... skipping ' int2str(ka2) ' so far']);
                                                        
                                                    end
                                                end
                                                %                                                     if any(isnan(tmp_array_after))
                                                %                                                         a=1;
                                                %                                                     end
                                            end
                                        else
                                            if ~isfield(IOI.color,'contrast') || (isfield(IOI.color,'contrast') && ~(IOI.color.eng(c1)==IOI.color.contrast))
                                                disp(['Skipped session ' int2str(s1) ' for color ' int2str(c1) ' for ROI ' int2str(r1)]);
                                            end
                                        end
                                    end
                                end
                                GMb{r2,m1}{c1} = tmp_array_before/kb; %global mean before
                                GMa{r2,m1}{c1} = tmp_array_after/ka; %global mean after
                                GDa{r2,m1}{c1} = std(GSa{r2,m1}{c1},0,1)/sqrt(ka); %SEM
                                GDb{r2,m1}{c1} = std(GSb{r2,m1}{c1},0,1)/sqrt(kb); %SEM
                                if ka2>0 || kb2 > 0
                                    disp(['Skipped ' int2str(ka2) ' segments after onsets and ' int2str(kb2) ' segments before onsets']);
                                end
                            end
                        end
                    end
                end
            end
            
            %                 %loop over onset type
            %                 for m1=1:length
            %                     %loop over ROIs
            %                     for r1=1:length(ROI)
            %                         if all_ROIs || sum(r1==selected_ROIs)
            %                             %loop over sessions
            %                             for s1=1:length(IOI.sess_res)
            %                                 if all_sessions || sum(s1==selected_sessions)
            %                                     %loop over available colors
            %                                     for c1=1:length(IOI.color.eng)
            %
            %
            %                                     end
            %                                 end
            %                             end
            %                         end
            %                     end
            %
            
            
            
            IOI.res.meanOK = 1;
            IOI.res.GMa = GMa;
            IOI.res.GMb = GMb;
            IOI.res.GDa = GDa;
            IOI.res.GDb = GDb;
            %arrays GSa and GSb not saved...
            ctotal = [];
            h1 = 0;
            if generate_figures || save_figures
                %line specification - ROI x color
                lp1{1} = '-'; lp1{2} = ':'; lp1{3} = '--'; lp1{4} = '-.';
                if isfield(IOI.color,'HbO')
                    lp2{IOI.color.eng==IOI.color.HbO} = 'r'; %HbO
                    lp2{IOI.color.eng==IOI.color.HbR} = 'b'; %HbR
                    ctotal = [ctotal find(IOI.color.eng==IOI.color.HbO) ...
                        find(IOI.color.eng==IOI.color.HbR)];
                end
                %                     if isfield(IOI.color,'flow')
                %                         lp2{IOI.color.eng==IOI.color.flow} = 'k'; %Flow
                %                         ctotal = [ctotal find(IOI.color.eng==IOI.color.flow)];
                %                     end
                %                     if isfield(IOI.color,'contrast')
                %                         lp2{IOI.color.eng==IOI.color.contrast} = 'y'; %contrast (=flow up to rescaling)
                %                         ctotal = [ctotal find(IOI.color.eng==IOI.color.contrast)];
                %                     end
                
                % %                     %global figures, only for 1st onset
                % %                     if exist('GMa','var')
                % %                         if length(GMa) <= length(lp1) %plot identifiable series
                % %                             h1 = h1 + 1;
                % %                             h(h1) = figure;
                % %                             ls = linspace(0,job.window_after,window_after);
                % %                             for r1 = 1:length(GMa)
                % %                                 for c1 = ctotal
                % %                                     if ~add_error_bars
                % %                                         plot(ls,GMa{r1,1}{c1},[lp1{r1} lp2{c1}]); hold on
                % %                                     else
                % %                                         if r1 == 1
                % %                                             errorbar(ls(2:end-1),GMa{r1,1}{c1}(2:end-1),GDa{r1,1}{c1}(2:end-1),[lp1{r1} lp2{c1}]); hold on
                % %                                         else
                % %                                             plot(ls,GMa{r1,1}{c1},[lp1{r1} lp2{c1}]); hold on
                % %                                         end
                % %                                     end
                % %                                 end
                % %                             end
                % %                             if save_figures
                % %                                 filen = fullfile(dir_fig,['Mean_ROI.tiff']); %save as .tiff
                % %                                 print(h(h1), '-dtiffn', filen);
                % %                             end
                % %                             if ~generate_figures, close(h(h1)); end
                % %                         else %plot all series with random colors, skip error bars
                % %                             ColorSet = varycolor(size(GMa,1));
                % %                             for c1 = ctotal
                % %                                 h1 = h1 + 1;
                % %                                 h(h1) = figure;
                % %                                 ls = linspace(0,job.window_after,window_after);
                % %                                 %                                 a1 = zeros(size(GMa,1),size(GMa{r1,1}{c1},2));
                % %                                 %                                 for r1=1:size(GMa,1)
                % %                                 %                                     a1(r1,:) = GMa{r1,1}{c1};
                % %                                 %                                 end
                % %                                 set(gca, 'ColorOrder', ColorSet);
                % %                                 hold all
                % %                                 for r1=1:size(GMa,1)
                % %                                     plot(ls,GMa{r1,1}{c1}); %hold on
                % %                                 end
                % %                                 legend off
                % %                                 set(gcf, 'Colormap', ColorSet);
                % %                                 colorbar
                % %                                 if save_figures
                % %                                     filen = fullfile(dir_fig,['Mean_' IOI.color.eng(c1) '_allROI.tiff']); %save as .tiff
                % %                                     print(h(h1), '-dtiffn', filen);
                % %                                 end
                % %                                 if ~generate_figures, close(h(h1)); end
                % %                             end
                % %                         end
                % %                     else
                % %
                % %                     end
            end
            if isfield(job.IOImatCopyChoice,'IOImatCopy')
                newDir = job.IOImatCopyChoice.IOImatCopy.NewIOIdir;
                newDir = fullfile(dir_ioimat,newDir);
                if ~exist(newDir,'dir'),mkdir(newDir); end
                IOImat = fullfile(newDir,'IOI.mat');
            end
            save(IOImat,'IOI');
        %end
        
        toc
        disp(['Subject ' int2str(SubjIdx) ' complete']);
        out.IOImat{SubjIdx} = IOImat;
    catch exception
        disp(exception.identifier)
        disp(exception.stack(1))
        out.IOImat{SubjIdx} = job.IOImat{SubjIdx};
    end
end