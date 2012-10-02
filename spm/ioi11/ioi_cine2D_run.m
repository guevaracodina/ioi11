function out = ioi_cine2D_run(job)
%Step 1: Get various fields (one or a few functions)

[all_sessions selected_sessions] = ioi_get_sessions(job);
%filters
HPF = ioi_get_HPF(job);
LPF = ioi_get_LPF(job);
[shrinkage_choice SH] = ioi_get_shrinkage_choice(job);

if isfield(job.stim_choice,'manual_onsets')
    onset_time = job.stim_choice.manual_onsets.onset_time;
    manual_stim_choice = 1;
else
    manual_stim_choice = 0;
end
which_onset_type = job.which_onset_type;
group_onset_types = job.group_onset_types; %overrides which_onset_type
IC = job.IC;

normalize_choice = job.normalize_choice;
show_movie = job.show_movie;
%*************by Cong on 12/08/28
% interactive_mode = job.interactive_mode;
show_images=job.generate_images;
save_images=job.save_images;
%**************end
dir_cine = 'Cine_movies';
% high_limit = job.high_limit/100;
% low_limit = job.low_limit/100;
skip_overlap = job.skip_overlap;

%Step 2: loop over subjects:
%get size of windows before and after stimulation to average on, in data points
for SubjIdx=1:length(job.IOImat)
    try
        tic
        clear onsets_list M
        %Load IOI.mat information
        [IOI IOImat dir_ioimat]= ioi_get_IOI(job,SubjIdx);
        
        if ~isfield(IOI.res,'cineOK') || job.force_redo
            cineDir = fullfile(dir_ioimat,dir_cine);
            if ~exist(cineDir,'dir'), mkdir(cineDir); end
            %get stimulation information - Careful, here onset duration is ignored!
            if IC.include_HbT
                if ~isfield(IOI.color,'HbT')
                    IOI.color.HbT = 'T';
                    IOI.color.eng = [IOI.color.eng IOI.color.HbT];
                end
            end
            if ~isfield(IOI,'dev')
                IOI.dev.TR = 0.2;
            end
            
            %Step 3: another function: get windows 
            WI = ioi_get_windows(IOI.dev.TR,job);
            
            
            %to allow at least one data point after window_before
            if isfield(job.stim_choice,'manual_onsets')
                if onset_time == 0
                    onset_time = 2*IOI.dev.TR;
                end
            end
            
            %save shrunk images
            if shrinkage_choice
                IOI = ioi_save_shrunk_images(IOI,job,SH,dir_ioimat);
            end
            %Overwrite IOI
            save(IOImat,'IOI');
            
            %loop over sessions
            for s1=1:length(IOI.sess_res);
                if all_sessions || sum(s1==selected_sessions)
                    if manual_stim_choice
                        onsets_list{s1} = {onset_time};
                    else
                        onsets_list{s1} = IOI.sess_res{s1}.onsets;
                    end
                    
                    if group_onset_types
                        tmp = [];
                        for m1=1:length(onsets_list{s1})
                            %if any(m1==which_onset_type)
                            tmp = [tmp onsets_list{s1}{m1}];
                            %    end
                        end
                        onsets_list{s1} = {};
                        onsets_list{s1}{1} = sort(tmp);
                    end
                    
                    %Function to skip_overlapping onsets
                    [onsets_list{s1} skipped] = ioi_skip_overlapping(skip_overlap, onsets_list{s1},...
                        job.which_onset_type,job.window_after,job.window_before);
                    %remove onsets that are too close together
                    %count onsets
                    cnt = [];
                    for m1=1:length(onsets_list{s1})
                        if any(m1==which_onset_type)
                            cnt = [cnt length(onsets_list{s1}{m1})];
                        end
                    end
                    disp(['Onset counts, session ' int2str(s1)]);
                    disp(cnt)
                    IOI.cine_onset_count{s1} = cnt;
                    IOI.cine_onsets_list{s1} = onsets_list{s1};
                    IOI.cine_skipped{s1} = skipped;
                    disp(['Skipped:' int2str(skipped)]);
                    %loop over colors
                    interactive_mode_all_color=0;
                    for c1=1:length(IOI.color.eng);
                        %Get the proper do_color function 
%                       doColor = ioi_doColor(IOI,c1,IC)
                        do_color = 0;
                        if IC.include_flow
                            if isfield(IOI.color,'flow')
                                if IOI.color.eng(c1)==IOI.color.flow
                                    do_color = 1;
                                end
                            end
                        end
                        if IC.include_OD
                            if any(IOI.color.eng(c1)==[IOI.color.green IOI.color.red IOI.color.yellow])
                                do_color = 1;
                            end
                        end
                        if isfield(IOI.color,'HbO')
                            if IC.include_HbT
                                HbColors = [IOI.color.HbO IOI.color.HbR IOI.color.HbT];
                            else
                                HbColors = [];
                            end
                            if IC.include_HbO
                                if IOI.color.eng(c1)==IOI.color.HbO
                                    do_color = 1;
                                end
                            end
                            if IC.include_HbR
                                if IOI.color.eng(c1)==IOI.color.HbR
                                    do_color = 1;
                                end
                            end
                            if any(IOI.color.eng(c1)==HbColors)
                                do_color = 1;
                            end
                        end
                        if do_color
                            m2 = 0;
                            %loop over onset types
                            for m1=1:length(onsets_list{s1})
                                if any(m1==which_onset_type)
                                    m2 = m2+1;
                                    arrays_assigned = 0;
                                    kb = 0; %counter of segments before onsets
                                    ka = 0; %counter of segments after onsets
                                    kb2 = 0; %counter of skipped segments before onsets
                                    ka2 = 0; %counter of skipped segments after onsets
                                    %loop over onsets for that session
                                    U = round(onsets_list{s1}{m1}/IOI.dev.TR)-WI.window_offset; %in data points
                                    
                                    %add a function for averaging
                                    
                                    for u1=1:length(U)
                                        %Frames to get:
                                        fr_before = U(u1)-WI.window_before:U(u1)-1;
                                        fr_after = U(u1):U(u1)+WI.window_after-1;
                                        %Load images before and after stimulus
                                        Yb = ioi_get_images(IOI,fr_before,c1,s1,dir_ioimat,shrinkage_choice);
                                        Ya = ioi_get_images(IOI,fr_after,c1,s1,dir_ioimat,shrinkage_choice);
                                        if ~isempty(Ya) && ~isempty(Yb)
                                            if ~arrays_assigned %assign arrays
                                                tmp_array_before = zeros(size(Yb));
                                                tmp_array_after = zeros(size(Ya));
                                                arrays_assigned = 1;
                                            end
                                            
                                            switch normalize_choice
                                                case 1
                                                    tmp_median = median(Yb,3);
                                                case 2
                                                    tmp_median = squeeze(Yb(:,:,end));
                                                case 3
                                                    tmp_median = mean(Yb,3);
                                                case 4
                                                    tmp_median = zeros(size(Yb,1),size(Yb,2));
                                            end
                                            tmp_median_b = repmat(tmp_median,[1 1 size(Yb,3)]);
                                            tmp_median_a = repmat(tmp_median,[1 1 size(Ya,3)]);
                                            tmp_array_before = tmp_array_before + Yb - tmp_median_b;
                                            tmp_array_after = tmp_array_after + Ya - tmp_median_a;
                                            ka = ka+1;
                                            kb = kb+1;
                                        else
                                            if isempty(Ya)
                                                ka2 = ka2+1;
                                                if ka2<3
                                                    disp(['Could not include segment due to incomplete data after onset ' int2str(u1) ' at time ' num2str(U(u1)*IOI.dev.TR) ...
                                                        ' for session ' int2str(s1) ...
                                                        ' for color ' int2str(c1) ' for onset type ' int2str(m1) ...
                                                        '... skipping ' int2str(ka2) ' so far']);
                                                    
                                                end
                                            end
                                            if isempty(Yb)
                                                kb2 = kb2+1;
                                                if kb2 < 3
                                                    disp(['Could not include segment due to incomplete data before onset ' int2str(u1) ' at time ' num2str(U(u1)*IOI.dev.TR) ...
                                                        ' for session ' int2str(s1) ...
                                                        ' for color ' int2str(c1) ' for onset type ' int2str(m1) ...
                                                        '... skipping ' int2str(kb2) ' so far']);
                                                end
                                            end
                                        end
                                    end
                                    %divide by number of found segments (ka = kb always in this module)
                                    %save movie
                                    d = tmp_array_after/ka;
                                    time=0;
                                    if job.downFact > 1
                                        nT = round(size(d,3)/job.downFact);
                                        tmpd = zeros(size(d,1),size(d,2),nT); 
                                        interactive_mode = 1;
                                        for j0=1:nT
                                            time=j0*job.downFact*0.2-job.window_offset;  %**************caculate the time for each frame
                                            tmpd(:,:,j0) = mean(d(:,:,(1:job.downFact)+job.downFact*(j0-1)),3);
                                            %*******by Cong on 12/08/28
                                            if show_images
                                                Images=figure; imagesc(tmpd(:,:,j0)); colorbar;
                                                tit=['Cine images' gen_num_str(s1,2) ' ' IOI.color.eng(c1) ' onset' int2str(m1)];
                                                title(tit);
                                                close(Images);
                                            end                                           
                                            if save_images
                                              
                                                if interactive_mode_all_color
                                                    clims = [clims1 clims2];
                                                    Images=figure;                                                    
                                                    imagesc(tmpd(:,:,j0),clims);
                                                    colorbar;
                                                    tit1=['Cine images' gen_num_str(s1,2) ' ' IOI.color.eng(c1) ' onset' int2str(m1) ' frame' int2str(j0) ' time ' int2str(time) 's'];
                                                    title(tit1);
                                                    tit=['Cine_images' gen_num_str(s1,2) '_' IOI.color.eng(c1) '_onset' int2str(m1) '_frame' int2str(j0) '_time_' int2str(time) 's'];
                                                    filename=[cineDir,tit];
                                                    print(Images, '-dpng', [filename '.png'], '-r300');
                                                    %saveas(Images,filename, 'tif');
                                                    saveas(Images,filename, 'fig');
                                                    close(Images);
                                                else

                                                if j0 == 1
                                                    keep_going = 1;
                                                    first_pass = 1;
                                                end
                                                other_first_pass = 1;
                                                %                                                 while (keep_going && (interactive_mode || first_pass ))|| (other_first_pass && j0 > 1)||((~keep_going )&& (interactive_mode))                                                
                                                flag=0;
                                                flag_first_pass = 0;
                                                while (keep_going && (interactive_mode || first_pass || other_first_pass) ) || (other_first_pass && j0 > 1)||flag
                                                    Images=figure;
                                                    if first_pass
                                                        imagesc(tmpd(:,:,j0));
                                                        first_pass = 0;
                                                        flag_first_pass = 1;
                                                    else
                                                        imagesc(tmpd(:,:,j0),clims);
                                                        other_first_pass = 0;
                                                    end
                                                    colorbar;
                                                    tit1=['Cine images' gen_num_str(s1,2) ' ' IOI.color.eng(c1) ' onset' int2str(m1) ' frame' int2str(j0) ' time ' int2str(time) 's'];
                                                    title(tit1);
                                                    tit=['Cine_images' gen_num_str(s1,2) '_' IOI.color.eng(c1) '_onset' int2str(m1) '_frame' int2str(j0) '_time_' int2str(time) 's'];
                                                    filename=[cineDir,tit];
                                                    print(Images, '-dpng', [filename '.png'], '-r300');
                                                    %saveas(Images,filename, 'tif');
                                                    saveas(Images,filename, 'fig');
                                                    if keep_going && interactive_mode
                                                        keep_going = spm_input('Change colorbar limits ',1,'y/n',[1 0]);
                                                        if keep_going
                                                            clims1 = spm_input('Enter min value ','+1');
                                                            clims2 = spm_input('Enter max value ','+1');
                                                            clims = [clims1 clims2];
                                                        end
                                                        interactive_mode = ~spm_input('Leave interactive mode for this color ','+1','y/n',[1 0]);
                                                        interactive_mode_all_color = ~spm_input('Leave interactive mode for all colors ','+1','y/n',[0 1]);
                                                    end 
                       
                                                    if keep_going && (~interactive_mode) && flag_first_pass
                                                        flag=flag+1;
                                                        if flag>1
                                                            flag=0;
                                                        end
                                                        
                                                        flag_first_pass = 0;
                                                    else
                                                        flag = 0;
                                                    end
%                                                     if ~interactive_mode_all_color
%                                                        clims = [clims1 clims2];
%                                                     end
                                                    close(Images);
                                                end
                                                end
                                            end
                                            %************end
                                        end
                                        d = tmpd;
                                    end
                                    m0 = min(d(:));
                                    M0 = max(d(:));
                                    %best is to save d
                                    save(fullfile(cineDir,['cine_S' gen_num_str(s1,2) '_' IOI.color.eng(c1) '_onset' int2str(m1) '.Mdat']),'d','-v7.3');
                                    dM = M0-m0;
                                    Nf = size(d,3);
                                    F(Nf) = struct('cdata',[],'colormap',[]);
                                    h0 = figure;
                                    %if IOI.color.eng(c1) == IOI.color.HbR
                                    %    clims = [m0+dM*(1-5*high_limit) M0+dM*(1-low_limit)];
                                    %else
%                                     clims = [m0+dM*low_limit M0+dM*high_limit];
                                    %end
                                    try
                                        try
                                            fname_movie = fullfile(cineDir,['cine_S' gen_num_str(s1,2) '_' IOI.color.eng(c1) '_onset' int2str(m1) '.avi']);
                                            vidObj = VideoWriter(fname_movie);
                                            open(vidObj);
                                            VideoOK = 1;
                                        catch
                                            VideoOK = 0;
                                            fname_movie = fullfile(cineDir,['cine_S' gen_num_str(s1,2) '_' IOI.color.eng(c1) '_onset' int2str(m1) '.mat']);
                                        end
                                        %set(gca,'NextPlot','replacechildren');
                                        for i0=1:size(d,3); %20
                                            try
                                                imagesc(squeeze(d(:,:,i0)),clims);
                                            catch %in case there are Inf in d
                                                imagesc(squeeze(d(:,:,i0)));
                                            end
                                            F(i0) = getframe;
                                            if VideoOK
                                                writeVideo(vidObj,F(i0));
                                            end
                                        end
                                        if VideoOK
                                            close(vidObj);
                                        else
                                            save(fname_movie,'F','-v7.3');
                                        end
                                        %movie(h0,F,1,1/IOI.dev.TR); %,[0 0 0 0]);
                                        try close(h0); end
                                        IOI.cine{s1,m1}{c1}.fname_movie = fname_movie;
                                    end
                                end
                            end
                            if (ka2>0 || kb2 > 0)
                                disp(['Skipped ' int2str(ka2+kb2) ' segments']);
                            end
                        end
                    end
                end
                
            end
            IOI.res.cineOK = 1;
            save(IOImat,'IOI');
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
end