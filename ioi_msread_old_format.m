function IOI = ioi_msread_old_format(IOI,job)
try
    if isfield(job.treatment_mode,'expeditive_mode')
        expedite = 1;
    else
        expedite = 0;
    end
    %recall dir structure
    subj_name = IOI.subj_name;
    dir_subj_raw = IOI.dir.dir_subj_raw;
    dir_subj_res = IOI.dir.dir_subj_res;
    PartialRedo2 = job.PartialRedo2;
    subj_OK = IOI.subj_OK; %always 1 at this stage
    block_size = 500; %number of images per block, if save_choice = 2 selected
    %Hard coded the 6 frame steps for different colors, not sure it is always true
    frames_per_cycle = 6;
    %Color names - could be put in user interface?
    str_anat = 'G';
    str_laser = 'L'; %Optional
    str_red = 'R'; %Varying number of colors allowed
    str_green = 'G';
    str_yellow = 'Y';
    str_contrast = 'C'; %Optional - proportional to flow, derived from laser speckle
    %If adding or removing colors, make sure the string of English and French
    %colors have the same length and are a one-to-one mapping (code loops over colors
    %to find images):
    str_color = [str_red str_green str_yellow str_laser str_contrast]; %red, green, yellow and laser speckle
    str_color_french = [str_red 'V' 'J' str_laser str_contrast];
    OD_label = 'OD'; %label for optical density images
    suffix_for_anat_file = 'anat'; %to build anatomical image name
    sess_label = 'Sess'; %prefix for name of directories for each session
    short_el_label = 'el'; %short name for output electro file
    %leave voxel size in arbitrary units for now, for anatomical image
    vx_anat = [1 1 1];
    %Variable changed meaning here - more convenient for user to specify
    %duration in seconds
    min_session_duration = job.sess_min_image_files;
    if isfield(job,'acq_freq')
        temp_TR = job.color_number/job.acq_freq;
    else
        temp_TR = 0.2;
    end
    temp_ImNum = 80;
    nzero_padding = 5;
    
    [shrinkage_choice SH] = ioi_get_shrinkage_choice(job);
    %select a subset of sessions
    [all_sessions selected_sessions] = ioi_get_sessions(job);
    
    %choose saving mode
    try
        save_choice = job.save_choice;
    catch
        save_choice = 1;
    end
    forceProcessingOn = job.forceProcessingOn;
    
    if ~job.PartialRedo2
        %partially fill IOI structure
        IOI.save_choice = save_choice;
        switch save_choice
            case 1
                IOI.save_choice_method = 'One file per session';
            case 2
                IOI.save_choice_method = 'One file per block';
            case 3
                IOI.save_choice_method = 'One file per image (frame)';
        end
        %Store color information
        IOI.color.eng = str_color;
        IOI.color.french = str_color_french;
        IOI.color.anat = str_anat;
        IOI.color.laser = str_laser;
        IOI.color.red = str_red;
        IOI.color.green = str_green;
        IOI.color.yellow = str_yellow;
        IOI.color.contrast = str_contrast;
    end
    %resume files - subject level
    [files_txt,dummy] = cfg_getfile('FPList',dir_subj_raw,'.txt');
    if length(files_txt) == 2 %Must have exactly two text files
        %careful, the order of the files here could change if their
        %names were to be changed
        IOI.info.expt = files_txt{1,:}; %resume_manip
        IOI.info.comments = files_txt{2,:};%resumecommentaire
    else
        disp(['Wrong number of text files at subject level in ' dir_subj_raw]);
        subj_OK = 0;
    end
    [dummy,dirs] = cfg_getfile('FPListRec',dir_subj_raw,'.bin');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %1- check that raw directory structure is consistent
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if job.PartialRedo2
        sess_raw = IOI.sess_raw; %Keep previously found sessions
    else
        if ~expedite
            sess_raw = {};
            for i=1:(length(dirs)/2)
                sess = [];
                %Boolean to decide whether a session can be kept
                sess_OK = 1;
                %skip sessions not selected by user - careful, as this
                %session order corresponds to raw session order, while
                %in later modules, session order corresponds to
                %processed session order
                if all_sessions || sum(i==selected_sessions)
                    if ~strcmp([dirs{2*i-1,:} '_images'],dirs{2*i,:})
                        IOI = disp_msg(IOI,['Missing directory: ' dirs{2*i,:}]);
                    else
                        %check that text files are present
                        [files_txt,dummy] = cfg_getfile('FPListRec',dirs{2*i-1,:},'.txt');
                        if length(files_txt) == 3
                            %Three files are required for each session,
                            %and they need to be found and assigned in
                            %the following order:
                            sess.Frame = files_txt{1,:};
                            sess.TTL = files_txt{2,:}; %TTLseuil
                            sess.info = files_txt{3,:};
                            [files_bin,dummy] = cfg_getfile('FPListRec',dirs{2*i,:},'.bin');
                            sess.images = files_bin; %all images
                            for c1=1:length(str_color)
                                if sess_OK
                                    %str1 = str_color(c1);
                                    strf = str_color_french(c1);
                                    %separate images into various colors
                                    [files_bin,dummy] = cfg_getfile('FPListRec',dirs{2*i,:},[strf strf]);
                                    sess.fnames{c1} = files_bin;
                                    if c1 ==1
                                        image_count = length(files_bin);
                                        %approximate computation of length
                                        %of session in seconds
                                        if image_count*temp_ImNum*temp_TR < min_session_duration
                                            sess_OK = 0;
                                            IOI = disp_msg(IOI,['Insufficient number of images for ' dirs{2*i-1,:} ' ...skipping.']);
                                        end
                                    else
                                        %weaker criterion - number of image
                                        %files might differ by 1 because
                                        %one only one extra image
                                        if abs(length(files_bin)- image_count)>= 1
                                            if ~strcmp(strf,str_contrast) %OK if contrast files missing
                                                sess_OK = 0;
                                                IOI = disp_msg(IOI,['Problem with number of ' strf strf ' images for ' dirs{2*i-1,:}]);
                                            end
                                        end
                                    end
                                end
                            end
                            [files_bin,dummy] = cfg_getfile('FPListRec',dirs{2*i,:},'electro');
                            sess.electro = files_bin;
                            if ~isempty(files_bin)
                                sess.electroOK = 1;
                            else
                                sess.electroOK = 0;
                                IOI = disp_msg(IOI,['Problem with electro file for ' dirs{2*i-1,:}]);
                            end
                            %add found session to session list
                            if forceProcessingOn
                                sess_OK = 1; %ignore warnings, and include session anyway
                            end
                            if sess_OK, sess_raw = [sess_raw sess]; end
                        else
                            IOI = disp_msg(IOI,['Problem with text files in ' dirs{2*i,:}]);
                        end
                    end
                end
            end
        else
            sess_raw = {}; sess_OK = 1;
            for i=1:(length(dirs)/2)
                [files_bin,dummy] = cfg_getfile('FPListRec',dirs{2*i,:},'electro');
                sess.electro = files_bin;
                if ~isempty(files_bin)
                    sess.electroOK = 1;
                else
                    sess.electroOK = 0;
                    IOI = disp_msg(IOI,['Problem with electro file for ' dirs{2*i-1,:}]);
                end
                if sess_OK, sess_raw = [sess_raw sess]; end
            end
        end
        if ~isempty(sess_raw)
            IOI.sess_raw = sess_raw;
        else
            disp(['No sessions found... skipping subject ' int2str(job.SubjIdx) ' (namely: ' job.top_bin_dir{job.SubjIdx} ')']);
            subj_OK = 0;
            %return
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %2- Extract various info and stimuli if available
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if subj_OK
        %extracted info and stimuli, by session
        if PartialRedo2
            sess_res = IOI.sess_res;
        else
            sess_res  = {};
        end
        for s1=1:length(sess_raw)
            %we normally process all the images in sess_raw, unless we are
            %recalculating a subset of the sessions, using option
            %Partial in force_redo, leading to PartialRedo2 == 1.
            if ~expedite
                if ~PartialRedo2 || all_sessions || sum(s1==selected_sessions)
                    % First read all information concerning the experiment
                    try
                        [scan_info, physio]=ioi_extract_info(IOI.info.expt,...
                            sess_raw{s1}.info,sess_raw{s1}.TTL,sess_raw{s1}.Frame);
                    catch exception
                        disp(exception.identifier)
                        disp(exception.stack(1))
                        IOI = disp_msg(IOI,['Problem with raw session ' int2str(i) ': could not extract info']);
                        disp('Session will be kept, but this could lead to various problems later');
                        disp('Best is to find out the problem with this session, and fix it or remove it from the raw data folder');
                    end
                    if s1==1 %assume acquisition frequency is unchanged for later sessions
                        % Identifiy times for each frame and each color separately
                        %try
                        IOI.dev.acq_freq_hz=(scan_info.Frame(end,1)-scan_info.Frame(1,1))/scan_info.Frame(end,3);
                        % Fix TR for everyone, the factor 6 comes from the fact that we need 6
                        % camera frames to return to the same colors
                        IOI.dev.TR=frames_per_cycle/IOI.dev.acq_freq_hz;
                        %                         catch
                        %                             IOI.dev.acq_freq_hz = 5;
                        %                             IOI.dev.TR=0.1999;
                        %                         end
                    end
                    % 6 images make one frame, each frame is indexed, partial frames are
                    % counted as one frame - this is only approximate
                    n_frames = ceil((scan_info.Frame(end,1)-scan_info.Frame(1,1)+1)/frames_per_cycle);
                    sess = [];
                    sess.scan_info = scan_info;
                    sess.physio = physio;
                    sess.n_frames = n_frames;
                    %stimuli information
                    list_stim = {};
                    %for each type of onsets
                    for stim_index=1:size(scan_info.stim3,1)
                        stim = [];
                        %right hand side of == is number representing onset type
                        %column 1 of stim2 is frame number sight after stimulation
                        %select entries in column 1 of stim2
                        %last column of stim1 is
                        tmp = scan_info.stim1(:,end)==scan_info.stim3(stim_index,1);
                        stim.onset_frame=scan_info.stim2(...
                            tmp(1:size(scan_info.stim2,1)),1);
                        stim.onset_frame=ceil(stim.onset_frame/frames_per_cycle);
                        %Next 4 entries: Approximate info - do not use (??):
                        stim.intensity=scan_info.stim3(stim_index,2);
                        stim.train_period=scan_info.stim3(stim_index,3);
                        stim.n_stim_in_train=scan_info.stim3(stim_index,4);
                        stim.train_duration=scan_info.stim3(stim_index,5);
                        list_stim = [list_stim stim];
                    end
                    clear names onsets durations parameters
                    %Converted to seconds, rather than frame number
                    for stim_index=1:length(list_stim)
                        names{stim_index}=['Stim_',num2str(stim_index)];
                        onsets{stim_index}=list_stim{stim_index}.onset_frame'*IOI.dev.TR;
                        durations{stim_index}=list_stim{stim_index}.train_duration;
                        parameters{stim_index}=stim_index;
                    end
                    %Store onset information
                    sess.list_stim = list_stim;
                    try
                        sess.names = names; %not quite the SPM format, need to put in Sess.U format
                        sess.onsets = onsets; %in seconds
                        sess.durations = durations; %in seconds
                        sess.parameters = parameters;
                    catch
                        disp('No onsets found');
                    end
                    %frames for each color
                    for c1=1:length(str_color)
                        sess.camera{c1}=scan_info.Frame(scan_info.FrameCouleur==str_color_french(c1),1);
                        sess.fname{c1} = {}; %for later, to add nifti file names of images
                    end
                    if PartialRedo2
                        sess_res{s1} = sess;
                    else
                        sess_res =[sess_res sess];
                    end
                end
            end
            IOI.sess_res = sess_res;
            IOI.res.shrinkageOn = shrinkage_choice;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %3- Anatomical image
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %use first green colored image as the anatomy
        if ~expedite
            image_total=ioi_read_single_bin_image_old_format(sess_raw{1}.fnames{str_anat == str_color}{1,:});
            fname = fullfile(dir_subj_res, [subj_name '_' suffix_for_anat_file]);
            ioi_save_images(image_total, fname, vx_anat,[],'Anatomical image')
            %ioi_save_nifti(image_total, fname, vx_anat);
            IOI.res.file_anat=[fname '.nii'];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %4- Functional images and electrophysiology
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % For each file to treat, reassemble and save to nifti
        % Optical density saved as this is more convenient for further processing
        % It is normalized by median in time direction
        %loop over sessions
        for s1=1:length(sess_raw)
            if ~PartialRedo2 || all_sessions || sum(s1==selected_sessions)
                %tic
                if ~expedite
                    %create dir for this session
                    if s1<10, str_s1 = '0'; else str_s1 = ''; end
                    str_s1 = [str_s1 int2str(s1)];
                    sess_str = [sess_label str_s1];
                    dir_sess_res = fullfile(dir_subj_res,sess_str);
                    if ~exist(dir_sess_res,'dir'),mkdir(dir_sess_res); end
                    n_frames = sess_res{s1}.n_frames;
                    %list (to be constructed) of colors available
                    hasRGY = [];
                    %loop over colors
                    for c1=1:length(str_color)
                        need_resize_message = 1;
                        str1 = str_color(c1);
                        if ~(str1==str_contrast && isempty(sess_raw{s1}.fnames{c1})) %exclude when contrast images are absent
                            fcount=0; %image counter (of actual frames recorded)
                            im_count = 0; %true frame counter (all images including interpolated ones)
                            %loop over image files for this session
                            image_total = [];
                            %Do NOT shrink laser speckle! - this will
                            %be done later, in the flow calculation module
                            vx = [1 1 1];
                            if ~(str1 == str_laser) %OK to shrink contrast images though
                                if shrinkage_choice
                                    %Keep same whatever value of PartialRedo2
                                    IOI.res.shrink_x=SH.shrink_x;
                                    IOI.res.shrink_y=SH.shrink_y;
                                    vx=[SH.shrink_x SH.shrink_y 1];
                                end
                            end
                            %loop over binary files for this session and color
                            for f1=1:length(sess_raw{s1}.fnames{c1})
                                %need to keep track of special cases when missing
                                %indices are the first index of the first file or the
                                %last index of the last file
                                if f1 == 1
                                    first_file = 1;
                                else
                                    first_file = 0;
                                end
                                if f1 == length(sess_raw{s1}.fnames{c1})
                                    last_file = 1;
                                else
                                    last_file = 0;
                                end
                                %Index for color to use in case of contrast
                                if (str1 == str_contrast)
                                    tc1 = find(str_color==str_laser);
                                else
                                    tc1 = c1;
                                end
                                %Read one binary file at a time (typically contains about 80 images)
                                [image_part fcount]=ioi_read_bin_image_old_format(sess_raw{s1}.fnames{c1}{f1,:}, ...
                                    fcount,ceil(sess_res{s1}.camera{tc1}/frames_per_cycle),vx,n_frames,first_file,last_file,im_count);
                                %image_part is of shape [nx, ny, 1, time] and has been
                                %interpolated for missing frames
                                [nx ny nz nt] = size(image_part); %Here nz = 1, always
                                im_count = im_count+nt;
                                if f1==1
                                    nx0 = nx;
                                    ny0 = ny;
                                end
                                if f1==1
                                    %create temporary file
                                    fid = fopen('tmp_file.dat','w');
                                    fwrite(fid, zeros([nx ny 1 n_frames]), 'single');
                                    fclose(fid);
                                    im_obj = memmapfile('tmp_file.dat','format',...
                                        {'single' [nx ny 1 n_frames] 'image_total'},...
                                        'writable',true);
                                end
                                %catch case when there is a resizing of the images during acquisition
                                %Ideally, this should never happen, but it has
                                if f1>1
                                    if ~(nx==nx0 && ny==ny0)
                                        if need_resize_message %Boolean, to avoid repeating warning message
                                            need_resize_message = 0; %first pass
                                            warning_message = ['Image resizing required for session ' ...
                                                int2str(s1) ' file ' int2str(f1) ' color ' int2str(c1)];
                                            IOI = disp_msg(IOI,warning_message);
                                        end
                                        image_part = ioi_imresize(image_part,1,nx0,ny0,1,1);
                                    end
                                end
                                %Build large image_total array, for later calculating the median
                                im_obj.Data.image_total(:,:,:,im_count-nt+1:im_count) = image_part;
                            end %end of loop over binary files
                            if ~(im_count==n_frames)
                                disp('Problem with image interpolation - wrong final number');
                                IOI = disp_msg(IOI,['Session ' int2str(s1) ',color ' str1 ...
                                    'n_frames = ' int2str(n_frames) '; final number images = ' int2str(im_count)]);
                            end
                            %check that color order is OK
                            [sts i0] = ioi_check_color_order(im_obj.Data.image_total,str1,str_laser);
                            IOI.bad_frames{s1,c1} = i0;
                            if ~sts
                                try
                                    warning_message = ['Possible problem with color order for frame ' int2str(i0.bfr(1)) ...
                                        ' and ' int2str(length(i0.bfr)-1) ' other frames, for session ' int2str(s1) ', for color ' int2str(c1)];
                                    IOI = disp_msg(IOI,warning_message);
                                end
                            end
                            %calculate median for whole image in time direction
                            image_median = median(im_obj.Data.image_total,4);
                            if ~(str1 == str_laser || str1 == str_contrast) %not laser speckle nor contrast
                                im_obj.Data.image_total = -log(im_obj.Data.image_total./repmat(image_median,[1 1 1 n_frames]));
                            end
                            
                            %Save the median
                            IOI.sess_res{s1}.fname_median{c1} = fullfile(dir_subj_res,sess_str, ...
                                [subj_name '_' OD_label '_median_' str1 '_' sess_str]);
                            %ioi_save_nifti(single(image_median),IOI.sess_res{s1}.fname_median{c1},vx);
                            tit0 = [subj_name ' ' OD_label ' median ' str1 ' ' sess_str];
                            ioi_save_images(single(image_median),IOI.sess_res{s1}.fname_median{c1},vx,[],tit0);
                            try
                                min_image = min(im_obj.Data.image_total,[],4);
                                max_image = max(im_obj.Data.image_total,[],4);
                                tenthpctle_image = prctile(im_obj.Data.image_total,10,4);
                                ninetiethpctle_image = prctile(im_obj.Data.image_total,90,4);
                                change = single(max_image) ./single(min_image);
                                change_90_10 = single(ninetiethpctle_image) ./ single(tenthpctle_image);
                                sess.fname_min{c1} = fullfile(dir_subj_res,sess_str, ...
                                    [subj_name '_' OD_label '_min_' str1 '_' sess_str]);
                                sess.fname_max{c1} = fullfile(dir_subj_res,sess_str, ...
                                    [subj_name '_' OD_label '_max_' str1 '_' sess_str]);
                                sess.fname_10pctle{c1} = fullfile(dir_subj_res,sess_str, ...
                                    [subj_name '_' OD_label '_10pctle_' str1 '_' sess_str]);
                                sess.fname_90pctle{c1} = fullfile(dir_subj_res,sess_str, ...
                                    [subj_name '_' OD_label '_90pctle_' str1 '_' sess_str]);
                                sess.fname_change{c1} = fullfile(dir_subj_res,sess_str, ...
                                    [subj_name '_' OD_label '_change_' str1 '_' sess_str]);
                                sess.fname_change_90_10{c1} = fullfile(dir_subj_res,sess_str, ...
                                    [subj_name '_' OD_label '_change_90_10_' str1 '_' sess_str]);
                                tit1 = [subj_name ' ' OD_label ' min ' str1 ' ' sess_str];
                                tit2 = [subj_name ' ' OD_label ' max ' str1 ' ' sess_str];
                                tit3 = [subj_name ' ' OD_label ' 10pctle ' str1 ' ' sess_str];
                                tit4 = [subj_name ' ' OD_label ' 90pctle ' str1 ' ' sess_str];
                                tit5 = [subj_name ' ' OD_label ' change min to max ' str1 ' ' sess_str];
                                tit6 = [subj_name ' ' OD_label ' change 90 to 10 ' str1 ' ' sess_str];
                                %                             ioi_save_nifti(single(min_image),sess.fname_min{c1},vx);
                                %                             ioi_save_nifti(single(max_image),sess.fname_max{c1},vx);
                                %                             ioi_save_nifti(single(change),sess.fname_change{c1},vx);
                                %                             ioi_save_nifti(single(tenthpctle_image),sess.fname_10pctle{c1},vx);
                                %                             ioi_save_nifti(single(ninetiethpctle_image),sess.fname_90pctle{c1},vx);
                                %                             ioi_save_nifti(single(change_90_10),sess.fname_change_90_10{c1},vx);
                                ioi_save_images(single(min_image),sess.fname_min{c1},vx,[],tit1);
                                ioi_save_images(single(max_image),sess.fname_max{c1},vx,[],tit2);
                                ioi_save_images(single(change),sess.fname_change{c1},vx,[],tit3);
                                ioi_save_images(single(tenthpctle_image),sess.fname_10pctle{c1},vx,[],tit4);
                                ioi_save_images(single(ninetiethpctle_image),sess.fname_90pctle{c1},vx,[],tit5);
                                ioi_save_images(single(change_90_10),sess.fname_change_90_10{c1},vx,[],tit6);
                            end
                            %vx currently not used in ioi_save_nifti
                            %and ioi_write_nifti, despite apparent dependency
                            switch save_choice
                                case 1
                                    fname = fullfile(dir_subj_res,sess_str, [subj_name '_' OD_label '_' str1 '_S' str_s1 '.nii']);
                                    IOI.sess_res{s1}.fname{c1} = [IOI.sess_res{s1}.fname{c1}; fname];
                                    ioi_save_nifti(im_obj.Data.image_total,fname,vx);
                                case 2
                                    %save per block of size size_block
                                    for j1=1:block_size:n_frames
                                        if j1 <= n_frames-block_size
                                            last_index = j1+block_size-1;
                                        else
                                            last_index = n_frames;
                                        end
                                        first_image_str = gen_num_str(j1,nzero_padding);
                                        last_image_str = gen_num_str(last_index,nzero_padding);
                                        fname = fullfile(dir_subj_res,sess_str, [subj_name '_' OD_label '_' str1 '_S' str_s1 '_' first_image_str '_to_' last_image_str  '.nii']);
                                        IOI.sess_res{s1}.fname{c1} = [IOI.sess_res{s1}.fname{c1}; fname];
                                        ioi_save_nifti(im_obj.Data.image_total(:,:,:,j1:last_index), fname, vx);
                                    end
                                case 3
                                    %save each frame separately
                                    for j1=1:n_frames
                                        image_str = gen_num_str(j1,nzero_padding);
                                        fname = fullfile(dir_subj_res,sess_str, [subj_name '_' OD_label '_' str1 '_S' str_s1 '_' image_str '.nii']);
                                        IOI.sess_res{s1}.fname{c1} = [IOI.sess_res{s1}.fname{c1}; fname];
                                        ioi_save_nifti(im_obj.Data.image_total(:,:,:,j1),fname,vx);
                                    end
                            end
                            clear im_obj;
                            delete('tmp_file.dat');
                            %Add found color - used later in concentration calculation module
                            if ~(str1 == str_laser || str1 == str_contrast)
                                hasRGY = [hasRGY str1];
                            end
                            disp(['Done processing session: ' int2str(s1) ', color: ' str1]);
                        end
                    end %for c1
                    IOI.sess_res{s1}.hasRGY = hasRGY;
                end
                %electrophysiology
                if isfield(IOI.sess_raw{s1},'electroOK')
                    IOI = extract_electro(IOI,s1,IOI.sess_raw{s1}.electro{1});
                    
                    dir_elfig = fullfile(dir_subj_res,'fig_el');
                    if ~exist(dir_elfig,'dir'), mkdir(dir_elfig); end
                    %if expedite
                    el = load(IOI.res.el{s1});
                    fname_el = fullfile(dir_elfig,[short_el_label '_' sess_label gen_num_str(s1,2)]);
                    ioi_plot_LFP(IOI,el.el,s1,1,1,fname_el);
                    disp(['Done processing session: ' int2str(s1) ', electrophysiology']);
                end
                %toc
                if ~expedite
                    disp(['Done processing session ' int2str(s1) ' images (' int2str(n_frames) 'images)']);
                end
            end
        end %for s1
    end
catch exception
    disp(exception.identifier)
    disp(exception.stack(1))
end