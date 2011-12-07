function IOI = ioi_msread_old_format(IOI,job)
try
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
    %leave voxel size in arbitrary units for now, for anatomical image
    vx_anat = [1 1 1];
    %Variable changed meaning here - more convenient for user to specify
    %duration in seconds
    min_session_duration = job.sess_min_image_files;
    temp_TR = 0.2;
    temp_ImNum = 80;
    nzero_padding = 5;
    
    %shrinkage configuration
    if isfield(job.configuration_choice,'configuration_shrink')
        shrinkageOn = 1;
        shrink_x = job.configuration_choice.configuration_shrink.shrink_x;
        shrink_y = job.configuration_choice.configuration_shrink.shrink_y;
    else
        shrinkageOn = 0;
    end
    %select a subset of sessions
    if isfield(job.session_choice,'select_sessions')
        all_sessions = 0;
        selected_sessions = job.session_choice.select_sessions.selected_sessions;
    else
        all_sessions = 1;
    end
    %choose saving mode
    memmapfileOn = job.memmapfileOn;
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
        if ~isempty(sess_raw)
            IOI.sess_raw = sess_raw;
        else
            disp(['No sessions found... skipping subject ' int2str(SubjIdx) ' (namely: ' job.top_bin_dir{SubjIdx} ')']);
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
                    IOI.dev.acq_freq_hz=(scan_info.Frame(end,1)-scan_info.Frame(1,1))/scan_info.Frame(end,3);
                    % Fix TR for everyone, the factor 6 comes from the fact that we need 6
                    % camera frames to return to the same colors
                    IOI.dev.TR=frames_per_cycle/IOI.dev.acq_freq_hz;
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
                clear names onsets durations
                %Converted to seconds, rather than frame number
                for stim_index=1:length(list_stim)
                    names{stim_index}=['Stim_',num2str(stim_index)];
                    onsets{stim_index}=list_stim{stim_index}.onset_frame'*IOI.dev.TR;
                    durations{stim_index}=list_stim{stim_index}.train_duration;
                end
                %Store onset information
                sess.list_stim = list_stim;
                sess.names = names; %not quite the SPM format, need to put in Sess.U format
                sess.onsets = onsets; %in seconds
                sess.durations = durations; %in seconds
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
        IOI.res.shrinkageOn = shrinkageOn;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %3- Anatomical image
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %use first green colored image as the anatomy
        image_total=private_read_single_bin_image(sess_raw{1}.fnames{str_anat == str_color}{1,:});
        fname = fullfile(dir_subj_res, [subj_name '_' suffix_for_anat_file '.nii']);
        ioi_save_nifti(image_total, fname, vx_anat);
        IOI.res.file_anat=fname;
        
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
                            if shrinkageOn
                                %Keep same whatever value of PartialRedo2
                                IOI.res.shrink_x=shrink_x;
                                IOI.res.shrink_y=shrink_y;
                                vx=[shrink_x shrink_y 1];
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
                            [image_part fcount]=private_read_bin_image(sess_raw{s1}.fnames{c1}{f1,:}, ...
                                fcount,ceil(sess_res{s1}.camera{tc1}/frames_per_cycle),vx,n_frames,first_file,last_file,im_count);
                            %image_part is of shape [nx, ny, 1, time] and has been
                            %interpolated for missing frames
                            [nx ny nz nt] = size(image_part); %Here nz = 1, always
                            im_count = im_count+nt;
                            if f1==1
                                nx0 = nx;
                                ny0 = ny;
                            end
                            if f1==1 && memmapfileOn
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
                            if ~memmapfileOn
                                image_total = cat(4,image_total,image_part);
                            else
                                im_obj.Data.image_total(:,:,:,im_count-nt+1:im_count) = image_part;
                            end
                        end %end of loop over binary files
                        if ~(im_count==n_frames)
                            disp('Problem with image interpolation - wrong final number');
                            IOI = disp_msg(IOI,['Session ' int2str(s1) ',color ' str1 ...
                                'n_frames = ' int2str(n_frames) '; final number images = ' int2str(im_count)]);
                        end
                        %calculate median for whole image in time direction
                        if ~(str1 == str_laser || str1 == str_contrast) %not laser speckle nor contrast
                            if ~memmapfileOn
                                image_median = median(image_total,4);
                                image_total = -log(image_total./repmat(image_median,[1 1 1 n_frames]));
                            else
                                image_median = median(im_obj.Data.image_total,4);
                                im_obj.Data.image_total = -log(im_obj.Data.image_total./repmat(image_median,[1 1 1 n_frames]));
                            end
                        end
                        %Save the median
                        IOI.sess_res{s1}.fname_median{c1} = fullfile(dir_subj_res,sess_str, ...
                            [subj_name '_' OD_label '_median_' str1 '_' sess_str '.nii']);
                        ioi_save_nifti(single(image_median),IOI.sess_res{s1}.fname_median{c1},vx);
                        try
                            if ~memmapfileOn
                                %min and max and relative change
                                min_image = min(image_total,[],4);
                                max_image = max(image_total,[],4);
                                %10th and 90th percentiles and relative change
                                tenthpctle_image = prctile(image_total,10,4);
                                ninetiethpctle_image = prctile(image_total,90,4);
                            else
                                min_image = min(im_obj.Data.image_total,[],4);
                                max_image = max(im_obj.Data.image_total,[],4);
                                tenthpctle_image = prctile(im_obj.Data.image_total,10,4);
                                ninetiethpctle_image = prctile(im_obj.Data.image_total,90,4);
                            end
                            change = max_image ./min_image;
                            change_90_10 = ninetiethpctle_image ./ tenthpctle_image;
                            sess.fname_min{c1} = fullfile(dir_subj_res,sess_str, ...
                                [subj_name '_' OD_label '_min_' str1 '_' sess_str '.nii']);
                            sess.fname_max{c1} = fullfile(dir_subj_res,sess_str, ...
                                [subj_name '_' OD_label '_max_' str1 '_' sess_str '.nii']);
                            sess.fname_10pctle{c1} = fullfile(dir_subj_res,sess_str, ...
                                [subj_name '_' OD_label '_10pctle_' str1 '_' sess_str '.nii']);
                            sess.fname_90pctle{c1} = fullfile(dir_subj_res,sess_str, ...
                                [subj_name '_' OD_label '_90pctle_' str1 '_' sess_str '.nii']);
                            sess.fname_change{c1} = fullfile(dir_subj_res,sess_str, ...
                                [subj_name '_' OD_label '_change_' str1 '_' sess_str '.nii']);
                            sess.fname_change_90_10{c1} = fullfile(dir_subj_res,sess_str, ...
                                [subj_name '_' OD_label '_change_90_10_' str1 '_' sess_str '.nii']);
                            ioi_save_nifti(single(min_image),sess.fname_min{c1},vx);
                            ioi_save_nifti(single(max_image),sess.fname_max{c1},vx);
                            ioi_save_nifti(single(change),sess.fname_change{c1},vx);
                            ioi_save_nifti(single(tenthpctle_image),sess.fname_10pctle{c1},vx);
                            ioi_save_nifti(single(ninetiethpctle_image),sess.fname_90pctle{c1},vx);
                            ioi_save_nifti(single(change_90_10),sess.fname_change_90_10{c1},vx);
                        end
                        %vx currently not used in ioi_save_nifti
                        %and ioi_write_nifti, despite apparent dependency
                        switch save_choice
                            case 1
                                fname = fullfile(dir_subj_res,sess_str, [subj_name '_' OD_label '_' str1 '_S' str_s1 '.nii']);
                                IOI.sess_res{s1}.fname{c1} = [IOI.sess_res{s1}.fname{c1}; fname];
                                if ~memmapfileOn
                                    ioi_save_nifti(image_total,fname,vx);
                                else
                                    ioi_save_nifti(im_obj.Data.image_total,fname,vx);
                                end
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
                                    if ~memmapfileOn
                                        ioi_save_nifti(image_total(:,:,:,j1:last_index), fname, vx);
                                    else
                                        ioi_save_nifti(im_obj.Data.image_total(:,:,:,j1:last_index), fname, vx);
                                    end
                                end
                            case 3
                                %save each frame separately
                                for j1=1:n_frames
                                    image_str = gen_num_str(j1,nzero_padding);
                                    fname = fullfile(dir_subj_res,sess_str, [subj_name '_' OD_label '_' str1 '_S' str_s1 '_' image_str '.nii']);
                                    IOI.sess_res{s1}.fname{c1} = [IOI.sess_res{s1}.fname{c1}; fname];
                                    if ~memmapfileOn
                                        ioi_save_nifti(image_total(:,:,:,j1),fname,vx);
                                    else
                                        ioi_save_nifti(im_obj.Data.image_total(:,:,:,j1),fname,vx);
                                    end
                                end
                        end
                        if memmapfileOn %clean up
                            clear im_obj;
                            delete('tmp_file.dat');
                        end
                        %Add found color - used later in concentration calculation module
                        if ~(str1 == str_laser || str1 == str_contrast)
                            hasRGY = [hasRGY str1];
                        end
                        disp(['Done processing session: ' int2str(s1) ', color: ' str1]);
                    end
                end %for c1
                IOI.sess_res{s1}.hasRGY = hasRGY;
                %electrophysiology
                if isfield(IOI.sess_raw{s1},'electroOK')
                    IOI = extract_electro(IOI,s1,IOI.sess_raw{s1}.electro{1});
                    IOI.res.electroOK = 1; %a bit early, but should be OK
                    disp(['Done processing session: ' int2str(s1) ', electrophysiology']);
                end
                %toc
                disp(['Done processing session ' int2str(s1) ' images (' int2str(n_frames) 'images)']);
            end
        end %for s1
    end
catch exception
    disp(exception.identifier)
    disp(exception.stack(1))
end



% Internal private function to read binary images
function [image_total fcount_out]= private_read_bin_image(fnames,...
    fcount,indices,vx,n_frames,first_file,last_file,im_count)
%If frames are missing between this fnames file and the next, by a convention
%defined here, they will be added at the end of this block of frames (not the next).
%fcount and fcount_out are the number of recorded frames before and after
%processing of the current block
read_format = 'int16';
% Open first file to see which format is used (new or old)
% Describe image formats here to make code clearer
%
fidA = fopen(fnames);
I = fread(fidA,2,read_format);
if all(I==[1;0]) % New format including image number
    noImageNouveauFormat=fread(fidA,1,read_format);
    I = fread(fidA,4,read_format);
else
    noImageNouveauFormat=[];
    I =[I; fread(fidA,2,read_format)];
end
% Image dimensions read from header (first few int16)
s1= I(3);
s2= I(1);
fclose(fidA);

% Now read all files with the dimensions determined above, variable
% image_total will contain the reconstituted images
nx=size(vx(1):vx(1):s1,2);
ny=size(vx(2):vx(2):s2,2);
fidA = fopen(fnames);
image_part=int16((fread(fidA,read_format)));
fclose(fidA);
% Reshape to images
if ~isempty(noImageNouveauFormat) % new format
    sz_add = 3;
else
    sz_add = 0;
end
image_part=reshape(image_part,(s1*s2+sz_add+4),round(length(image_part)/(s1*s2+sz_add+4)));
list_frames = fcount+(1:size(image_part,2)); %this list does not know about missing frames
fcount_out = list_frames(end);
last_missing_is_last_index = 0;
%first_missing_is_first_index = 0;
diff0 = diff(indices)~=1;
%missing_indices_locations=find(diff0); %location of missing frames in indices list
missing_indices=indices(diff0); %actual first missing frame numbers

%catch in case we are at the last list of frames
try
    missing_indices_select = missing_indices < indices(list_frames(end)+1) & ...
        missing_indices >= indices(list_frames(1)); %?
catch
    missing_indices_select = missing_indices <= indices(list_frames(end)) & ...
        missing_indices >= indices(list_frames(1));
    if indices(end) < n_frames
        indices(end+1) = n_frames+1;
        %redo previous steps
        diff0 = diff(indices)~=1;
        %missing_indices_locations=find(diff0); %location of missing frames in indices list
        missing_indices=indices(diff0);
        missing_indices_select = missing_indices < indices(list_frames(end)+1) & ...
            missing_indices >= indices(list_frames(1));
    end
end
diff1 = diff(indices);
diff_index = diff1(diff0); %gap size
diff_index = diff_index(missing_indices_select);
missing_indices = missing_indices(missing_indices_select);
%missing_indices_locations = missing_indices_locations(missing_indices_select);

if ~isempty(missing_indices) && (missing_indices(end) >= indices(list_frames(end))) %actual frames
    last_index = indices(list_frames(end)+1)-indices(list_frames(1));
    last_missing_is_last_index = 1;
else
    last_index = indices(list_frames(end))-indices(list_frames(1))+1;
end

temp = zeros(size(image_part,1),last_index);
temp(:,indices(list_frames)-indices(list_frames(1))+1) = image_part(:,list_frames-fcount);
% At each location fill the gaps naively
shft = indices(list_frames(1))-1;
for idx=1:length(missing_indices)
    b=missing_indices(idx)+1;
    e=b+diff_index(idx)-2;
    if ~last_missing_is_last_index || ~(idx == length(missing_indices))
        % Find gap width
        for n=b:e
            %left side: actual frames (with first index reset to 1)
            %right side:
            temp(:,n-shft)=0.5*(temp(:,b-1-shft)+temp(:,e+1-shft));
        end
    else
        %in case last_missing_is_last_index, then fill all missing with
        %last available image in this block.
        for n=b:e
            temp(:,n-shft)=temp(:,b-1-shft);
        end
    end
end
if first_file == 1
    %check whether first image is missing
    first_rec_frame = indices(list_frames(1));
    if first_rec_frame > 1
        %augment temp and last_index
        last_index = last_index+first_rec_frame-1;
        temp = [repmat(temp(:,1),[1 first_rec_frame-1])  temp];
    end
end
output_debugging_info = 0;
if output_debugging_info
    %Debugging info
    format compact
    disp(['Fr processed earlier: ' int2str(fcount) ', Frames now: ' int2str(fcount_out)]);
    disp(['Fr available: ' int2str(fcount_out-fcount) ', Images added: ' int2str(last_index)]);
    disp('missing_indices_locations:');
    disp(missing_indices_locations');
    disp('missing_indices:');
    disp(missing_indices');
    disp([int2str(list_frames(1)) ': list_frames(1), ' int2str(indices(list_frames(1))) ': indices(list_frames(1))']);
    disp([int2str(list_frames(end)) ': list_frames(end), ' int2str(indices(list_frames(end))) ': indices(list_frames(end))']);
    disp('**********');
end
%in the hopefully rare instances when there is still a mismatch between
%n_frames and the total number of interpolated images, enforce the n_frames
%value by adding or subtracting
if last_file
    im_diff = n_frames-im_count-last_index;
    if im_diff>0
        %missing some frames, add some
        temp = [temp repmat(temp(:,end),[1 im_diff])];
        last_index = last_index + im_diff;
        disp(['Warning: missing ' int2str(im_diff) ' images at the end, filled in to get n_frames']);
    else
        if im_diff<0
            %too many frames, drop some
            temp = temp(:,1:(end+im_diff));
            last_index = last_index + im_diff;
            disp(['Warning (minor): too many images ' int2str(abs(im_diff)) ' at the end, removed some to get n_frames']);
        end
    end
end
%Put here in 4D format
image_total=reshape(temp((5+sz_add):end,:),[s1 s2 1 last_index]);
if vx(1) > 1 || vx(2) > 1
    image_total = ioi_imresize(image_total,0,nx,ny,vx(1),vx(2));
end
image_total = single(image_total);

% Internal private function to read anatomical binary image
function image_total = private_read_single_bin_image(fname)
% Open file containing a single volume
read_format = 'int16';
fidA = fopen(fname);
I = fread(fidA,2,read_format);
if all(I==[1;0]) % New format including image number
    noImageNouveauFormat=fread(fidA,1,read_format);
    I = fread(fidA,4,read_format);
else
    noImageNouveauFormat=[];
    I =[I; fread(fidA,2,read_format)];
end
% Image dimensions read from header (first few int16)
nx= I(3);
ny= I(1);
fclose(fidA);

fidA = fopen(fname);
image_total=int16((fread(fidA,read_format)));
fclose(fidA);
% Reshape to images
if ~isempty(noImageNouveauFormat) % new format
    temp(:,:)=squeeze(reshape(image_total,(nx*ny+3+4),round(length(image_total)/(nx*ny+3+4))));
else
    temp(:,:)=squeeze(reshape(image_total,(nx*ny+4),round(length(image_total)/(nx*ny+4))));
end
clear image_total
if ~isempty(noImageNouveauFormat)  % new format
    image_total=reshape(temp(8:end,:),[nx ny]);
else
    image_total=reshape(temp(5:end,1),[nx ny]);
end
image_total = single(image_total);

