function out = ioi_msioi_run(job)
%Note that for the OD data, time is saved as the z dimension here, which may be
%more convenient for viewing. Later, in ioi_concentrations_run and
%ioi_flow_run, time will be saved as the 4th dimension.
% TODO:
%   -Calculate median over whole session instead of over block
%   -Hard coded the 6 frame steps for different colors, not sure it is
%   always true
%   -Image shrinkage should be implemented by interpolation, not by
%   choosing one point out of two
%   -Need to implement memory management better, follow examples from SPM
%
%   -if memmapfile works well, consider removing save_choice==2 && ==3
%   entirely, remove all loops over f1 in all modules, and simplify image
%   interpolation (bring back Fred's code)

%IOI structure: one per subject, good for multi-sessions
%IOI.dir: various directories
%IOI.sess_raw and IOI.sess_res: session specific raw and processed information
%IOI.res: results
%IOI.dev: device-level recording settings
%block_size = 100; %number of images per block, if save_choice = 2 selected

out.IOImat = {};
%Color names
str_anat = 'G'; %Varying number of colors allowed
str_laser = 'L';
str_red = 'R'; 
str_green = 'G';
str_yellow = 'Y';
str_contrast = 'C'; %Optional
str_color = [str_red str_green str_yellow str_laser str_contrast]; %red, green, yellow and laser speckle
str_color_french = [str_red 'V' 'J' str_laser str_contrast];

%minimum number of image files per session
sess_min_image_files = job.sess_min_image_files;
%               Directory Structure
%                   dir_group_all
%             /                     \
%     dir_group_raw              dir_group_res
%           |                       |
%     dir_subj_raw               dir_subj_res
%
dir_group_res_default = 'Res';
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
    save_choice = 1; %other choices no longer maintained, superseded by
end
forceProcessingOn = job.forceProcessingOn;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Big loop over subjects
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for SubjIdx=1:length(job.subj.top_bin_dir)
    try
        tic
        clear IOI sess_raw sess_res
        subj_OK = 1; %Boolean to exit loop if subject data is too corrupted
        %Directory structure
        dir_subj_raw = job.subj.top_bin_dir{SubjIdx}; %location of raw data for this subject
        tmpsep = strfind(dir_subj_raw,filesep);
        dir_group_raw = dir_subj_raw(1:tmpsep(end-1)); %group level dir of raw data
        dir_group_all = dir_subj_raw(1:tmpsep(end-2)); %group level dir for all data
        subj_name = dir_subj_raw((tmpsep(end-1)+1):(tmpsep(end)-1)); %name of subject
        %path configuration
        if isfield(job.output_path_choice, 'output_path_select')
            dir_group_res = job.output_path_choice.output_path_select.output_path; %understood as group output path
            %dir_group_default = 0;
        else
            %dir_group_default = 1;
            dir_group_res = fullfile(dir_group_all,dir_group_res_default); %group level dir of processed data
        end
        dir_subj_res = fullfile(dir_group_res,subj_name);
        if (~exist(fullfile(dir_subj_res,'IOI.mat'),'file') || job.subj.force_redo)
            %Use existing IOI.mat in case force_redo==2 was chosen
            if (exist(fullfile(dir_subj_res,'IOI.mat'),'file') && job.subj.force_redo==2)
                load(fullfile(dir_subj_res,'IOI.mat'));
                PartialRedo2 = 1;
            else
                PartialRedo2 = 0;
            end
            if ~PartialRedo2
                %fill IOI structure
                IOI.warning = {};
                IOI.save_choice = save_choice;
                switch save_choice
                    case 1
                        IOI.save_choice_method = 'One file per session';
                    case 2
                        IOI.save_choice_method = 'One file per block';
                    case 3
                        IOI.save_choice_method = 'One file per image (frame)';
                end
                IOI.subj_name = subj_name;
                IOI.dir.dir_group_all = dir_group_all;
                IOI.dir.dir_group_raw = dir_group_raw;
                IOI.dir.dir_group_res = dir_group_res;
                IOI.dir.dir_subj_raw = dir_subj_raw;
                IOI.dir.dir_subj_res = dir_subj_res;
                IOI.color.eng = str_color;
                IOI.color.french = str_color_french;
                IOI.color.anat = str_anat;
                IOI.color.laser = str_laser;
                IOI.color.red = str_red;
                IOI.color.green = str_green;
                IOI.color.yellow = str_yellow;
                IOI.color.contrasts = str_contrast;
                %Create required directories
                if ~exist(dir_group_res,'dir'),mkdir(dir_group_res); end
                if ~exist(dir_subj_res,'dir'),mkdir(dir_subj_res); end
            end
            %resume files - subject level
            [files_txt,dummy] = cfg_getfile('FPList',dir_subj_raw,'.txt');
            if length(files_txt) == 2
                %careful, the order of the files here could change if their
                %names were to be changed
                IOI.info.expt = files_txt{1,:}; %resume_manip
                IOI.info.comments = files_txt{2,:};%resumecommentaire
            else
                disp(['Wrong number of text files at subject level in ' dir_subj_raw]);
            end
            %collect files to analyze
            %generate output folders
            [dummy,dirs] = cfg_getfile('FPListRec',dir_subj_raw,'.bin');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %1- check that raw directory structure is consistent
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if PartialRedo2
                sess_raw = IOI.sess_raw;
            else
                sess_raw = {};
                for i=1:(length(dirs)/2)
                    sess = []; sess_OK = 1;
                    %skip sessions not selected by user
                    if all_sessions || sum(i==selected_sessions)
                        if ~strcmp([dirs{2*i-1,:} '_images'],dirs{2*i,:})
                            disp(['Missing directory: ' dirs{2*i,:}]);
                        else
                            %check that text files are present
                            [files_txt,dummy] = cfg_getfile('FPListRec',dirs{2*i-1,:},'.txt');
                            if length(files_txt) == 3
                                sess.Frame = files_txt{1,:};
                                sess.TTL = files_txt{2,:}; %TTLseuil
                                sess.info = files_txt{3,:};
                                [files_bin,dummy] = cfg_getfile('FPListRec',dirs{2*i,:},'.bin');
                                sess.images = files_bin; %all images
                                for c1=1:length(str_color)
                                    %str1 = str_color(c1);
                                    strf = str_color_french(c1);
                                    %separate images into various colors
                                    [files_bin,dummy] = cfg_getfile('FPListRec',dirs{2*i,:},[strf strf]);
                                    sess.fnames{c1} = files_bin;
                                    if c1 ==1
                                        image_count = length(files_bin);
                                        if image_count < sess_min_image_files
                                            sess_OK = 0;
                                            disp(['Insufficient number of images for ' dirs{2*i-1,:} ' ...skipping.']);
                                        end
                                    else
                                        %weaker criterion - number of image
                                        %files might differ by 1 because
                                        %one only one extra image
                                        if abs(length(files_bin)- image_count)>= 1
                                            if ~strcmp(strf,str_contrast) %OK if contrast files missing
                                                sess_OK = 0;
                                                disp(['Problem with number of ' strf strf ' images for ' dirs{2*i-1,:}]);
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
                                    disp(['Problem with electro file for ' dirs{2*i-1,:}]);
                                end
                                %add found session to session list
                                if forceProcessingOn
                                    sess_OK = 1; %ignore warnings, and include session anyway
                                end
                                if sess_OK, sess_raw = [sess_raw sess]; end
                            else
                                disp(['Problem with text files in ' dirs{2*i,:}]);
                            end
                        end
                    end
                end
                if ~isempty(sess_raw)
                    IOI.sess_raw = sess_raw;
                else
                    disp('No sessions found... skipping this subject');
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
                        [scan_info, physio]=ioi_extract_info(IOI.info.expt,...
                            sess_raw{s1}.info,sess_raw{s1}.TTL,sess_raw{s1}.Frame);
                        if s1==1
                            % Identifiy times for each frame and each color separately
                            IOI.dev.acq_freq_hz=(scan_info.Frame(end,1)-scan_info.Frame(1,1))/scan_info.Frame(end,3);
                            % Fix TR for everyone, the factor 6 comes from the fact that we need 6
                            % camera frames to return to the same colors
                            IOI.dev.TR=6/IOI.dev.acq_freq_hz;
                        end
                        % 6 images make one frame, each frame is indexed, partial frames are
                        % counted as one frame - this is only approximate
                        n_frames = ceil((scan_info.Frame(end,1)-scan_info.Frame(1,1)+1)/6);
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
                            stim.onset_frame=ceil(stim.onset_frame/6);
                            %Next 4 entries: Approximate info - do not use:
                            stim.intensity=scan_info.stim3(stim_index,2);
                            stim.train_period=scan_info.stim3(stim_index,3);
                            stim.n_stim_in_train=scan_info.stim3(stim_index,4);
                            stim.train_duration=scan_info.stim3(stim_index,5);
                            %IOI.dev.stim{stim_index}=stim;
                            list_stim = [list_stim stim];
                        end
                        clear names onsets durations
                        % Everything is in scan number here
                        for stim_index=1:length(list_stim)
                            names{stim_index}=['Stim_',num2str(stim_index)];
                            onsets{stim_index}=list_stim{stim_index}.onset_frame';
                            durations{stim_index}=round(list_stim{stim_index}.train_duration/IOI.dev.TR);
                        end
                        sess.list_stim = list_stim;
                        sess.names = names; %not quite the SPM format, need to put in Sess.U format
                        sess.onsets = onsets;
                        sess.durations = durations;
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
                
                %leave voxel size in arbitrary units for now
                vx = [1 1 1];
                % Do and save anatomy file
                %use first green colored image as the anatomy
                image_total=private_read_single_bin_image(sess_raw{1}.fnames{str_anat == str_color}{1,:});
                fname = fullfile(dir_subj_res, [subj_name '_anat.nii']);
                ioi_save_nifti(image_total, fname, vx);
                IOI.res.file_anat=fname;
                
                % For each file to treat, reassemble and save to nifti
                % We save optical density as this is more convenient for further processing.
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %4- Functional images and electrophysiology
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %loop over sessions
                for s1=1:length(sess_raw)
                    if ~PartialRedo2 || all_sessions || sum(s1==selected_sessions)
                        %tic
                        %create dir for this session
                        if s1<10, str_s1 = '0'; else str_s1 = ''; end
                        str_s1 = [str_s1 int2str(s1)];
                        sess_str = ['Sess' str_s1];
                        dir_sess_res = fullfile(dir_subj_res,sess_str);
                        if ~exist(dir_sess_res,'dir'),mkdir(dir_sess_res); end
                        n_frames = sess_res{s1}.n_frames;
                        hasRGY = [];
                        %loop over colors
                        for c1=1:length(str_color)
                            need_resize_message = 1;
                            str1 = str_color(c1);
                            if ~(str1==str_contrast && isempty(sess_raw{s1}.fnames{c1})) %exclude when contrast images are absent
                            fcount=0; %image counter (actual frames recorded)
                            im_count = 0; %image counter (all images including interpolated ones)
                            %loop over image files for this session
                            image_total = [];
                            %tic
                            %Do NOT shrink laser speckle!
                            if ~(str1 == str_laser)
                                if shrinkageOn
                                    %if ~PartialRedo2
                                    %With PartialRedo2 == 1, keep old choice
                                    %of shrink_x -- this makes this option
                                    %useful if images for different
                                    %sessions were saved with different
                                    %resolutions.
                                    IOI.res.shrink_x=shrink_x;
                                    IOI.res.shrink_y=shrink_y;
                                    %end
                                    vx=[shrink_x shrink_y 1];
                                end
                            else
                                vx = [1 1 1];
                            end
                            
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
                                %image_part is of shape [nx, ny, time] and has been
                                %interpolated for missing frames
                                if (str1 == str_contrast)
                                    tc1 = find(str_color==str_laser);
                                else
                                    tc1 = c1;
                                end
                                [image_part fcount]=private_read_bin_image(sess_raw{s1}.fnames{c1}{f1,:}, ...
                                    fcount,ceil(sess_res{s1}.camera{tc1}/6), vx,n_frames, first_file, last_file,im_count);
                                image_part = single(image_part);
                                [nx ny nt] = size(image_part);
                                im_count = im_count+nt;
                                if ~(str1 == str_laser || str1 == str_contrast) %not laser speckle
                                    %calculate median for this block (about 80 images) along time direction
                                    if ~(save_choice==1)
                                        image_median = median(image_part,3);
                                        image_part = -log(image_part./repmat(image_median,[1 1 nt]));
                                    end
                                end
                                if f1==1 && memmapfileOn
                                    %create temporary file
                                    fid = fopen('tmp_file.dat','w');
                                    fwrite(fid, zeros([nx ny n_frames]), 'single');
                                    fclose(fid);
                                    im_obj = memmapfile('tmp_file.dat','format',...
                                        {'single' [nx ny n_frames] 'image_total'},...
                                        'writable',true);
                                end
                                %catch case when there is a resizing of the images during acquisition
                                if f1>1
                                    if ~(size(image_part,1)==size(image_total,1) && size(image_part,2)==size(image_total,2))
                                        if need_resize_message
                                            need_resize_message = 0; %first pass
                                            warning_message = ['Image resizing required for session ' ...
                                                int2str(s1) ' file ' int2str(f1) ' color ' int2str(c1)];
                                            IOI.warning = [IOI.warning warning_message];
                                            disp(warning_message);
                                        end
                                        image_part = ioi_imresize(image_part,1,size(image_total,1),size(image_total,2),1,1);
                                    end
                                end
                                switch save_choice
                                    case 1
                                        if ~memmapfileOn
                                            image_total = cat(3,image_total, image_part);
                                        else
                                            im_obj.Data.image_total(:,:,im_count-nt+1:im_count) = image_part;
                                        end
                                    case 2
                                        %save block
                                        first_image_str = gen_num_str(im_count-nt+1);
                                        last_image_str = gen_num_str(im_count);
                                        fname = fullfile(dir_subj_res,sess_str, [subj_name '_OD_' str1 '_S' str_s1 '_' first_image_str '_to_' last_image_str '.nii']);
                                        IOI.sess_res{s1}.fname{c1} = [IOI.sess_res{s1}.fname{c1}; fname];
                                        ioi_save_nifti(image_part, fname, vx);
                                    case 3
                                        %save each frame of block separately
                                        for j1=1:nt
                                            image_str = gen_num_str(im_count-nt+j1);
                                            fname = fullfile(dir_subj_res,sess_str, [subj_name '_OD_' str1 '_S' str_s1 '_' image_str '.nii']);
                                            IOI.sess_res{s1}.fname{c1} = [IOI.sess_res{s1}.fname{c1}; fname];
                                            ioi_save_nifti(squeeze(image_part(:,:,j1)), fname, vx);
                                        end
                                end
                            end
                            if ~(im_count==n_frames)
                                disp('Problem with image interpolation - wrong final number');
                                disp(['n_frames = ' int2str(n_frames) '; final number images = ' int2str(im_count)]);
                                IOI.warning = [IOI.warning ['Session ' int2str(s1) ',color ' str1 ...
                                    'n_frames = ' int2str(n_frames) '; final number images = ' int2str(im_count)]];
                                
                            end
                            if save_choice == 1
                                if ~(str1 == str_laser || str1 == str_contrast) %not laser speckle
                                    %calculate median for whole image time direction
                                    if ~memmapfileOn
                                        image_median = median(image_total,3);
                                        image_total = -log(image_total./repmat(image_median,[1 1 size(image_total,3)]));
                                    else
                                        image_median = median(im_obj.Data.image_total,3);
                                        for i3=1:n_frames
                                            im_obj.Data.image_total(:,:,i3) = -log(im_obj.Data.image_total(:,:,i3)./image_median);
                                        end
                                    end
                                end
                                fname = fullfile(dir_subj_res,sess_str, [subj_name '_OD_' str1 '_S' str_s1 '.nii']);
                                IOI.sess_res{s1}.fname{c1} = [IOI.sess_res{s1}.fname{c1}; fname];
                                if ~memmapfileOn
                                    ioi_save_nifti(image_total, fname, vx);
                                else
                                    ioi_save_nifti(im_obj.Data.image_total, fname, vx);
                                    clear im_obj;
                                    delete('tmp_file.dat');
                                end
                            end
                            if ~(str1 == str_laser || str1 == str_contrast)
                                hasRGY = [hasRGY str1];
                            end
                            disp(['Done processing session: ' int2str(s1) ', color: ' str1]);
                            %toc
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
            % Write out IOI in .mat file format at the subject results path location
            % All (links to) data are preserved in this structure
            IOImat = fullfile(dir_subj_res,'IOI.mat');
            save(IOImat,'IOI');
        else
            %IOI.mat already exists, and rerun is not enforced - therefore just
            %skip this module and pass IOImat structure to the next module
            IOImat = fullfile(dir_subj_res,'IOI.mat');
        end
        
        out.IOImat = [out.IOImat IOImat];
        toc
        disp(['Subject ' int2str(SubjIdx) ' complete']);
    catch exception
        disp(exception.identifier)
        disp(exception.stack(1))
    end
end
%END

function str_out = gen_num_str(num_in)
%utility to pad with zeros in front, so that files are listed in
%consecutive image frame order in Windows Explorer
%string length desired
str_len = 5;
if num_in > 0
    num_dig = floor(log10(num_in));
else
    num_dig = -1;
end
str = repmat('0',[1 str_len-num_dig-1]);
str_out = [str int2str(num_in)];

% Internal private function to read binary images, facilitates reading the
% code above.
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
image_total=reshape(temp((5+sz_add):end,:),[s1 s2 last_index]);
if vx(1) > 1 || vx(2) > 1
    image_total = ioi_imresize(image_total,0,nx,ny,vx(1),vx(2));
end
image_total = single(image_total);

% Internal private function to read binary images, facilitates reading the
% code above.
function image_total = private_read_single_bin_image(fname)
% Open file containing a single volume and return it
%
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
% To correct, this should be an interpolation instead of of simple
% We add a single z dimension to be compatible with nifti
if ~isempty(noImageNouveauFormat)  % new format
    %image_total=reshape(temp(8:end,:),[nx ny 1 1]);
    image_total=reshape(temp(8:end,:),[nx ny]);
else
    %image_total=reshape(temp(5:end,1),[nx ny 1 1]);
    image_total=reshape(temp(5:end,1),[nx ny]);
end
image_total = single(image_total);