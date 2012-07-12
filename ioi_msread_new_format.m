function IOI = ioi_msread_new_format(IOI,job)
try
    image_anat_first_pass = 1;
    %select a subset of sessions
    [all_sessions selected_sessions] = ioi_get_sessions(job);

    if isfield(job.treatment_mode,'expeditive_mode')
        expedite = 1;
    else
        expedite = 0;
    end
    %recall dir structure
    subj_name = IOI.subj_name;
    dir_subj_raw = IOI.dir.dir_subj_raw;
    dir_subj_res = IOI.dir.dir_subj_res;
    
    %Color names - could be put in user interface?
    str_anat = 'G';
    str_laser = 'L'; %Optional
    str_red = 'R'; %Varying number of colors allowed
    str_green = 'G';
    str_yellow = 'Y';
    str_green_french = 'V'; %vert
    str_yellow_french = 'J'; %jaune
    %Code will read info.txt file for each session to find actual color
    %order, and then will swap the order to get str_color order
    str_color_default = [str_green str_yellow str_laser str_red]; %VJLR -- this is the default at time of recording
    %If adding or removing colors, make sure the string of English and French
    %colors have the same length and are a one-to-one mapping (code loops over colors
    %to find images):
    str_color = [str_red str_green str_yellow str_laser]; %red, green, yellow and laser speckle -- this is the desired color order for processed images
    hasRGY = [str_red str_green str_yellow];
    nColors = length(str_color);
    OD_label = 'OD'; %label for optical density images
    suffix_for_anat_file = 'anat'; %to build anatomical image name
    sess_label = 'S'; %prefix for name of directories for each session
    el_label = 'electro'; %only use for raw electro file
    el_info = 'info';
    short_el_label = 'el'; %short name for output electro file
    el_suffix = 'tdms';
    nzero_padding = 5;
    %leave voxel size in arbitrary units for now, for anatomical image
    vx_anat = [1 1 1];
    vx = [1 1 1];
    % Is it not dangerous to let the user choose the acq_freq? //EGC
    if isfield(job,'acq_freq')
        TR = job.color_number/job.acq_freq;
    else
        TR = 0.2;
    end
    if isfield(job,'stim_cutoff')
        stim_cutoff = job.stim_cutoff;
    else
        stim_cutoff = 5;
    end
    %subj_OK = 1; %Boolean to exit loop if subject data is too corrupted
    if ~job.PartialRedo2
        IOI.color.eng = str_color;
        IOI.color.red = str_red;
        IOI.color.green = str_green;
        IOI.color.yellow = str_yellow;
        IOI.color.laser = str_laser;
        IOI.res.shrinkageOn = 0; %no shrinkage of images
        IOI.dev.TR = TR;
    end
    
    if ~job.PartialRedo2
        [dummy,dirs] = cfg_getfile('FPListRec',dir_subj_raw,'.bin');
        sC = 0; %good session counter
        %loop over sessions
        for s1=1:length(dirs)
            if all_sessions || sum(s1==selected_sessions)
                clear sess
                sC = sC + 1;
                %session directory -- relabel
                sess_str = [sess_label gen_num_str(s1,2)];
                dir_sess = fullfile(dir_subj_res,sess_str);
                if ~exist(dir_sess,'dir'), mkdir(dir_sess); end
                %read info file
                try
                    [file_info,dummy] = cfg_getfile('FPList',dirs{s1},'info.txt');
                    fid = fopen(file_info{1},'r');
                    t0 = fread(fid,'uint8');
                    fclose(fid);
                    t0 = char(t0');
                    t0 = deblank(t0);
                    color_order_french = t0(end-3:end);
                    if strcmp(color_order_french,'flat')
                        color_order_french = 'RVJL';
                    end
                    %color_order_french = color_order_french(1:4);
                    color_order = color_order_french;
                    %get English colors
                    for i0 = 1:length(color_order_french)
                        if strcmp(color_order_french(i0),str_yellow_french)
                            color_order(i0) = str_yellow;
                        end
                        if strcmp(color_order_french(i0),str_green_french)
                            color_order(i0) = str_green;
                        end
                    end
                catch
                    color_order = str_color_default;
                    IOI = disp_msg(IOI,['Problem with colors in info.txt for session ' int2str(s1) ';  defaulting to ' str_color_default]);
                end
                %Read electro
                [files_el,dummy] = cfg_getfile('FPList',dirs{s1},[el_label '.' el_suffix]);
                [ConvertedData,ConvertVer,ChanNames,GroupNames,ci]=convertTDMS(1,files_el{1});
                fil_loc = fullfile(ConvertedData.FileFolder,[el_label '.mat']);
                %Get default stimulations -- what if there are several
                %types of stimulations???
                stims = ConvertedData.Data.MeasuredData(6).Data;
                stims_dt= ConvertedData.Data.MeasuredData(6).Property(3).Value;
                ons0 = stims_dt*find(stims>stim_cutoff);
                ons1 = find(diff(ons0)>=1);
                if ~isempty(ons1)
                    ons = ons0([1 1+ons1']);
                else
                    ons = []; %in case we have a rest session
                end
                clear names onsets durations
                %Converted to seconds, rather than frame number
                for stim_index=1:1
                    names{stim_index}=['Stim_',num2str(stim_index)];
                    onsets{stim_index}=ons;
                    durations{stim_index}=ones(1,length(ons));
                end
                %Store onset information
                %sess.list_stim = list_stim;
                sess.names = names; %not quite the SPM format, need to put in Sess.U format
                sess.onsets = onsets; %in seconds
                sess.durations = durations; %in seconds
                
                %copy and remove
                %IOI.sess_res{sC} = {};
                elfile_info = fullfile(dir_subj_res,[short_el_label el_info '_' sess_label gen_num_str(sC,2) '.mat']);
                elfile = fullfile(dir_subj_res,[short_el_label '_' sess_label gen_num_str(sC,2) '.mat']);
                IOI.res.el{sC} = elfile;
                %IOI.res.elinfo{sC} = elfile_info; %not used in later modules as did not exist in earlier reading module
                copyfile(fil_loc, elfile_info);
                delete(fil_loc);
                %save electrophysiological data
                el = ConvertedData.Data.MeasuredData(5).Data;
                save(elfile,'el');
                IOI.res.electroOK = 1; %not by session...
                IOI.res.sfel = 1/stims_dt;
                if expedite
                    ioi_plot_LFP(IOI,el,s1);
                end
                if ~expedite
                    %extract ttl
                    ttl=TDMS2ttl(ConvertedData);
                    %save(elfile,el(1,:)); %incorrect.
                    %loop over image files for this session
                    %[imfiles dummy] = cfg_getfile('FPList',dirs{s1},'image');
                    %read frame present
                    [Imout frameout frameReadOut fileNo] = LectureImageMultiFast(dirs{s1},'image',-1);
                    fileNo = sort(fileNo);
                    missing_frames = [];
                    n_frames = ceil(max(frameReadOut(:,1))/nColors);
                    sess.n_frames = n_frames;
                    for i0=1:(n_frames*nColors)
                        tmp = find(i0==frameReadOut(:,1), 1);
                        if isempty(tmp)
                            missing_frames = [missing_frames i0];
                        end
                    end
                    iC = 0; %counter for number of images so far
                    
                    for c1=1:nColors
                        sess.fname{c1} = {};
                    end
                    sess.si = {};
                    sess.ei = {};
                    %Get image size
                    frames = frameReadOut(frameReadOut(1,3)==fileNo(1),:);
                    images=LectureImageMultiFast(dirs{s1},'image',frames);
                    %Check if images have duplicated pixels -- an early bug at
                    %acquisition, corrected for new acquisitions
                    [nx ny0] = size(images);
                    nxt = round(nx/2); nyt = round(ny0/2); %center of image
                    %two tests should be enough... -- weaker test now
                    if (images(nxt,nyt) == images(nxt,nyt+1)) || (images(nxt,nyt+1) == images(nxt,nyt+2))
                        %(images(nxt,nyt) == images(nxt,nyt+1)) && (images(nxt+3,nyt+3) == images(nxt+3,nyt+4))
                        %if mod(ny,2) == 0
                        dupOn = 1;
                        ny = length(images(1,1:2:end));
                        %else %This happens sometimes, that's OK
                        %    disp('Odd number of pixels in y direction, yet pixels are duplicated!')
                        %end
                    else
                        dupOn = 0;
                        ny = ny0;
                    end
                    
                    %shrink if required
                    [shrinkage_choice SH] = ioi_get_shrinkage_choice(job); 
                    IOI.res.shrinkageOn = shrinkage_choice; 
                    IOI.res.shrink_x = SH.shrink_x;
                    IOI.res.shrink_y = SH.shrink_y;                    
                    NX = nx;
                    NY = ny;
                    if shrinkage_choice
                        nx = floor(NX/SH.shrink_x);
                        if dupOn
                            ny = floor(NY/(SH.shrink_y/2));
                        else
                            ny = floor(NY/SH.shrink_y);
                        end
                        K.radius = SH.shrink_x;
                        K.k1 = NX;
                        K.k2 = NY;
                        K = ioi_spatial_LPF('set', K);
                    end
                    
                    %create large memmapfile for Y, R, G images (several GB)
                    fmem_name = 'all_images.dat';
                    %                 if exist(fmem_name,'file')
                    %                     try
                    %                         delete(fmem_name);
                    %                     catch
                    %                         disp([fmem_name ' : File could not be closed']);
                    %                     end
                    %                 end
                    fid = fopen(fmem_name,'w');
                    fwrite(fid, zeros([nx ny 1 n_frames nColors]), 'int16');
                    fclose(fid);
                    im_obj = memmapfile(fmem_name,'format',...
                        {'int16' [nx ny 1 n_frames nColors] 'image_total'},...
                        'writable',true);
                    % Initialize progress bar
                    spm_progress_bar('Init', length(fileNo), sprintf('Multispectral reading. Subject %s, Session %d\n',IOI.subj_name,s1), 'Files');
                    % Loop over files
                    for f1 = 1:length(fileNo)
                        frames = frameReadOut(frameReadOut(:,3)==fileNo(f1),:);
                        % extract certain frame
                        images=LectureImageMultiFast(dirs{s1},'image',frames);
                        if dupOn
                            for g1=1:length(images) %remove the duplicated pixels
                                images{g1} = images{g1}(:,1:2:end);
                            end
                        end
                        nImages = ceil(frames(end,1)/nColors)-iC;
                        %image_total{c1} = zeros(nx,ny,nz,nImages);
                        si = iC+1; %start index
                        ei = iC+nImages; %end index
                        %create an array for laser
                        if shrinkage_choice
                            laser_array = zeros(NX,NY,nImages);
                        end
                        sess.si = [sess.si si];
                        sess.ei = [sess.ei ei];
                        image_str_start = gen_num_str(si,nzero_padding);
                        image_str_last = gen_num_str(ei,nzero_padding);
                        image_str = [image_str_start 'to' image_str_last];
                        %separate images of different colors
                        for c1=1:nColors
                            str1 = str_color(c1);
                            fname{c1} = fullfile(dir_subj_res,sess_str, ...
                                [subj_name '_' OD_label '_' str1 '_' sess_str '_' image_str '.nii']);
                            sess.fname{c1} = [sess.fname{c1} fname{c1}];
                        end
                        
                        for fr1=1:(nColors*nImages)
                            tcol0 = mod(fr1,nColors); %which color for this image
                            if tcol0 == 0
                                tcol0 = nColors;
                            end
                            %which color tcol0 corresponds to:
                            for i0=1:length(color_order)
                                if strcmp(color_order(tcol0),str_color(i0))
                                    tcol = i0;
                                end
                            end
                            tframe = ceil(fr1/nColors); %which image for this frame
                            if ~(any(fr1+iC*nColors == missing_frames))
                                if ~(IOI.color.eng(tcol)==str_laser)
                                    tmp_image{tcol} = images{frames(frames(:,1)==(fr1+iC*nColors),2)};
                                    if shrinkage_choice
                                        % ioi_MYimresize replaces imresize
                                        tmp_image{tcol} = ioi_MYimresize(tmp_image{tcol}, [nx ny]);
                                        %                                         try
%                                             % Added by EGC: downsampling below
%                                             % does not always work, e.g. for
%                                             % NX=792, nx=396, but
%                                             % size(tmp_image{tcol},1) = 395 and
%                                             % there is a mismatch when assigning
%                                             % dimensions in the large memory
%                                             % mapped file im_obj
%                                             tmp_image{tcol} = imresize(tmp_image{tcol},[nx ny]);
%                                         catch exception
%                                             disp(exception.identifier)
%                                             disp(exception.stack(1))
%                                             % first spatially filter the images
%                                             tmp_image{tcol} = ioi_spatial_LPF('lpf', K, tmp_image{tcol});
%                                             % then downsample
%                                             tmp_image{tcol} = tmp_image{tcol}((SH.shrink_x/2+1):SH.shrink_x:(end-(SH.shrink_x/2)),(SH.shrink_y/2+1):SH.shrink_y:(end-(SH.shrink_y/2)));
%                                         end
                                    end
                                else
                                    %laser -- do not shrink now
                                    tmp_image{tcol} = images{frames(frames(:,1)==(fr1+iC*nColors),2)};
                                end
                            else
                                %missing frame: tmp_image{tcol} does not change, so the
                                %previous image will be used
                            end
                            if ~shrinkage_choice || (shrinkage_choice && ~(IOI.color.eng(tcol)==str_laser))
                                % im_obj.Data.image_total(:,:,:,tframe+iC,tcol) = tmp_image{tcol};
                                % Correct assignment for the 3rd dimension //EGC
                                im_obj.Data.image_total(:,:,1,tframe+iC,tcol) = tmp_image{tcol};
                            else
                                %store laser data
                                laser_array(:,:,tframe) = tmp_image{tcol}; %tframe???
                            end
                            %image_total{tcol}(:,:,:,tframe) = tmp_image{tcol};
                            if ~image_anat_first_pass && f1==1 && (fr1 <= nColors) && strcmp(str_color(tcol),str_anat)
                                image_anat = images{fr1}; %problem with fr1 sometimes...
                                image_anat_first_pass = 0;
                            end
                        end
                        iC = iC + nImages;
                        if shrinkage_choice
                            %save laser images
                            ioi_save_nifti(single(laser_array),sess.fname{IOI.color.eng == str_laser}{f1},vx);
                        end
                        % Update progress bar
                        spm_progress_bar('Set', f1);
                    end % Files loop
                    % Clear progress bar
                    spm_progress_bar('Clear');
                    
                    
                    
                    %Compute median
                    for c1=1:nColors
                        %skip laser
                        if ~(str_color(c1)==str_laser)
                            fprintf('Computing median. Subject %s, Session %d, Color %d...\n',IOI.subj_name,s1,c1);
                            median0{c1} = median(single(im_obj.Data.image_total(:,:,:,:,c1)),4);
                            %Save the median
                            str1 = str_color(c1);
                            sess.fname_median{c1} = fullfile(dir_subj_res,sess_str, ...
                                [subj_name '_' OD_label '_median_' str1 '_' sess_str '.nii']);
                            %ioi_save_nifti(single(median0{c1}),sess.fname_median{c1},vx);
                            tit0 = [subj_name ' ' OD_label ' median ' str1 ' ' sess_str];
                            ioi_save_images(single(median0{c1}),sess.fname_median{c1},vx,[],tit0);
                            try
                                %min and max and relative change
                                min_image = min(im_obj.Data.image_total(:,:,:,:,c1),[],4);
                                max_image = max(im_obj.Data.image_total(:,:,:,:,c1),[],4);
                                change = single(max_image) ./single(min_image);
                                %10th and 90th percentiles and relative change
                                tenthpctle_image = prctile(im_obj.Data.image_total(:,:,:,:,c1),10,4);
                                ninetiethpctle_image = prctile(im_obj.Data.image_total(:,:,:,:,c1),90,4);
                                change_90_10 = single(ninetiethpctle_image) ./ single(tenthpctle_image);
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
                        end
                    end % color loop
                    fprintf('Saving files. Subject %s, Session %d...\n',IOI.subj_name,s1);
                    for f1 = 1:length(fileNo)
                        ind0 = sess.si{f1}:sess.ei{f1};
                        %save images in nifti format
                        for c1=1:nColors
                            if ~(str_color(c1)==str_laser)
                                ioi_save_nifti(-log(single(im_obj.Data.image_total(:,:,:,ind0,c1))./repmat(single(median0{c1}),[1 1 1 length(ind0)])),sess.fname{c1}{f1},vx);
                            else
                                if ~shrinkage_choice
                                    ioi_save_nifti(single(im_obj.Data.image_total(:,:,:,ind0,c1)),sess.fname{c1}{f1},vx);
                                end
                            end
                        end
                    end
                    %disp(['Done processing session: ' int2str(s1) ', color: ' str1]);
                    IOI.sess_res{sC} = sess;
                    IOI.sess_res{sC}.hasRGY = hasRGY;
                    %sess_res =[sess_res sess];
                    disp(['Done processing session ' int2str(s1) ' images (' int2str(iC) ' images)']);
                    clear im_obj;
                    delete('all_images.dat');
                end
            end
        end % sessions for
        %IOI.sess_res = sess_res;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %4- Anatomical image
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %use first green colored image as the anatomy
        if ~expedite
            fname = fullfile(dir_subj_res, [subj_name '_' suffix_for_anat_file '.nii']);
            ioi_save_images(image_anat, fname, vx_anat,[],sprintf('%s Anatomical image',IOI.subj_name))
            %ioi_save_nifti(image_anat, fname, vx_anat);
            IOI.res.file_anat=fname;
        end
    end
catch exception
    disp(exception.identifier)
    disp(exception.stack(1))
end
