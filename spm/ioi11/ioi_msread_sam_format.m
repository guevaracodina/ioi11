function IOI = ioi_msread_sam_format(IOI,job)
% Reads and converts sam data into NIfTI format

try
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
    wsize=job.window_size;
    if (mod(wsize,2)==0)
        disp(['Contrast Computation: Need to specify an odd window size, changing to: ' num2str(wsize+1)]);
        wsize=wsize+1;
        job.window_size=wsize;  % For reference
    end
    % XXX: This is potentially an error: No assignment to IOI should occur if we
    % are not redoing
    IOI.res.flow.window_size=wsize;
    
    if isfield(job,'stim_cutoff')
        stim_cutoff = job.stim_cutoff;
    else
        stim_cutoff = 5;
    end
    if isfield(job,'minTimeBetweenStim')
        minTimeBetweenStim = job.minTimeBetweenStim;
    else
        minTimeBetweenStim = 1.001;
    end
    
    %subj_OK = 1; %Boolean to exit loop if subject data is too corrupted
    if ~job.PartialRedo2
        IOI.color.eng = str_color;
        IOI.color.red = str_red;
        IOI.color.green = str_green;
        IOI.color.yellow = str_yellow;
        IOI.color.laser = str_laser;
        IOI.res.shrinkageOn = 0; %no shrinkage of images
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
                % Read Sam images
                [frames, data] = ioi_read_sam_data(dirs{s1});
                %read info file
                try
                    [file_info,dummy] = cfg_getfile('FPList',dirs{s1},'info.txt');
                    fid = fopen(file_info{1},'r');
                    t0 = fread(fid,'uint8');
                    fclose(fid);
                    t0 = char(t0');
                    t0 = deblank(t0);
                    % Get acquisition frequency directly from the file //EGC
                    [startIndex endIndex]= regexp(t0, 'Frequency :', 'once');
                    tmpString = t0(endIndex:end);
                    % Get acquisition frequency as a number //EGC
                    [startIndex endIndex] = regexp(tmpString,'\d+','once');
                    acq_freq=str2double(tmpString(startIndex:endIndex));
                    % Find the phrase color order //EGC
                    [startIndex endIndex] = regexp(t0, 'Color order : ', 'once');
                    % Look for the remaining characters (more robust) //EGC
                    color_order_french = t0(endIndex+1:end);
                    % color_order_french = t0(end-3:end);
                    % if strcmp(color_order_french,'flat')
                    switch color_order_french
                        case 'flat'
                            nColors=4;
                            try
                                sp0 = strfind(t0,'Light : ');
                                sr0 = t0((sp0+8):(sp0+19));
                                color_order_french = [sr0(1) sr0(4) sr0(7) sr0(10)];
                                %These are the number of frames for each color
                                %not used currently, but will be for calcium
                                %imaging
                                color_numbers = [sr0(2) sr0(5) sr0(8) sr0(11)];
                            catch
                                color_order_french = 'RVJL';
                            end
                        case 'version1'
                            % This means that the following lines have been
                            % manually replaced in the info.txt file
                            % Light : R1,V1,J1,L1
                            % Color order : flat
                            % By
                            % Light : V1,J1,L1,R1
                            % Color order : version1
                        case 'Version2, Txt flash'
                            % It is supossed to be the correct LabView code
                            try
                                % Find the phrase Light //EGC
                                [startIndex endIndex] = regexp(t0, 'Light : ', 'once');
                                tmpString = t0(endIndex:end);
                                % Find the next end of line
                                eolIndex = regexp(tmpString, '\n', 'once');
                                % Get the colors string if \r eolIndex-1
                                tmpString = strtrim(tmpString(1:eolIndex));
                                % Find the letters
                                letterIndex = regexp(tmpString, '[A-Z]');
                                % Find the numbers
                                numberIndex = regexp(tmpString, '[0-9]');
                                color_order_french = tmpString(letterIndex);
                                %These are the number of frames for each color
                                %not used currently, but will be for calcium
                                %imaging
                                color_numbers = tmpString(numberIndex);
                                tmpString = [];
                                for iString = 1:numel(color_order_french)
                                    tmpString = [ tmpString, repmat(color_order_french(iString), [1 str2double(color_numbers(iString))]) ];
                                end
                                % Now each character represents a color
                                color_order_french = tmpString;
                                % Get the real number of colors
                                nColors = numel(color_order_french);
                            catch
                                color_order_french = 'RVJL';
                                nColors=4;
                            end
                    end
                    % Save repetition time in IOI //EGC
                    IOI.dev.TR = nColors/acq_freq;

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
                    nColors = numel(unique(frames.Color));
                    color_order = frames.Color(1:nColors);
%                     color_order = str_color_default;
                    IOI.dev.TR = 1/frames.FreqR;
                    IOI = disp_msg(IOI,['No info.txt for session ' int2str(s1) ';TR set to ' num2str(IOI.dev.TR) ' and  color order: ' color_order]);    
                end
                %Read electro
                [files_el,dummy] = cfg_getfile('FPList',dirs{s1},[el_label '.' el_suffix]);
                if ~isempty(files_el)
                    [ConvertedData,ConvertVer,ChanNames,GroupNames,ci]=convertTDMS(1,files_el{1});
                    fil_loc = fullfile(ConvertedData.FileFolder,[el_label '.mat']);
                    %Get default stimulations -- what if there are several
                    %types of stimulations???
                    stims = ConvertedData.Data.MeasuredData(6).Data;
                    stims_dt= ConvertedData.Data.MeasuredData(6).Property(3).Value;
                    
                    ons0 = stims_dt*find(stims>stim_cutoff);
                    par0 = stims(stims>stim_cutoff);
                    ons1 = find(diff(ons0)>=minTimeBetweenStim);
                    %par1 = par0(ons1);
                    if ~isempty(ons1)
                        ons = ons0([1 1+ons1']);
                        par = par0([1 1+ons1']);
                    else
                        ons = []; %in case we have a rest session
                        par = [];
                    end
                    clear names onsets durations parameters
                    %Converted to seconds, rather than frame number
                    for stim_index=1:1
                        names{stim_index}=['Stim_',num2str(stim_index)];
                        onsets{stim_index}=ons;
                        durations{stim_index}=ones(1,length(ons));
                        try
                            parameters{stim_index}=par;
                        catch
                            parameters{stim_index}=[];
                        end
                    end
                    %Store onset information
                    %sess.list_stim = list_stim;
                    sess.names = names; %not quite the SPM format, need to put in Sess.U format
                    sess.onsets = onsets; %in seconds
                    sess.durations = durations; %in seconds
                    sess.parameters = parameters;
                end
                % Dummy electrophysiology files
                elfile_info = fullfile(dir_subj_res,[short_el_label el_info '_' sess_label gen_num_str(sC,2) '.mat']);
                elfile = fullfile(dir_subj_res,[short_el_label '_' sess_label gen_num_str(sC,2) '.mat']);
                IOI.res.el{sC} = elfile;
                IOI.res.elinfo{sC} = elfile_info; %not used in later modules as did not exist in earlier reading module
                elfile_info2 = fullfile(dir_subj_res,[short_el_label '2' el_info '_' sess_label gen_num_str(sC,2) '.mat']);
                elfile2 = fullfile(dir_subj_res,[short_el_label '2_' sess_label gen_num_str(sC,2) '.mat']);
                %save 2nd channel of electrophysiology
                IOI.res.el2{sC} = elfile2;
                
                if ~isempty(files_el)
                    copyfile(fil_loc, elfile_info);
                    delete(fil_loc);
                    %save electrophysiological data
                    el = ConvertedData.Data.MeasuredData(5).Data;
                    save(elfile,'el');
                    try
                        elfile_info2 = fullfile(dir_subj_res,[short_el_label '2' el_info '_' sess_label gen_num_str(sC,2) '.mat']);
                        elfile2 = fullfile(dir_subj_res,[short_el_label '2_' sess_label gen_num_str(sC,2) '.mat']);
                        %save 2nd channel of electrophysiology
                        %*****************************************************
                        IOI.res.el2{sC} = elfile2;
                        %*****************************************************
                        el2 = ConvertedData.Data.MeasuredData(12).Data;
                        save(elfile2,'el2');
                    end
                    IOI.res.sfel = 1/stims_dt;
                end
               
                %if expedite
                dir_elfig = fullfile(dir_subj_res,'fig_el');
                if ~exist(dir_elfig,'dir'), mkdir(dir_elfig); end
                fname_el = fullfile(dir_elfig,[short_el_label '1_' sess_label gen_num_str(sC,2)]);
                fname_el2 = fullfile(dir_elfig,[short_el_label '1&2_' sess_label gen_num_str(sC,2)]);
                fname_el3 = fullfile(dir_elfig,[short_el_label '2_' sess_label gen_num_str(sC,2)]);
                % This saves figures of electrophysiology to inspect data
                % later (try is for safery in case they don't exist)
                try ioi_plot_LFP(IOI,el,s1,1,1,fname_el); end;
                try ioi_plot_LFP(IOI,el2,s1,1,1,fname_el3); end;
                try ioi_plot_LFP2(IOI,el,s1,el2,1,1,fname_el2); end;
                %end
                if ~expedite
                    if ~isempty(files_el)
                        %extract ttl
                        ttl=TDMS2ttl(ConvertedData);
                    end
                    %loop over image files for this session
                    
                    fileNo = 1;
                    missing_frames = [];
                    IOI.res.sfel = 10000;
                    
                    frames.LUT_Red      = frames.LUT_Red(frames.FrameToSkip_Start+1:end-frames.FrameToSkip_End);
                    frames.LUT_Green    = frames.LUT_Green(frames.FrameToSkip_Start+1:end-frames.FrameToSkip_End);
                    frames.LUT_Yellow   = frames.LUT_Yellow(frames.FrameToSkip_Start+1:end-frames.FrameToSkip_End);
                    
                    n_frames = min([numel(frames.LUT_Red) numel(frames.LUT_Green) numel(frames.LUT_Yellow)]);
                    
                    frames.LUT_Red      = frames.LUT_Red(1:n_frames);
                    frames.LUT_Green    = frames.LUT_Green(1:n_frames);
                    frames.LUT_Yellow   = frames.LUT_Yellow(1:n_frames);
                    
                    sess.n_frames = n_frames;

                    for c1=1:nColors
                        sess.fname{c1} = {};
                    end
                    sess.si = {};
                    sess.ei = {};
                    %Get image size
                    nx = data.format{1,2}(1);
                    ny = data.format{1,2}(2);
                    
                    %shrink if required
                    [shrinkage_choice SH] = ioi_get_shrinkage_choice(job); 
                    IOI.res.shrinkageOn = shrinkage_choice; 
                    if shrinkage_choice
                        % Only if shrinkage is chosen //EGC
                        IOI.res.shrink_x = SH.shrink_x;
                        IOI.res.shrink_y = SH.shrink_y;    
                        vx =[IOI.res.shrink_x IOI.res.shrink_y 1];
                    end
                    
                    NX = nx;
                    NY = ny;
                    if shrinkage_choice
                        nx = floor(NX/SH.shrink_x);
                        ny = floor(NY/SH.shrink_y);
                        K.radius = SH.shrink_x;
                        K.k1 = NX;
                        K.k2 = NY;
                        K = ioi_spatial_LPF('set', K);
                    end
                    
                    %create large memmapfile for Y, R, G L images (several GB)
                    fmem_name = 'all_images.dat';
                    fid = fopen(fmem_name,'w');
                    fwrite(fid, zeros([nx ny 1 n_frames nColors]), 'int16');
                    fclose(fid);
                    im_obj = memmapfile(fmem_name,'format',...
                        {'int16' [nx ny 1 n_frames nColors] 'image_total'},...
                        'writable',true);
                    % Initialize progress bar
%                     spm_progress_bar('Init', length(fileNo), sprintf('Multispectral reading. Subject %s, Session %d\n',IOI.subj_name,s1), 'Files');
                    % Loop over files
                    iC = 0; %counter for number of images so far
                    
                    
                    for f1 = 1:length(fileNo)
                        
                        image_str = [gen_num_str(1,nzero_padding) 'to' gen_num_str(sess.n_frames,nzero_padding)];
                        
                        for c1=1:nColors
                            str1 = str_color(c1);
                            if ~(str1==str_laser)
                                fname{c1} = fullfile(dir_subj_res,sess_str, ...
                                    [subj_name '_' OD_label '_' str1 '_' sess_str '_' image_str '.nii']);
                            else
                                fname{c1} = fullfile(dir_subj_res,sess_str, ...
                                    [subj_name '_C_' sess_str '_' image_str '.nii']);
                            end
                            sess.fname{c1} = [sess.fname{c1} fname{c1}];
                            
                        end
                        
                        % Initialize progress bar
                        spm_progress_bar('Init', sess.n_frames, ...
                            sprintf('Multispectral reading. Subject %s, Session %d\n',IOI.subj_name,s1), 'Frames');
                        
                        % color_order = 'RGY' //EGC
                        for iFrame = 1:sess.n_frames
                            if shrinkage_choice
                                im_obj.Data.image_total(:,:,1,iFrame,1) = ioi_MYimresize(squeeze(data.Data(frames.LUT_Red(iFrame)).framej), [nx ny]);
                                im_obj.Data.image_total(:,:,1,iFrame,2) = ioi_MYimresize(squeeze(data.Data(frames.LUT_Green(iFrame)).framej), [nx ny]);
                                im_obj.Data.image_total(:,:,1,iFrame,3) = ioi_MYimresize(squeeze(data.Data(frames.LUT_Yellow(iFrame)).framej), [nx ny]);
                            else
                                im_obj.Data.image_total(:,:,1,iFrame,1) = squeeze(data.Data(frames.LUT_Red(iFrame)).framej);
                                im_obj.Data.image_total(:,:,1,iFrame,2) = squeeze(data.Data(frames.LUT_Green(iFrame)).framej);
                                im_obj.Data.image_total(:,:,1,iFrame,3) = squeeze(data.Data(frames.LUT_Yellow(iFrame)).framej);
                            end
                            % Update progress bar
                            spm_progress_bar('Set', iFrame);
                        end
                        % Clear progress bar
                        spm_progress_bar('Clear');
                        iC = iC + sess.n_frames;
                    end % Files loop
              
                    %Compute median for colors other than laser
                    for c1=1:nColors
                        %skip laser
                        if ~(str_color(c1)==str_laser)
                            fprintf('Computing median. Subject %s, Session %d, Color %d...\n',IOI.subj_name,s1,c1);
                            median0{c1} = median(single(im_obj.Data.image_total(:,:,:,:,c1)),4);
                            %Save the median
                            str1 = str_color(c1);
                            sess.fname_median{c1} = fullfile(dir_subj_res,sess_str, ...
                                [subj_name '_' OD_label '_median_' str1 '_' sess_str]);
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
                                ioi_save_images(single(min_image),sess.fname_min{c1},vx,[],tit1);
                                ioi_save_images(single(max_image),sess.fname_max{c1},vx,[],tit2);
                                ioi_save_images(single(change),sess.fname_change{c1},vx,[],tit3);
                                ioi_save_images(single(tenthpctle_image),sess.fname_10pctle{c1},vx,[],tit4);
                                ioi_save_images(single(ninetiethpctle_image),sess.fname_90pctle{c1},vx,[],tit5);
                                ioi_save_images(single(change_90_10),sess.fname_change_90_10{c1},vx,[],tit6);
                            end
                        end
                    end % color loop
                    c1 = find(str_color == str_green);
                    image_anat = single(data.Data(frames.LUT_Green(1)).framej);
                    fprintf('Saving files. Subject %s, Session %d...\n',IOI.subj_name,s1);
                    for f1 = 1:length(fileNo)
                        %save images in nifti format
                        for c1=1:nColors
                            if ~(str_color(c1)==str_laser)
                                ioi_save_nifti(-log(single(im_obj.Data.image_total(:,:,:,:,c1))./repmat(single(median0{c1}),[1 1 1 sess.n_frames])),sess.fname{c1}{f1},vx);
                            else
                                % Here the 10000 is because the INT16
                                % memmapfile is an issue with contrast
                                % between 0-1
                                ioi_save_nifti(single(im_obj.Data.image_total(:,:,:,ind0,c1))/10000,sess.fname{c1}{f1},vx);
                            end
                        end
                    end
                    IOI.sess_res{sC} = sess;
                    IOI.sess_res{sC}.hasRGY = hasRGY;
                    % Save available colors
                    IOI.sess_res{sC}.availCol = color_order;
                    % Stats on missing frames 
                    IOI.sess_res{sC}.missingFrames = numel(missing_frames);
                    disp(['Done processing session ' int2str(s1) ' images (' int2str(iC) ' images)']);
                    clear im_obj;
                    delete('all_images.dat');
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %4- Anatomical image - Save for each session
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Use first green colored image as the anatomy
                    anat_fname = fullfile(dir_subj_res,sess_str, ...
                                    [subj_name '_' suffix_for_anat_file '_' sess_str]);
                    ioi_save_images(image_anat, anat_fname, vx_anat,[],sprintf('%s Anatomical image',IOI.subj_name))
                    % It will always point to the last anatomical image
                    IOI.res.file_anat = [anat_fname '.nii'];
                end
            end
        end % sessions for
    end
catch exception
    disp(exception.identifier)
    disp(exception.stack(1))
end

% EOF
