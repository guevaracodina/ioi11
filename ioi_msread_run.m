function out = ioi_msread_run(job)
%IOI structure: one per subject, good for multi-sessions
%IOI.dir: various directories
%IOI.sess_raw and IOI.sess_res: session specific raw and processed information
%IOI.res: results
%IOI.dev: device-level recording settings
out.IOImat = {};
%Color names - could be put in user interface?
str_anat = 'G';
str_laser = 'L'; %Optional
str_red = 'R'; %Varying number of colors allowed
str_green = 'G';
str_yellow = 'Y';
%If adding or removing colors, make sure the string of English and French
%colors have the same length and are a one-to-one mapping (code loops over colors
%to find images):
str_color = [str_red str_green str_yellow str_laser ]; %red, green, yellow and laser speckle
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
TR = 0.1999;
if length(job.force_redo) == 1
    force_redo = ones(1,length(job.top_bin_dir))*job.force_redo;
else
    force_redo = job.force_redo;
end

%               Directory Structure
%                   dir_group_all
%             /                     \
%     dir_group_raw              dir_group_res
%           |                       |
%     dir_subj_raw               dir_subj_res
%
dir_group_res_default = 'Res';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Big loop over subjects
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for SubjIdx=1:length(job.top_bin_dir)
    try
        tic
        clear IOI sess_raw sess_res tmp_image
        subj_OK = 1; %Boolean to exit loop if subject data is too corrupted
        %Directory structure
        dir_subj_raw = job.top_bin_dir{SubjIdx}; %location of raw data for this subject
        tmpsep = strfind(dir_subj_raw,filesep);
        dir_group_raw = dir_subj_raw(1:tmpsep(end-1)); %group level dir of raw data
        dir_group_all = dir_subj_raw(1:tmpsep(end-2)); %group level dir for all data
        subj_name = dir_subj_raw((tmpsep(end-1)+1):(tmpsep(end)-1)); %name of subject
        %path configuration
        if isfield(job.output_path_choice, 'output_path_select')
            dir_group_res = job.output_path_choice.output_path_select.output_path; %understood as group output path
        else
            dir_group_res = fullfile(dir_group_all,dir_group_res_default); %group level dir of processed data
        end
        dir_subj_res = fullfile(dir_group_res,subj_name);
        if (~exist(fullfile(dir_subj_res,'IOI.mat'),'file') || force_redo(SubjIdx))
            IOI.subj_name = subj_name;
            %Store directory information
            IOI.dir.dir_group_all = dir_group_all;
            IOI.dir.dir_group_raw = dir_group_raw;
            IOI.dir.dir_group_res = dir_group_res;
            IOI.dir.dir_subj_raw = dir_subj_raw;
            IOI.dir.dir_subj_res = dir_subj_res;
            IOI.color.eng = str_color;
            IOI.color.red = str_red;   
            IOI.color.green = str_green;
            IOI.color.yellow = str_yellow;
            IOI.color.laser = str_laser;
            IOI.res.shrinkageOn = 0;
            IOI.dev.TR = TR;
            %Create required directories
            if ~exist(dir_group_res,'dir'),mkdir(dir_group_res); end
            if ~exist(dir_subj_res,'dir'),mkdir(dir_subj_res); end
            [dummy,dirs] = cfg_getfile('FPListRec',dir_subj_raw,'.bin');
            sC = 0; %good session counter
            %loop over sessions
            for s1=1:length(dirs)
                sC = sC + 1;
                %session directory -- relabel
                sess_str = [sess_label gen_num_str(s1,2)];
                dir_sess = fullfile(dir_subj_res,sess_str);
                if ~exist(dir_sess,'dir'), mkdir(dir_sess); end
                %Read electro
                [files_el,dummy] = cfg_getfile('FPList',dirs{s1},[el_label '.' el_suffix]);
                [ConvertedData,ConvertVer,ChanNames,GroupNames,ci]=convertTDMS(1,files_el{1});
                fil_loc = fullfile(ConvertedData.FileFolder,[el_label '.mat']);
                %copy and remove
                IOI.sess_res{sC} = {};
                elfile_info = fullfile(dir_subj_res,[short_el_label el_info '_' sess_label gen_num_str(sC,2) '.mat']);
                elfile = fullfile(dir_subj_res,[short_el_label '_' sess_label gen_num_str(sC,2) '.mat']);
                IOI.res.el{sC} = elfile;
                %IOI.res.elinfo{sC} = elfile_info; %not used in later modules as did not exist in earlier reading module               
                copyfile(fil_loc, elfile_info);
                delete(fil_loc);
                IOI.res.electroOK = 1; %not by session...
                %extract ttl
                ttl=TDMS2ttl(ConvertedData);
                %save(elfile,el(1,:)); %incorrect. 
                %loop over image files for this session
                [imfiles dummy] = cfg_getfile('FPList',dirs{s1},'image');                
                %read frame present
                [Imout frameout frameReadOut fileNo] = LectureImageMultiFast(dirs{s1},'image',-1);
                missing_frames = [];
                for i0=1:frameReadOut(end,1)
                    tmp = find(i0==frameReadOut(:,1));
                    if isempty(tmp)
                        missing_frames = [missing_frames i0];
                    end
                end
                iC = 0;
                sess = [];
                for c1=1:nColors
                    sess.fname{c1} = {};
                end
                for f1 = 1:length(fileNo)
                    frames = frameReadOut(frameReadOut(:,3)==fileNo(f1),:);
                    % extract certain frame
                    images=LectureImageMultiFast(dirs{s1},'image',frames);
                    
                    [nx ny nz] = size(images{1});
                    %separate images of different colors
                    for c1=1:nColors
                        nImages = ceil(frames(end,1)/nColors)-iC;
                        image_total{c1} = zeros(nx,ny,nz,nImages);
                        str1 = str_color(c1);
                        
                        image_str_start = gen_num_str(iC+1,nzero_padding);
                        image_str_last = gen_num_str(iC+nImages,nzero_padding);
                        image_str = [image_str_start 'to' image_str_last];
                        fname{c1} = fullfile(dir_subj_res,sess_str, ...
                            [subj_name '_' OD_label '_' str1 '_' sess_str '_' image_str '.nii']);
                        sess.fname{c1} = [sess.fname{c1} fname{c1}]; 
                    end
      
                    for fr1=((nColors*iC)+1):(nColors*nImages)
                        tcol = mod(fr1,nColors)+1;
                        tframe = ceil(fr1/nColors)-iC;
                        if ~(any(fr1 == missing_frames))                            
                            tmp_image{tcol} = images{frames(frames(:,1)==fr1,2)};
                        else
                            %missing frame: tmp_image{tcol} does not change, so the
                            %previous image will be used 
                        end
                        image_total{tcol}(:,:,:,tframe) = tmp_image{tcol};
                        if s1==1 && f1==1 && (fr1 <= nColors) && strcmp(str_color(fr1),str_anat)
                            image_anat = images{fr1};
                        end
                    end   
                    iC = iC + nImages;
                    %save images in nifti format
                    for c1=1:nColors
                        ioi_save_nifti(image_total{c1},fname{c1},vx);
                    end
                end
                
                %disp(['Done processing session: ' int2str(s1) ', color: ' str1]); 
                IOI.sess_res{sC} = sess;
                IOI.sess_res{sC}.hasRGY = hasRGY;            
                %sess_res =[sess_res sess];
                disp(['Done processing session ' int2str(s1) ' images (' int2str(iC) 'images)']);
                
            end %for s1
        
            %IOI.sess_res = sess_res;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %4- Anatomical image
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %use first green colored image as the anatomy
            fname = fullfile(dir_subj_res, [subj_name '_' suffix_for_anat_file '.nii']);
            ioi_save_nifti(image_anat, fname, vx_anat);
            IOI.res.file_anat=fname;
            
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
        %toc
        disp(['Subject ' int2str(SubjIdx) ' complete']);
    catch exception
        disp(exception.identifier)
        disp(exception.stack(1))
    end
end


function IOI = disp_msg(IOI,msg)
disp(msg);
IOI.warning = [IOI.warning msg];