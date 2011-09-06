function out = ioi_msread_run(job)
%IOI structure: one per subject, good for multi-sessions
%IOI.dir: various directories
%IOI.sess_raw and IOI.sess_res: session specific raw and processed information
%IOI.res: results
%IOI.dev: device-level recording settings
block_size = 500; %number of images per block, if save_choice = 2 selected
%Hard coded the 6 frame steps for different colors, not sure it is always true
frames_per_cycle = 6;
out.IOImat = {};
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
for SubjIdx=1:length(job.subj.top_bin_dir)
    try
        tic
        clear IOI sess_raw sess_res
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
            %Create required directories
            if ~exist(dir_group_res,'dir'),mkdir(dir_group_res); end
            if ~exist(dir_subj_res,'dir'),mkdir(dir_subj_res); end
            
            %loop over sessions
            for s1=1:length()
                %Read electro
                [files_el,dummy] = cfg_getfile('FPList',dir_subj_raw,'electro.tdms');
                [ConvertedData,ConvertVer,ChanNames,GroupNames,ci]=convertTDMS(1,'electro.tdms');
                %extract ttl
                out=TDMS2ttl(ConvertedData);
                %read frame rpesent
                [Imout frameout frameReadOut fileNo] = LectureImageMultiFast(dir_subj_raw,'image',-1);
                % extract certain frame
                images=LectureImageMultiFast(dir_subj_raw,'image',frameReadOut(30:33,:));
                %Add found color - used later in concentration calculation module
                if ~(str1 == str_laser || str1 == str_contrast)
                    hasRGY = [hasRGY str1];
                end
                disp(['Done processing session: ' int2str(s1) ', color: ' str1]);
                
                
                IOI.sess_res{s1}.hasRGY = hasRGY;
                %electrophysiology
                if isfield(IOI.sess_raw{s1},'electroOK')
                    IOI = extract_electro(IOI,s1,IOI.sess_raw{s1}.electro{1});
                    IOI.res.electroOK = 1; %a bit early, but should be OK
                    disp(['Done processing session: ' int2str(s1) ', electrophysiology']);
                end
                %toc
                disp(['Done processing session ' int2str(s1) ' images (' int2str(n_frames) 'images)']);
                
            end %for s1
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


function IOI = disp_msg(IOI,msg)
disp(msg);
IOI.warning = [IOI.warning msg];