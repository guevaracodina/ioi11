function out = ioi_msioi_run(job)
%IOI structure: one per subject, good for multi-sessions
%IOI.dir: various directories
%IOI.sess_raw and IOI.sess_res: session specific raw and processed
%information; sess_raw is not used for new data format
%IOI.res: results
%IOI.dev: device-level recording settings
out.IOImat = {};
%               Directory Structure
%                   dir_group_all
%             /                     \
%     dir_group_raw              dir_group_res
%           |                       |
%     dir_subj_raw               dir_subj_res
%
dir_group_res_default = 'Res';
if length(job.force_redo) == 1
    force_redo = ones(1,length(job.top_bin_dir))*job.force_redo;
else
    force_redo = job.force_redo;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Big loop over subjects
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for SubjIdx=1:length(job.top_bin_dir)
    try
        %check which data format -- old or new
        tic
        clear IOI 
        subj_OK = 1; %Boolean to exit loop if subject data is too corrupted
        %Directory structure
        dir_subj_raw = job.top_bin_dir{SubjIdx}; %location of raw data for this subject
        tmpsep = strfind(dir_subj_raw,filesep);
        dir_group_raw = dir_subj_raw(1:tmpsep(end-1)); %group level dir of raw data
        dir_group_all = dir_subj_raw(1:tmpsep(end-2)); %group level dir for all data
        subj_name = dir_subj_raw((tmpsep(end-1)+1):(tmpsep(end)-1)); %name of subject
        job.SubjIdx = SubjIdx;
        % path configuration
        if isfield(job.output_path_choice, 'output_path_select')
            dir_group_res = job.output_path_choice.output_path_select.output_path{:}; % understood as group output path
        else
            dir_group_res = fullfile(dir_group_all,dir_group_res_default); %group level dir of processed data
        end
        dir_subj_res = fullfile(dir_group_res,subj_name);
        if (~exist(fullfile(dir_subj_res,'IOI.mat'),'file') || force_redo(SubjIdx))
            %Use existing IOI.mat in case force_redo==2 (Partial) was chosen
            if (exist(fullfile(dir_subj_res,'IOI.mat'),'file') && force_redo(SubjIdx)==2)
                load(fullfile(dir_subj_res,'IOI.mat'));
                %Boolean for partial reprocessing of data - useful
                %if some raw data was found to be corrupted
                %but was somehow corrected -- cannot add sessions this way
                job.PartialRedo2 = 1;
            else
                job.PartialRedo2 = 0;
            end
            if ~job.PartialRedo2
                %fill IOI structure
                IOI.warning = {};
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
            end
            %resume files - subject level
            [files_txt,dummy] = cfg_getfile('FPList',dir_subj_raw,'.txt');
            if ~isempty(files_txt) && ~isempty(strfind(files_txt{1},'resume_manip.txt'))
                data_format = 0; %old format
            else
                if ~isempty(files_txt) && (~isempty(strfind(files_txt{1},'resume.txt')) || ...
                        ~isempty(strfind(files_txt{1},'resumeExp.txt')))
                    data_format = 1; %new format
                else
                    currDir = pwd; cd(dir_subj_raw);
                    isBLK = ~system('dir *.blk /s');
                    isSam = ~system('dir *.seq /s');
                    isMHI = ~system('dir *.tif /s');
                    
                    if isBLK
                        data_format = 2; %blk format (Casanova's lab)
                    else
                        if isSam
                            data_format = 3; %Sam format (Sam's setup)
                        else
                            if isMHI
                                data_format = 4; % MHI format (MHI setup single wavelnegth)
                            else
                                subj_OK = 0;
                                disp(['Data format unrecognized for Subject ' int2str(SubjIdx)]);
                                disp('Start by checking if resume.txt (new) or resume_manip.txt (old) file is missing');
                            end
                        end
                    end
                    cd(currDir)
                end
            end
            if subj_OK
                IOI.subj_OK = subj_OK;
                if data_format == 1
                    fprintf('Reading IOI new format (Montreal Heart Institute)...\n');
                    IOI = ioi_msread_new_format(IOI,job);
                elseif data_format == 2
                    % blk format (Casanova's lab)
                    fprintf('Reading IOI .BLK format (Optometry school)...\n');
                    IOI = ioi_msread_blk_format(IOI,job);
                elseif data_format == 3
                    % Sam's IOI setup
                    fprintf('Reading IOI .SEQ format (Sam''s setup)...\n');
                    IOI = ioi_msread_sam_format(IOI,job);
                elseif data_format == 4
                    % Sam's IOI setup
                    fprintf('Reading IOI .TIF format (MHI setup, single wavelength)...\n');
                    IOI = ioi_msread_mhi_format(IOI,job);
                else
                    fprintf('Reading IOI old format (Sacre-Coeur)...\n');
                    IOI = ioi_msread_old_format(IOI,job);
                end
                % Write out IOI in .mat file format at the subject results path location
                % All (links to) data are preserved in this structure
                IOImat = fullfile(dir_subj_res,'IOI.mat');
                save(IOImat,'IOI');
            end                 
        else
            %IOI.mat already exists, and rerun is not enforced - therefore just
            %skip this module and pass IOImat structure to the next module
            IOImat = fullfile(dir_subj_res,'IOI.mat');
        end
        if subj_OK
            out.IOImat = [out.IOImat IOImat];
        end
        disp(['Elapsed time: ' datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS')]);
        % disp(['Subject ' int2str(SubjIdx) ' complete']);
        disp(['Subject ' int2str(SubjIdx) ' (' subj_name ')' 'complete']);
    catch exception
        disp(exception.identifier)
        disp(exception.stack(1))
        out.IOImat = [out.IOImat dir_subj_raw];
    end
end

