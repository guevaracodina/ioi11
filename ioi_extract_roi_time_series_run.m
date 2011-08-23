function out = ioi_extract_roi_time_series_run(job)
%select a subset of sessions
if isfield(job.session_choice,'select_sessions')
    all_sessions = 0;
    selected_sessions = job.session_choice.select_sessions.selected_sessions;
else
    all_sessions = 1;
end
%select a subset of ROIs
if isfield(job.ROI_choice,'select_ROIs')
    all_ROIs = 0;
    selected_ROIs = job.ROI_choice.select_ROIs.selected_ROIs;
else
    all_ROIs = 1;
end
for SubjIdx=1:length(job.IOImat)
    try
        tic
        clear IOI ROI
        %Load IOI.mat information
        IOImat = job.IOImat{SubjIdx};
        load(IOImat);
        
        if ~isfield(IOI.res,'ROIOK')
            disp(['No ROI available for subject ' int2str(SubjIdx) ' ... skipping series extraction']);
        else
            if ~isfield(IOI.res,'seriesOK') || job.force_redo
                [dir_ioimat dummy] = fileparts(job.IOImat{SubjIdx});
                %loop over ROIs
                for r1=1:length(IOI.res.ROI)
                    if all_ROIs || sum(r1==selected_ROIs)
                        vol = spm_vol(IOI.res.ROI{r1}.fname);
                        tmp_mask = logical(spm_read_vols(vol));
                        %shrink mask to voxel size
                        if IOI.res.shrinkageOn
                            sz = size(tmp_mask);
                            %careful, this floor might not lead to the
                            %correct size - better to do a check on image
                            %size
                            mask{r1} = imresize(tmp_mask,[floor(sz(1)/IOI.res.shrink_x) floor(sz(2)/IOI.res.shrink_y)],'bicubic');
                        else
                            mask{r1} = tmp_mask;
                        end
                    end
                end
                
                %loop over sessions
                for s1=1:length(IOI.sess_res)
                    if all_sessions || sum(s1==selected_sessions)
                        %loop over available colors
                        for c1=1:length(IOI.sess_res{s1}.fname)
                            if ~(IOI.color.eng(c1)==IOI.color.laser)
                                %skip laser - only extract for flow
                                fname_list = IOI.sess_res{s1}.fname{c1};
                                
                                %initialize
                                for r1=1:length(IOI.res.ROI)
                                    if all_ROIs || sum(r1==selected_ROIs)
                                        ROI{r1}{s1,c1} = [];
                                    end
                                end
                                %loop over files
                                for f1=1:length(fname_list)
                                    fname = fname_list{f1};
                                    vols = spm_vol(fname);
                                    d = spm_read_vols(vols);
                                    [d1 d2 d3 d4] = size(d);
                                    %time dimension in 3rd dimension for colors
                                    %R, G, Y, but in 4th dimension for O, D, F
                                    %Loop over ROIs
                                    for r1=1:length(IOI.res.ROI)
                                        if all_ROIs || sum(r1==selected_ROIs)
                                            for i3=1:d3
                                                for i4=1:d4
                                                    %extracted data
                                                    %tmp_d = squeeze(d(:,:,i3,i4));
                                                    tmp_d = d(:,:,i3,i4);
                                                    %just take mean over mask for now
                                                    if ~(IOI.color.eng(c1)==IOI.color.contrasts)
                                                        e = mean(tmp_d(mask{r1}));
                                                    else
                                                        %contrast images will be smaller and need to be resized
                                                        tmask = imresize(mask{r1},[d1 d2]);
                                                        e = mean(tmp_d(tmask));
                                                    end
                                                    ROI{r1}{s1,c1} = [ROI{r1}{s1,c1} e];
                                                end
                                            end
                                        end
                                    end
                                end
                                disp(['ROIs for session ' int2str(s1) ' and color ' IOI.color.eng(c1) ' completed']);
                            end
                        end
                    end
                end
                
                IOI.res.seriesOK = 1;
                if isfield(job.IOImatCopyChoice,'IOImatCopy')
                    newDir = job.IOImatCopyChoice.IOImatCopy.NewIOIdir;
                    newDir = fullfile(dir_ioimat,newDir);
                    if ~exist(newDir,'dir'),mkdir(newDir); end
                    IOImat = fullfile(newDir,'IOI.mat');
                end
                [dir1 dummy] = fileparts(IOImat);
                ROIfname = fullfile(dir1,'ROI.mat');
                save(ROIfname,'ROI');
                IOI.ROI.ROIfname = ROIfname;
                save(IOImat,'IOI');
            end
        end
        toc
        disp(['Subject ' int2str(SubjIdx) ' complete']);
        out.IOImat{SubjIdx} = IOImat;
    catch exception
        disp(exception.identifier)
        disp(exception.stack(1))
    end
end