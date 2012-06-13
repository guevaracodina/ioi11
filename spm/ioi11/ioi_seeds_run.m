function out = ioi_seeds_run(job)
if isfield(job.AutoROIchoice,'AutoROI')
    ROIsize = job.AutoROIchoice.AutoROI.ArrayROI;
    autoROI = 1;
else
    autoROI = 0;
    if isfield(job.AutoROIchoice,'ManualROI')
        graphicalROI = 1;
    else
        graphicalROI = 0;
    end
end
for SubjIdx=1:length(job.IOImat)
    try
        clear IOI
        %Load IOI.mat information
        IOImat = job.IOImat{SubjIdx};
        [dir_ioimat dummy] = fileparts(job.IOImat{SubjIdx});
        if isfield(job.IOImatCopyChoice,'IOImatCopy')
            newDir = job.IOImatCopyChoice.IOImatCopy.NewIOIdir;
            newDir = fullfile(dir_ioimat,newDir);
            if ~exist(newDir,'dir'),mkdir(newDir); end
            IOImat = fullfile(newDir,'IOI.mat');
        else
            newDir = dir_ioimat;
        end
        try
            load(IOImat);
        catch
            load(job.IOImat{SubjIdx});
        end
        if ~isfield(IOI.res,'ROIOK') || job.force_redo
            if job.RemovePreviousROI
                try
                    IOI.res = rmfield(IOI.res,'ROIOK');
                end
                try
                    for i1=1:length(IOI.res.IOI)
                        %clean up: delete ROI mask files
                        delete(IOI.res.ROI{i1}.fname);
                    end
                end
                try
                    IOI.res = rmfield(IOI.res,'ROI');
                end
                index = 0;
            else
                if isfield(IOI.res,'ROI')
                    index = length(IOI.res.ROI);
                else
                    index = 0;
                end
            end
            %display anatomical image
            vol = spm_vol(IOI.res.file_anat);
            vx = [1 1 1];
            [dir1 fil1] = fileparts(vol.fname);
            im_anat = spm_read_vols(vol);
            if ~autoROI
                %Display images of changes from 10th to 90th percentile for all sessions
                try
                    for i0=1:length(IOI.sess_res)
                        hs{i0} = figure;
                        V = spm_vol(IOI.sess_res{i0}.fname_change_90_10{1}); %color green
                        tmp_image = spm_read_vols(V);
                        imagesc(tmp_image);
                        title(['Session ' int2str(i0) ': ratio of 90th to 10th percentile']);
                    end
                end
                p = 1;
                h2 = figure; spm_input(['Subject ' int2str(SubjIdx)],'-1','d');
                full_mask = ones(size(im_anat));
                linecount = 0;
                while p
                    figure(h2);
                    p = spm_input('Add an ROI?',2*linecount+2,'y/n');
                    if p == 'y', p = 1; else p = 0; end
                    if p
                        h1 = figure('Position',[20 50 3*size(im_anat,1) 3*size(im_anat,2)]);
                        colormap(gray); imagesc(im_anat.*full_mask);
                        if graphicalROI
                            title('Make ROI polygon, then double click in it to create ROI.');
                            mask = roipoly;
                        else
                            linecount = linecount + 1;
                            rc = spm_input('Enter [row,column] of center',2*linecount+1,'e',[],2);
                            linecount = linecount + 1;
                            radius = spm_input('Enter radius in pixels',2*linecount+1,'e',[],1);
                            radius = round(radius);
                            if radius < 0, radius = 0; end
                            mask = zeros(size(im_anat));
                            for x1=-radius:radius
                                for y1=-radius:radius
                                    if x1^2+y1^2 <= radius^2
                                        try %will skip pixels outside the image
                                            mask(rc(1)+y1,rc(2)+x1) = 1;
                                        end
                                    end
                                end
                            end
                        end
                        mask = single(mask);
                        full_mask = full_mask-mask;
                        index = index + 1; linecount = linecount + 1;
                        %imagesc(mask);
                        if job.select_names
                            figure(h2);
                            name = spm_input(['Enter name of ROI' int2str(index)],2*linecount+1,'s');
                        else
                            name = int2str(index);
                        end
                        IOI.res.ROI{index}.name = name;
                        if index < 10, str0 = '0'; else str0 = ''; end
                        str = [str0 int2str(index)];
                        fname_mask = fullfile(dir1,[fil1 '_ROI_' str '.nii']);
                        IOI.res.ROI{index}.fname = fname_mask;
                        ioi_save_nifti(mask, fname_mask, vx);
                    end
                end
                try close(h1); end
                try close(h2); end
                for i0=1:length(IOI.sess_res)
                    try close(hs{i0}); end
                end
            else
                %automatic ROIs
                %h1 = figure('Position',[20 50 3*size(im_anat,1) 3*size(im_anat,2)]);
                sz = size(im_anat);
                sz = sz(1:2); %remove 3rd component, which is one
                sz = floor(sz./ROIsize);
                for i1 = 1:ROIsize(1) %N
                    for i2 = 1:ROIsize(2) %M
                        %ROI created in following order
                        % 1      2       3 ...   ROIsize(2)=M
                        % M+1   M+2     M+3 ...  2*M
                        %  .    .        .       .
                        %  .    .        .       N*M
                        index = index+1;
                        %vertices in the order: UL, UR, LR, LL
                        mask = roipoly(im_anat,...
                            [1+sz(2)*(i2-1) 1+sz(2)*i2 1+sz(2)*i2 1+sz(2)*(i2-1)],...
                            [1+sz(1)*(i1-1) 1+sz(1)*(i1-1) 1+sz(1)*i1 1+sz(1)*i1]);
                        mask = single(mask);
                        if index < 10, str0 = '0'; else str0 = ''; end
                        str = [str0 int2str(index)];
                        fname_mask = fullfile(dir1,[fil1 '_ROI_' str '_' int2str(i1) 'x' int2str(i2) '.nii']);
                        IOI.res.ROI{index}.fname = fname_mask;
                        ioi_save_nifti(mask, fname_mask, vx);
                    end
                end
            end
            IOI.ROIname = {};
            for i0=1:length(IOI.res.ROI)
                if isfield(IOI.res.ROI{i0},'name')
                    IOI.ROIname = [IOI.ROIname; IOI.res.ROI{i0}.name];
                else
                    IOI.ROIname = [IOI.ROIname; ['ROI' gen_num_str(i0,3)]];
                end
            end
            IOI.res.ROIOK = 1;
            if isfield(job.IOImatCopyChoice,'IOImatCopy')
                newDir = job.IOImatCopyChoice.IOImatCopy.NewIOIdir;
                newDir = fullfile(dir_ioimat,newDir);
                if ~exist(newDir,'dir'),mkdir(newDir); end
                IOImat = fullfile(newDir,'IOI.mat');
            end
            save(IOImat,'IOI');
        end
        disp(['Subject ' int2str(SubjIdx) ' complete']);
        out.IOImat{SubjIdx} = IOImat;
    catch exception
        disp(exception.identifier)
        disp(exception.stack(1))
        out.IOImat{SubjIdx} = job.IOImat{SubjIdx};
    end
end
