function out = ioi_create_roi_run(job)
if isfield(job.AutoROIchoice,'AutoROI')
    ROIsize = job.AutoROIchoice.AutoROI.ArrayROI;
    autoROI = 1;
else
    autoROI = 0;
end
for SubjIdx=1:length(job.IOImat)
    try
        clear IOI
        %Load IOI.mat information
        IOImat = job.IOImat{SubjIdx};
        load(IOImat);
        
        if ~isfield(IOI.res,'ROIOK') || job.force_redo
            [dir_ioimat dummy] = fileparts(job.IOImat{SubjIdx});
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
                        title('Make ROI polygon, then double click in it to create ROI.');
                        mask = roipoly;
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
            IOI.res.ROIOK = 1;
            if isfield(job.IOImatCopyChoice,'IOImatCopy')
                newDir = job.IOImatCopyChoice.IOImatCopy.NewIOIdir;
                newDir = fullfile(dir_ioimat,newDir);
                if ~exist(newDir,'dir'),mkdir(newDir); end
                IOImat = fullfile(newDir,'IOI.mat');
            end
            save(IOImat,'IOI');
        end
        out.IOImat{SubjIdx} = IOImat;
    catch exception
        disp(exception.identifier)
        disp(exception.stack(1))
    end
end