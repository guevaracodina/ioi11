function out = ioi_create_roi_run(job)
try
    use_gray_contrast = job.use_gray_contrast;
catch
    use_gray_contrast = 0;
end
% Manual/automatic selection of ROIs/seeds (pointNclickROI ManualROIspline)

% if isfield(job.AutoROIchoice,'AutoROI')
%     ROIsize = job.AutoROIchoice.AutoROI.ArrayROI;
%     autoROI = 1;
% else
%     autoROI = 0;
%     if isfield(job.AutoROIchoice,'ManualROI')
%         graphicalROI = 1;
%     else
%         graphicalROI = 0;
%     end
% end
fNamesROIchoice = fieldnames(job.AutoROIchoice);
switch(fNamesROIchoice{1})
    case 'AutoROI'
        autoROI         = 1;
        graphicalROI    = 0;
        ManualROIspline = 0;
        pointNclickROI  = 0;
        ROIsize         = job.AutoROIchoice.AutoROI.ArrayROI;
    case 'ManualROI'
        autoROI         = 0;
        graphicalROI    = 1;
        ManualROIspline = 0;
        pointNclickROI  = 0;
    case 'ManualROIspline'
        autoROI         = 0;
        graphicalROI    = 0;
        ManualROIspline = 1;
        pointNclickROI  = 0;
    case 'pointNclickROI'
        autoROI         = 0;
        graphicalROI    = 0;
        ManualROIspline = 0;
        pointNclickROI  = 1;
        radius          = job.AutoROIchoice.pointNclickROI.ManualROIradius;
    case 'ManualEnterROI'
        autoROI         = 0;
        graphicalROI    = 0;
        ManualROIspline = 0;
        pointNclickROI  = 0;
    otherwise
        % Do nothing
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
            
            % Display anatomical image
            try
                vol = spm_vol(IOI.res.file_anat);
            catch
                disp('Could not find anatomical image');
                [t sts] = spm_select(1,'image','Select anatomical image','',dir_ioimat,'.*',1);
                IOI.res.file_anat = t;
                vol = spm_vol(IOI.res.file_anat);
            end
            vx = [1 1 1];
            [dir1 fil1] = fileparts(vol.fname);
            im_anat = spm_read_vols(vol);
            % Gray-scale colormap to enhance contrast
            if use_gray_contrast
                cmap = contrast(im_anat);
            end
            % Goto figures window
            spm_figure('GetWin', 'Graphics');
            spm_figure('Clear', 'Graphics');
            
            if isfield(job,'displayBrainmask')
                if job.displayBrainmask == 1
                    % Display only brain pixels mask
                    vol = spm_vol(IOI.fcIOS.mask.fname);
                    full_mask = logical(spm_read_vols(vol));
                else
                    % Display all the image
                    full_mask = ones(size(im_anat));
                end
            else
                full_mask = ones(size(im_anat));
            end
            
            if ~autoROI
                %Display images of changes from 10th to 90th percentile for all sessions
                try % <--- Needs a catch statement //EGC
                    nCols = ceil(sqrt(length(IOI.sess_res)));
                    nRows = ceil(length(IOI.sess_res) / nRows);
                    for i0=1:length(IOI.sess_res)
                        % Only open 1 figure for all the sessions
                        if i0 == 1,
                            hs = figure;
                            set(hs, 'Name', '10-90 percentile changes')
                        end
                        V = spm_vol(IOI.sess_res{i0}.fname_change_90_10{1}); %color green
                        tmp_image = spm_read_vols(V);
                        figure(hs);
                        subplot(nRows, nCols, i0);
                        imagesc(tmp_image); axis image
                        title(['Session ' int2str(i0) ': ratio of 90th to 10th percentile']);
                    end
                end
                oneMoreROI = 1;
                % Goto interactive window
                h2 = spm_figure('GetWin', 'Interactive');
                spm_input(['Subject ' int2str(SubjIdx) ' (' IOI.subj_name ')' ],'-1','d');
                linecount = 0;
                while oneMoreROI
                    figure(h2);
                    oneMoreROI = spm_input('Add an ROI?',2*linecount+2,'y/n');
                    if oneMoreROI == 'y', oneMoreROI = 1; else oneMoreROI = 0; end
                    if oneMoreROI
                        % h1 = figure('Position',[20 50 3*size(im_anat,1) 3*size(im_anat,2)]);
                        % Display anatomical image on SPM graphics window
                        spm_figure('GetWin', 'Graphics');
                        spm_figure('Clear', 'Graphics');
                        imagesc(im_anat .* full_mask);
                        if use_gray_contrast
                            colormap(cmap);
                        else
                            colormap(gray);
                        end
                        axis image;
                        if graphicalROI
                            % Specify polygonal region of interest (ROI)
                            title('Make ROI polygon, then double click in it to create ROI.');
                            mask = roipoly;
                        elseif ManualROIspline
                            % Start interactive ROI tool to choose spline
                            % ROI/seed
                            mask = ioi_roi_spline(im_anat,[],[],'Mark spline points, then right-click in it to create ROI/seed');
                        else
                            if pointNclickROI
                                % Specify center of circular ROI/seed with mouse
                                % point & click on the anatomical image
                                title('Click the center of circular ROI/seed')
                                % Circular seed setup
                                t = 0:pi/100:2*pi;
                                % Prompt user to point & click
                                p = ginput(1);
                                % Center of the seed coordinates
                                x0 = p(1); y0 = p(2);
                                % Parametric function for a circle
                                xi = radius * cos(t) + x0;
                                yi = radius * sin(t) + y0;
                                LineHandler = line(xi,yi,'LineWidth',3,'Color',[.8 0 0]);
                                % Create ROI/seed mask
                                mask = poly2mask(xi, yi, size(im_anat,1), size(im_anat,2));
                                % Save coordinates of seed for later display.
                                % NOTE: Row is 1st coordinate, Column is 2nd
                                IOI.res.ROI{index+1}.center = [y0 x0];
                                IOI.res.ROI{index+1}.radius = radius;
                            else
                                % Manual ROI coordinate entry
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
                                % Save coordinates of seed for later display
                                % NOTE: Row is 1st coordinate, Column is 2nd
                                IOI.res.ROI{index+1}.center = rc;
                                IOI.res.ROI{index+1}.radius = radius;
                            end
                        end
                        
                        mask = single(mask);
                        % Update ROI's display
                        full_mask = full_mask - mask;
                        index = index + 1; linecount = linecount + 1;
                        if job.select_names
                            figure(h2);
                            name = spm_input(['Enter name of ROI' int2str(index)],2*linecount+1,'s');
                        else
                            name = int2str(index);
                        end
                        IOI.res.ROI{index}.name = name;
                        if index < 10, str0 = '0'; else str0 = ''; end
                        str = [str0 int2str(index)];
                        % Save nifti files in ROI sub-folder //EGC
                        fname_mask = fullfile(newDir,[fil1 '_ROI_' str '.nii']);
                        IOI.res.ROI{index}.fname = fname_mask;
                        ioi_save_nifti(mask, fname_mask, vx);
                    end
                end
                % try close(h1); end
                % try close(h2); end
                try close(hs); end
                
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
                        % Save nifti files in ROI sub-folder //EGC
                        fname_mask = fullfile(newDir,[fil1 '_ROI_' str '_' int2str(i1) 'x' int2str(i2) '.nii']);
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

% EOF
