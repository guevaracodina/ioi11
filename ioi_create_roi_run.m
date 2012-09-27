function out = ioi_create_roi_run(job)
use_gray_contrast = job.use_gray_contrast;
try
    SelectPreviousROI = job.SelectPreviousROI;
catch
    SelectPreviousROI = 0;
end
% Manual/automatic selection of ROIs/seeds (pointNclickROI ManualROIspline)
activMask = 0; %Boolean to use activation mask
fNamesROIchoice = fieldnames(job.AutoROIchoice);
autoROI         = 0;
graphicalROI    = 0;
ManualROIspline = 0;
pointNclickROI  = 0;
pointNclickROIsquare  = 0;
switch(fNamesROIchoice{1})
    case 'AutoROI'
        autoROI         = 1;
        ROIsize         = job.AutoROIchoice.AutoROI.ArrayROI;
    case 'ManualROI'
        graphicalROI    = 1;
    case 'ManualROIspline'
        ManualROIspline = 1;
    case 'pointNclickROI'
        pointNclickROI  = 1;
        radius          = job.AutoROIchoice.pointNclickROI.ManualROIradius;
    case 'pointNclickSquare'
        pointNclickROIsquare  = 1;
        ManualROIwidth  = job.AutoROIchoice.pointNclickSquare.ManualROIwidth;
        ManualROIheight = job.AutoROIchoice.pointNclickSquare.ManualROIheight;
    case 'ManualEnterROI'
        % Do nothing
    otherwise
        % Do nothing
end

for SubjIdx=1:length(job.IOImat)
    try
        %Load IOI.mat information
        [IOI IOImat dir_ioimat]= ioi_get_IOI(job,SubjIdx);
        
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
            
            if isfield(job.activMask_choice,'activMask')
                try
                    mask_image = job.activMask_choice.activMask.mask_image{1};
                    hM=hgload(mask_image);
                    close(hM);
                    activMask = 1;
                catch
                    disp('Could not mask by specified mask -- ROI creation failed')
                end
            end
            
            % Goto figures window
            if ~activMask
                spm_figure('GetWin', 'Graphics');
                spm_figure('Clear', 'Graphics');
            end
            if isfield(job,'displayBrainmask')
                if job.displayBrainmask == 1 && isfield(IOI,'fcIOS') && isfield(IOI.fcIOS,'mask') && isfield(IOI.fcIOS.mask,'fname')
                    % Display only brain pixels mask
                    vol = spm_vol(IOI.fcIOS.mask.fname);
                    full_mask = logical(spm_read_vols(vol));
                else
                    % Display all the image
                    full_mask = ones(size(im_anat));
                    if job.displayBrainmask
                        disp('Could not find brain mask')
                    end
                end
            else
                full_mask = ones(size(im_anat));
            end
            if SelectPreviousROI
                % Goto interactive window
                h2 = spm_figure('GetWin', 'Interactive');
                spm_input(['Subject ' int2str(SubjIdx) ' (' IOI.subj_name ')' ],'-1','d');
                SelPrevROI = spm_input('Select a previous list of ROIs?','+1','y/n');
                if SelPrevROI == 'y'
                    [tPrevROI stsPrevROI] = spm_select(1,'mat','Select IOI.mat structure containing information on desired ROIs','',dir_ioimat,'IOI.mat',1);
                    if stsPrevROI
                        IOI0 = IOI; %Store current IOI
                        try
                            load(tPrevROI);
                            IOI_withROIs = IOI;
                            IOI = IOI0;
                            try
                                if ~isfield(IOI.res,'ROI')
                                    IOI.ROIname = IOI_withROIs.ROIname;
                                    IOI.res.ROI = IOI_withROIs.res.ROI;
                                else
                                    IOI.ROIname = [IOI.ROIname; IOI_withROIs.ROIname];
                                    IOI.res.ROI = [IOI.res.ROI IOI_withROIs.res.ROI];
                                end
                                index = index+length(IOI.res.ROI);
                                clear IOI_withROIs
                            catch
                                IOI = IOI0;
                                disp('Specified IOI.mat structure does not contain valid ROI information')
                            end
                        catch
                            disp('Could not load IOI.mat structure containing desired ROIs')
                        end
                    end
                end
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
                        if ~activMask
                            imagesc(tmp_image); axis image
                        else
                            hM = hgload(mask_image);
                            movegui(hM,'center');
                        end
                        title(['Session ' int2str(i0) ': ratio of 90th to 10th percentile']);
                    end
                end
                oneMoreROI = 1;
                % Goto interactive window
                h2 = spm_figure('GetWin', 'Interactive');
                spm_input(['Subject ' int2str(SubjIdx) ' (' IOI.subj_name ')' ],'-1','d');
                linecount = 0;
                if activMask
                    if length(job.activMask_choice.activMask.mask_image) == 2
                        oxyOn = spm_input('Y = oxy, N = deoxy)',1,'y/n');
                        if oxyOn == 'y'
                            mask_image = job.activMask_choice.activMask.mask_image{1};
                        else
                            mask_image = job.activMask_choice.activMask.mask_image{2};
                        end
                    end
                end
                while oneMoreROI
                    figure(h2);
                    oneMoreROI = spm_input('Add an ROI?',2*linecount+2,'y/n');
                    if oneMoreROI == 'y', oneMoreROI = 1; else oneMoreROI = 0; end
                    if oneMoreROI
                        % h1 = figure('Position',[20 50 3*size(im_anat,1) 3*size(im_anat,2)]);
                        % Display anatomical image on SPM graphics window
                        if ~activMask
                            spm_figure('GetWin', 'Graphics');
                            spm_figure('Clear', 'Graphics');
                            imagesc(im_anat .* full_mask);
                            if use_gray_contrast
                                colormap(cmap);
                            else
                                colormap(gray);
                            end
                            axis image;
                        else
                            hM = hgload(mask_image);
                            movegui(hM,'center');
                        end
                        
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
                                if pointNclickROIsquare
                                    % Specify center of a square ROI/seed with mouse
                                    % point & click on the anatomical image
                                    title('Click the center of rectangular ROI/seed')
                                    % Square setup
                                    hR = imrect(gca,[0 0 ManualROIwidth ManualROIheight]);
                                    pR = wait(hR);
                                    xi = [pR(1) pR(1)+pR(3) pR(1)+pR(3) pR(1)];
                                    yi = [pR(2) pR(2)       pR(2)+pR(4) pR(2)+pR(4)];
                                    pRr = round(pR);
                                    disp(['Width = ' int2str(pRr(3)) ', Height = ' int2str(pRr(4)) ...
                                        ', Center(x) = ' int2str(round(pRr(1)+pRr(3)/2)) ...
                                        ', Center(y) = ' int2str(round(pRr(2)+pRr(4)/2))])
                                    % Create the mask
                                    mask = poly2mask(xi, yi, size(im_anat,1), size(im_anat,2));
                                  
                                    % Save coordinates of seed for later display.
                                    % NOTE: Row is 1st coordinate, Column is 2nd
                                    %IOI.res.ROI{index+1}.center = [y0 x0];
                                    %IOI.res.ROI{index+1}.radius = radius;
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
                        fname_mask = fullfile(dir_ioimat,[fil1 '_ROI_' str '_' name '.nii']);
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
                        fname_mask = fullfile(dir_ioimat,[fil1 '_ROI_' str '_' int2str(i1) 'x' int2str(i2) '.nii']);
                        IOI.res.ROI{index}.fname = fname_mask;
                        ioi_save_nifti(mask, fname_mask, vx);
                    end
                end
            end
            if activMask, try close(hM); end, end
            IOI.ROIname = {};
            for i0=1:length(IOI.res.ROI)
                if isfield(IOI.res.ROI{i0},'name')
                    IOI.ROIname = [IOI.ROIname; IOI.res.ROI{i0}.name];
                else
                    IOI.ROIname = [IOI.ROIname; ['ROI' gen_num_str(i0,3)]];
                end
            end
            
            IOI.res.ROIOK = 1;
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
