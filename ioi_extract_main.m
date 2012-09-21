function [IOI ROI] = ioi_extract_main(IOI,ROI,job,d,d3,d4,c1,s1,colorOK,mask,Amask)
[all_ROIs selected_ROIs] = ioi_get_ROIs(job);
msg_ColorNotOK = 1;

nROI = 1:length(IOI.res.ROI); % All the ROIs
if isfield(job,'extractingBrainMask')
    if job.extractingBrainMask
        nROI = 1; % Only 1 brain mask
    end
end
if ~isempty(Amask)
    mask = mask .* Amask;
end

for r1 = nROI,
    if all_ROIs || sum(r1==selected_ROIs)
        tmp_mask_done = 0;
        for i3=1:d3
            for i4=1:d4
                %extracted data
                %tmp_d = squeeze(d(:,:,i3,i4));
                try tmp_d = d(:,:,i3,i4); end
                %just take mean over mask for now
                
                if ~isfield(IOI.color,'contrast') || (isfield(IOI.color,'contrast') && ~(IOI.color.eng(c1)==IOI.color.contrast))
                    try
                        e = mean(tmp_d(mask{r1}));
                    catch
                        if msg_ColorNotOK
                            msg = ['Problem extracting for color ' int2str(c1) ', session ' int2str(s1) ...
                                ',region ' int2str(r1) ': size mask= ' int2str(size(mask{r1},1)) 'x' ...
                                int2str(size(mask{r1},2)) ', but size image= ' int2str(size(tmp_d,1)) 'x' ...
                                int2str(size(tmp_d,2))];
                            IOI = disp_msg(IOI,msg);
                            msg_ColorNotOK = 0;
                        end
                        if colorOK
                            try
                                %try to resize mask - but only attempt to do it once
                                if ~tmp_mask_done
                                    % tmp_mask = imresize(mask{r1},size(tmp_d));
                                    % ioi_MYimresize works with no image
                                    % processing toolbox //EGC
                                    tmp_mask = ioi_MYimresize(mask{r1},size(tmp_d));
                                    tmp_mask_done = 1;
                                end
                                e = mean(tmp_d(tmp_mask));
                            catch
                                msg = ['Unable to extract color ' int2str(c1) ', session ' int2str(s1)];
                                IOI = disp_msg(IOI,msg);
                                colorOK = 0;
                            end
                        end
                    end
                else
                    %contrast images will be smaller and need to be resized
                    % tmask = imresize(mask{r1},[d1 d2]);
                    % ioi_MYimresize works with no image processing toolbox
                    % //EGC
                    tmask = ioi_MYimresize(mask{r1},[d1 d2]);
                    e = mean(tmp_d(tmask));
                end
                if colorOK
                    ROI{r1}{s1,c1} = [ROI{r1}{s1,c1} e];
                end
            end
        end
    end
end

% EOF

