%read binary images
function [image_total fcount_out]= ioi_read_bin_image_old_format(fnames,...
    fcount,indices,vx,n_frames,first_file,last_file,im_count)
%Example of usage to read a single volume:
% [image_total fcount_out]= ioi_read_bin_image_old_format('JJ000000.bin',0,1:83,[1 1 1],83,1,20,0);

%If frames are missing between this fnames file and the next, by a convention
%defined here, they will be added at the end of this block of frames (not the next).
%fcount and fcount_out are the number of recorded frames before and after
%processing of the current block
read_format = 'int16';
% Open first file to see which format is used (new or old)
% Describe image formats here to make code clearer
%
fidA = fopen(fnames);
I = fread(fidA,2,read_format);
if all(I==[1;0]) % New format including image number
    noImageNouveauFormat=fread(fidA,1,read_format);
    I = fread(fidA,4,read_format);
else
    noImageNouveauFormat=[];
    I =[I; fread(fidA,2,read_format)];
end
% Image dimensions read from header (first few int16)
s1= I(3);
s2= I(1);
fclose(fidA);

% Now read all files with the dimensions determined above, variable
% image_total will contain the reconstituted images
nx=size(vx(1):vx(1):s1,2);
ny=size(vx(2):vx(2):s2,2);
fidA = fopen(fnames);
image_part=int16((fread(fidA,read_format)));
fclose(fidA);
% Reshape to images
if ~isempty(noImageNouveauFormat) % new format
    sz_add = 3;
else
    sz_add = 0;
end
image_part=reshape(image_part,(s1*s2+sz_add+4),round(length(image_part)/(s1*s2+sz_add+4)));
list_frames = fcount+(1:size(image_part,2)); %this list does not know about missing frames
fcount_out = list_frames(end);
last_missing_is_last_index = 0;
%first_missing_is_first_index = 0;
diff0 = diff(indices)~=1;
%missing_indices_locations=find(diff0); %location of missing frames in indices list
missing_indices=indices(diff0); %actual first missing frame numbers

%catch in case we are at the last list of frames
try
    missing_indices_select = missing_indices < indices(list_frames(end)+1) & ...
        missing_indices >= indices(list_frames(1)); %?
catch
    missing_indices_select = missing_indices <= indices(list_frames(end)) & ...
        missing_indices >= indices(list_frames(1));
    if indices(end) < n_frames
        indices(end+1) = n_frames+1;
        %redo previous steps
        diff0 = diff(indices)~=1;
        %missing_indices_locations=find(diff0); %location of missing frames in indices list
        missing_indices=indices(diff0);
        missing_indices_select = missing_indices < indices(list_frames(end)+1) & ...
            missing_indices >= indices(list_frames(1));
    end
end
diff1 = diff(indices);
diff_index = diff1(diff0); %gap size
diff_index = diff_index(missing_indices_select);
missing_indices = missing_indices(missing_indices_select);
%missing_indices_locations = missing_indices_locations(missing_indices_select);

if ~isempty(missing_indices) && (missing_indices(end) >= indices(list_frames(end))) %actual frames
    last_index = indices(list_frames(end)+1)-indices(list_frames(1));
    last_missing_is_last_index = 1;
else
    last_index = indices(list_frames(end))-indices(list_frames(1))+1;
end

temp = zeros(size(image_part,1),last_index);
temp(:,indices(list_frames)-indices(list_frames(1))+1) = image_part(:,list_frames-fcount);
% At each location fill the gaps naively
shft = indices(list_frames(1))-1;
for idx=1:length(missing_indices)
    b=missing_indices(idx)+1;
    e=b+diff_index(idx)-2;
    if ~last_missing_is_last_index || ~(idx == length(missing_indices))
        % Find gap width
        for n=b:e
            %left side: actual frames (with first index reset to 1)
            %right side:
            temp(:,n-shft)=0.5*(temp(:,b-1-shft)+temp(:,e+1-shft));
        end
    else
        %in case last_missing_is_last_index, then fill all missing with
        %last available image in this block.
        for n=b:e
            temp(:,n-shft)=temp(:,b-1-shft);
        end
    end
end
if first_file == 1
    %check whether first image is missing
    first_rec_frame = indices(list_frames(1));
    if first_rec_frame > 1
        %augment temp and last_index
        last_index = last_index+first_rec_frame-1;
        temp = [repmat(temp(:,1),[1 first_rec_frame-1])  temp];
    end
end
output_debugging_info = 0;
if output_debugging_info
    %Debugging info
    format compact
    disp(['Fr processed earlier: ' int2str(fcount) ', Frames now: ' int2str(fcount_out)]);
    disp(['Fr available: ' int2str(fcount_out-fcount) ', Images added: ' int2str(last_index)]);
    disp('missing_indices_locations:');
    disp(missing_indices_locations');
    disp('missing_indices:');
    disp(missing_indices');
    disp([int2str(list_frames(1)) ': list_frames(1), ' int2str(indices(list_frames(1))) ': indices(list_frames(1))']);
    disp([int2str(list_frames(end)) ': list_frames(end), ' int2str(indices(list_frames(end))) ': indices(list_frames(end))']);
    disp('**********');
end
%in the hopefully rare instances when there is still a mismatch between
%n_frames and the total number of interpolated images, enforce the n_frames
%value by adding or subtracting
if last_file
    im_diff = n_frames-im_count-last_index;
    if im_diff>0
        %missing some frames, add some
        temp = [temp repmat(temp(:,end),[1 im_diff])];
        last_index = last_index + im_diff;
        disp(['Warning: missing ' int2str(im_diff) ' images at the end, filled in to get n_frames']);
    else
        if im_diff<0
            %too many frames, drop some
            temp = temp(:,1:(end+im_diff));
            last_index = last_index + im_diff;
            disp(['Warning (minor): too many images ' int2str(abs(im_diff)) ' at the end, removed some to get n_frames']);
        end
    end
end
%Put here in 4D format
image_total=reshape(temp((5+sz_add):end,:),[s1 s2 1 last_index]);
if vx(1) > 1 || vx(2) > 1
    image_total = ioi_imresize(image_total,0,nx,ny,vx(1),vx(2));
end
image_total = single(image_total);