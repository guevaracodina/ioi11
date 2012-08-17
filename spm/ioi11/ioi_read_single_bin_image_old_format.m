% read anatomical binary image -- old format -- this only reads the first frame
function image_total = ioi_read_single_bin_image_old_format(fname)
% Open file containing a single volume
read_format = 'int16';
fidA = fopen(fname);
I = fread(fidA,2,read_format);
if all(I==[1;0]) % New format including image number
    noImageNouveauFormat=fread(fidA,1,read_format);
    I = fread(fidA,4,read_format);
else
    noImageNouveauFormat=[];
    I =[I; fread(fidA,2,read_format)];
end
% Image dimensions read from header (first few int16)
nx= I(3);
ny= I(1);
fclose(fidA);

fidA = fopen(fname);
image_total=int16((fread(fidA,read_format)));
fclose(fidA);
% Reshape to images
if ~isempty(noImageNouveauFormat) % new format
    temp(:,:)=squeeze(reshape(image_total,(nx*ny+3+4),round(length(image_total)/(nx*ny+3+4))));
else
    temp(:,:)=squeeze(reshape(image_total,(nx*ny+4),round(length(image_total)/(nx*ny+4))));
end
clear image_total
if ~isempty(noImageNouveauFormat)  % new format
    image_total=reshape(temp(8:end,:),[nx ny]);
else
    image_total=reshape(temp(5:end,1),[nx ny]);
end
image_total = single(image_total);