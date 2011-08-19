function ioi_shrink_run(job)
%
%shrinkage configuration
shrink_x = job.shrink_x;
shrink_y = job.shrink_y;
vx=[shrink_x shrink_y 1];
files1 = job.files1;
for i1=1:length(files1)
    try
        private_read_write_bin_image(files1{i1,:},vx);
    catch exception
        disp(exception.identifier)
        disp(exception.stack(1))
    end
end
end
%END


% Internal private function to read binary images, facilitates reading the
% code above.
function private_read_write_bin_image(fname,vx)
read_format = 'int16';
% Open first file to see which format is used (new or old)
% Describe image formats here to make code clearer
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
s1= I(3);
s2= I(1);
fclose(fidA);

% Now read all files with the dimensions determined above, variable
% image_total will contain the reconstituted images
nx=size(vx(1):vx(1):s1,2);
ny=size(vx(2):vx(2):s2,2);
fidA = fopen(fname);
image_part=int16((fread(fidA,read_format)));
fclose(fidA);
% Reshape to images
if ~isempty(noImageNouveauFormat) % new format
    sz_add = 3;
else
    sz_add = 0;
end
nframes = round(length(image_part)/(s1*s2+sz_add+4));
temp=reshape(image_part,(s1*s2+sz_add+4),nframes); 
header = reshape(temp(1:(4+sz_add),:),[4 nframes]);
image_total=reshape(temp((5+sz_add):end,:),[s1 s2 nframes]);
if vx(1) > 1 || vx(2) > 1
    image_total = ioi_imresize(double(image_total),0,nx,ny,vx(1),vx(2));
end
% for i1=1:nframes
% if vx(1) > 1 || vx(2) > 1
%     image_total(:,:,i1) = ioi_imresize(squeeze(double(image_total(:,:,i1))),0,nx,ny,vx(1),vx(2));
% end
% end
image_total = int16(image_total);
%add header
for i1=1:nframes
    header(1,i1) = ny;
    header(3,i1) = nx; 
end
image_total = reshape(image_total,[nx*ny nframes]);
image_total = cat(1,header,image_total);
fidA = fopen(fname,'w');
fwrite(fidA,image_total(:),read_format);
fclose(fidA);
end