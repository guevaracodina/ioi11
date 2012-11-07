function image_out = ioi_imresize(image,mode,nx,ny,vx1,vx2)
% Resize 4-D data slice by slice.
%_______________________________________________________________________________
% Copyright (C) 2012 LIOM Laboratoire d'Imagerie Optique et Moleculaire
%                    Ecole Polytechnique de Montreal
%_______________________________________________________________________________
% mode = 0 is faster
% Preallocating memory for the output can improve speed in mode 1. //EGC
image_out = zeros(nx, ny, size(image,3), size(image,4));
if mode %3 times slower than naive shrink below
    for i3=1:size(image,3)
        for i4=1:size(image,4)
            % image_out(:,:,i3,i4) = imresize(squeeze(image(:,:,i3,i4)),[nx ny],'bicubic');
            image_out(:,:,i3,i4) = ioi_MYimresize(squeeze(image(:,:,i3,i4)),[nx ny],'bicubic');
        end
    end
else
    % image_out = zeros(nx,ny,size(image,3),size(image,4)); 
    % Preallocation is already done //EGC
    for i1=1:vx1
        for i2=1:vx2
            image_out = image_out + ...
                image(i1:vx1:(end-vx1+i1),i2:vx2:(end-vx2+i2),:,:);
        end
    end    
    image_out = image_out/(vx1*vx2);
end

