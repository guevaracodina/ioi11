function ioi_write_nifti(varargin)
% Writes a 3D or 4D multivolume image in NIFTI-1 format (.nii extension)
% SYNTAX:
% ioi_write_nifti(img,voxel_size,origin,datatype,description,file)
% INPUTS:
% img:          Usually, img is a 3D matrix [x y z], or a 4D
%               matrix with time series [x y z t]. NIfTI allows a maximum of 7D
%               matrix, however spm_lot curently support up to 4D data.
%  
% voxel_size:   [OPTIONAL]:	Voxel size in millimeter for each
%  				dimension. Default is [1 1 1 1].
%  
% origin:       [OPTIONAL]	The AC (anterior comissure) origin. Default is at
%               the middle of the central slice.
%  
% datatype:     [OPTIONAL]	Storage data type:
%               2 - uint8,  4 - int16,  8 - int32,  16 - float32,
%               32 - complex64,  64 - float64,  128 - RGB24,
%               256 - int8,  512 - uint16,  768 - uint32, 
%               1792 - complex128
%               Default will use the data type of 'img' matrix
%  
% description:  [OPTIONAL]	Any text you like as a description of data. 
%               Default is 'Created with spm_lot.'.
% file:         [OPTIONAL] Location and filename of .nii file, optional, but
%               strongly recommended
% OUTPUTS:
% none
%_______________________________________________________________________________
% Copyright (C) 2010 Laboratoire d'Imagerie Optique et Moleculaire
% Edgar Guevara
% 2010/04/13

%     FIELD         NOTES
%     -----------------------------------------------------
%     sizeof_hdr    must be 348
%     -----------------------------------------------------
%     dim           dim[0] and dim[1] are always required; 
%                   dim[2] is required for 2-D volumes, 
%                   dim[3] for 3-D volumes, etc.
%     -----------------------------------------------------
%     datatype      needed to specify type of image data
%     -----------------------------------------------------
%     bitpix        should correspond correctly to datatype
%     -----------------------------------------------------
%     pixdim        with the exception of pixdim[0] (which 
%                   is required when qform_code != 0), 
%                   pixdim[n] is required when dim[n] is 
%                   required
%     -----------------------------------------------------
%     vox_offset    required for an "n+1" header
%     -----------------------------------------------------
%     magic         must be "ni1\0" or "n+1\0"
%     -----------------------------------------------------

img = varargin{1};
dims = size(img);
dims = [dims ones(1,8)];
dims = dims(1:4);

voxel_size = ones([1 4]);
origin = zeros([1 3]);
descrip = 'Created with ioi_write_niftii.';

switch class(img)
    case 'uint8'
        datatype = 2;
    case 'int16'
        datatype = 4;
    case 'int32'
        datatype = 8;
    case 'single'
        datatype = 16;
    case 'double'
        datatype = 64;
    case 'int8'
        datatype = 256;
    case 'uint16'
        datatype = 512;
    case 'uint32'
        datatype = 768;
    otherwise
        error('Datatype is not supported.');
end

if nargin > 1 && ~isempty(varargin{2})
    voxel_size(1:3) = double(varargin{2});
    origin = -voxel_size(1:3).*dims(1:3)/2;
end

if nargin > 2 && ~isempty(varargin{3})
    origin(1:3) = double(varargin{3});
end

if nargin > 3 && ~isempty(varargin{4})
    datatype = double(varargin{4});
end

if nargin > 4 && ~isempty(varargin{5})
    descrip = varargin{5};
end

if nargin > 5 && ~isempty(varargin{6})
    file = varargin{6};
    [pathName fileName fileExt] = fileparts(file);
    if isempty(pathName)
        pathName =  @(val)ioi_get_defaults('dirTemplate', val{:});
    end
    if isempty(fileName)
        fileName = 'defaultNIFTI';
    end
    if isempty(fileExt)
        fileExt = '.nii';
    end
else
    pathName =  @(val)ioi_get_defaults('dirTemplate', val{:});
    fileName = 'defaultNIFTI';
    fileExt = '.nii';
end

if datatype == 128
    if ndims(img) > 5
        error('ioi_write_nifti only allows a maximum of 4 Dimension matrix.');
    end
%     dims(1) = dims(1)-1;
%     dims(5:8) = [dims(6:8) 1];
else
    if ndims(img) > 4
        error('ioi_write_nifti only allows a maximum of 4 Dimension matrix.');
    end
end
   
%currentDir = cd;                                % Current directory
%cd(pathName);                                   % Directory to save data
% Initialize file array object
dat         = file_array;                       % File array object
dat.fname   = fullfile(pathName,[fileName fileExt]);               % Filename
dat.dim     = dims;                             % Dimensions (in voxels)
switch datatype
    case 2
        dat.dtype = 'UINT8-BE';                 % Data type (Big Endian)
        img = uint8(img);                       % Data type conversion
    case 4
        dat.dtype = 'INT16-BE';                 % Data type (Big Endian)
        img = int16(img);                       % Data type conversion
    case 8
        dat.dtype = 'INT32-BE';                 % Data type (Big Endian)
        img = int32(img);                       % Data type conversion
    case 16
        dat.dtype = 'FLOAT32-BE';               % Data type (Big Endian)
        img = single(img);                      % Data type conversion
    case 64
        dat.dtype = 'FLOAT64-BE';               % Data type (Big Endian)
        img = double(img);                      % Data type conversion
    case 256
        dat.dtype = 'INT8-BE';                  % Data type (Big Endian)
        img = int8(img);                        % Data type conversion
    case 512
        dat.dtype = 'UINT16-BE';                % Data type (Big Endian)
        img = uint16(img);                      % Data type conversion
    case 768
        dat.dtype = 'UINT32-BE';                % Data type (Big Endian)
        img = uint32(img);                      % Data type conversion
    otherwise
        error('Datatype is not supported.');
end
dat.offset  = ceil(348/8)*8;                    % Voxel offset
% In a .nii file, the vox_offset field value is interpreted as the start
% location of the image data bytes in that file. 352 means data starts
% immediately after the header.
% disp(dat)                                       % Show file structure

% Create an empty NIFTI structure
N = nifti;                                      % N will contain the structure
% fieldnames(N)                                   % Dump fieldnames

% Creating all the NIFTI header stuff
% More information on:
% http://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/index_html
% ==============================================================================
switch length(dims)
    case 1,
    dat(:) = img;                               % Filling up with data
    case 2,
    dat(:,:) = img;                             % Filling up with data
    case 3,
    dat(:,:,:) = img;                           % Filling up with data
    case 4,
    dat(:,:,:,:) = img;                         % Filling up with data
    otherwise
        disp('Dimensions not supported in ioi_write_nifti')
end
N.dat = dat;                                    % 3D volume information
% NIFTI-1 uses a right-handed system for volume orientation and location in
% space. +x = Right  +y = Anterior  +z = Superior.
%  METHOD 3 (used when sform_code > 0):
% -----------------------------------
% The (x,y,z) coordinates are given by a general affine transformation
% of the (i,j,k) indexes:
% 
% x = srow_x[0] * i + srow_x[1] * j + srow_x[2] * k + srow_x[3] 
% y = srow_y[0] * i + srow_y[1] * j + srow_y[2] * k + srow_y[3] 
% z = srow_z[0] * i + srow_z[1] * j + srow_z[2] * k + srow_z[3]
% 
% i.e.
% | x |       | i |
% | y | =   M | j |
% | z |       | k |
% | 0 |       | 0 |
%
% Columns 1 to 3 represent a scale factor (voxel size), whereas column 4 is an
% offset (origin)
% The srow_* vectors are in the NIFTI_1 header. 
% Note that no use is made of pixdim[] in this method.
N.mat = [...
    % 1st row affine transform (x).
    voxel_size(1)   0               0               origin(1);
    % 2nd row affine transform (y).
    0               voxel_size(2)   0               origin(2);
    % 3rd row affine transform (z).
    0               0               voxel_size(3)   origin(3);
    % 4th row must always be [0 0 0 1]
    0               0               0               1];

% Fill up description field
N.descrip = descrip;
% ANOTHER AFFINE MATRIX??? //EGC
N.mat0 = N.mat;
%cd(pathName)
create(N);                                      % Writes hdr info

%disp(['NIFTI-1 image saved as: ' pathName filesep fileName fileExt])
%cd(currentDir)