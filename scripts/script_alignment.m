%% script_alignment
clc;
s1 = 1;
c1 = 5;
r1 = 1;
treatmentString = 'NC';

%% Load manually aligned images with ImageJ plugin TurboReg
clear imDataArray
topDir = 'D:\Edgar\Data\IOS_Carotid_Res\alignment';
currentDir = sprintf('%s_R%02dC%02d', treatmentString, r1, c1);
[dirAlignment, sts] = cfg_getfile(1,'dir','Select folder',{fullfile(topDir,currentDir)}, topDir, '.*');
[images2align, sts] = cfg_getfile(Inf,'image','Select folder',[], dirAlignment{1}, '.*');
V = spm_vol(images2align);
for iVols = 1:numel(V),
    [imData, imXYZ] = spm_read_vols(V{iVols});
    % Set aligned IMG/HDR images to the orientation of source image
    if ndims(imData) ~= 2
    	imData = squeeze(imData(:,:,1));
        imData = fliplr(imData);
    end
    imDataArray(:,:,iVols) = imData;
end

imAvg = median(imDataArray, 3);

%% Display average image
h = figure; imagesc(imAvg, [-1 1]); axis image; colorbar
figFolder = 'D:\Edgar\Documents\Dropbox\Docs\Carotid\Figures\aligned';
% Save as PNG
print(h, '-dpng', fullfile(figFolder,[currentDir '_avg']), '-r150');
