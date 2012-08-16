%% LPF script
vol = spm_vol('D:\Edgar\Data\IOS_Res\12_06_20,CO04\12_06_20,CO04_anat.nii');
im_anat = spm_read_vols(vol);
K.radius = 5;
K.k1 = size(im_anat,1);
K.k2 = size(im_anat,2);
nTimes = 1000;

%% SPM spatial LPF
t = 0;
K = ioi_spatial_LPF('set', K);
for iTimes=1:nTimes,
    tic
    im_filt_spm = ioi_spatial_LPF('lpf', K, im_anat);
    t = t + toc;
end
t = t / nTimes;
fprintf('Average SPM time: %f sec\n',t)


%% Native Matlab LPF
t = 0;
sigma   = K.radius/2;
H = fspecial('gaussian',2*[K.radius K.radius]+1, sigma);

for iTimes=1:nTimes,
    tic
    im_filt = imfilter(im_anat, H, 'replicate');
    t = t + toc;
end
t = t / nTimes;
fprintf('Average MATLAB time: %f sec\n',t)

%% Display results
close all; figure;
subplot(121); imagesc(im_filt_spm); axis image; colormap (gray); title('SPM LPF')
subplot(122); imagesc(im_filt); axis image; colormap (gray); title('Matlab LPF')
figure; imagesc(im_filt_spm - im_filt); axis image; colormap (gray); title('Difference')


