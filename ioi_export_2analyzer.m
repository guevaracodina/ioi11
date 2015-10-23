function OK = ioi_export_2analyzer (pathNameParent, pathNamePhysio)
% Export IOI data in analyzer format
% SYNTAX:
%  OK = ioi_export_2analyzer (pathNameParent, pathNamePhysio)
% INPUTS:
% pathNameParent	Parent directory with IOImat
% pathNamePhysio    Directory with physiological data (EKG, respiration, etc.)
% OUTPUTS:
% OK                Flag indicating a succesful execution
%_______________________________________________________________________________
% Copyright (C) 2012 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________
try
    %% Loading data after filtering and downsampling
    
    addpath(genpath('C:\spm8\toolbox\nirs10')); % After running spm8 shortcut
    % pathNameParent = 'C:\Edgar\Data\IOIResults\j03\ROItest';
    IOImat = fullfile(pathNameParent, 'IOI.mat');
    load(IOImat)
    if IOI.fcIOS.filtNdown.filtNdownOK
        % Load filtered & downsampled image time-course
        vol = spm_vol(IOI.fcIOS.filtNdown.fnameWholeImage{1});
        y = spm_read_vols(vol);
    end
    % Get image dimensions
    [nX, nY, nZ, nT] = size(y);
    % Reshape as 2-D matrix
    y = reshape(y,[nX*nY nT]);
    % Intrinsic Imaging frequency (after downsampling to 1 Hz)
    Fs = IOI.fcIOS.filtNdown.fs;
    Ts = 1/Fs;
    
    %% Load EKG & Respiration Data
    physio = ioi_Bin2Mat(pathNamePhysio);
    resp_freq = 250;        % sampling rate (respiration)
    resp_Ts = 1/resp_freq;  % sampling time (respiration) interval
    t = 0:resp_Ts:(numel(physio.Resp)-1)*resp_Ts;
    
    %% Check if filter coefficients are saved to a file
    if ~exist(fullfile(fileparts(mfilename('fullpath')), 'resp_filt_coeff.mat'), 'file')
        %%  Optimal Minimum Order Designs (BPF from 0.01 to 30 Hz)
        % Fp1 — frequency at the edge of the start of the pass band. Specified in
        % normalized frequency units. Also called Fpass1. }
        Fp1 = 0.03*2/resp_freq;
        % Fp2 — frequency at the edge of the end of the pass band. Specified in
        % normalized frequency units. Also called Fpass2.
        Fp2 = 28*2/resp_freq;
        % Fst1 — frequency at the edge of the start of the first stop band.
        % Specified in normalized frequency units. Also called Fstop1.
        Fst1 = 0.01*2/resp_freq;
        % Fst2 — frequency at the edge of the start of the second stop band.
        % Specified in normalized frequency units. Also called Fstop2.
        Fst2 = 30*2/resp_freq;
        % Ap — amount of ripple allowed in the pass band. Also called Apass.
        Ap = 1;
        % Ast1 — attenuation in the first stop band in decibels (the default
        % units). Also called Astop1.
        Ast1 = 30;
        Ast2 = Ast1;
        Hf = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2);
        % Equiripple designs result in the lowpass
        % filter with the smallest possible order to meet a set of specifications.
        Hd = design(Hf,'equiripple','systemobject',true);
        measure(Hd)
        hfvt = fvtool(Hd,'Color','White');
        legend(hfvt,'Equiripple design')
        
        % -------------------- Cosmetic changes -----------------------------------
        % My Red Color
        myColor = [204 0 0]/255;
        % set(gca,'FontSize', 12)
        legend({'Equiripple design'}, 'FontSize', 14)
        figureHandle = gcf;
        %# make all text in the figure to size 14 and bold
        set(findall(figureHandle,'type','text'),'fontSize',14)
        gcachild = get(gca,'Children');
        set(gca, 'Xcolor', 'k');
        set(gca, 'Ycolor', 'k');
        set(gcachild, 'Linewidth', 2);
        % Change plot marker
        hchildren = get(hfvt,'children');
        haxes = hchildren(strcmpi(get(hchildren,'type'),'axes'));
        hline = get(haxes,'children');
        set(hline(1), 'Color', [0 0 204]/255, 'LineWidth',2);
        set(hline(2), 'Color',  myColor, 'LineWidth',2);
        title(IOI.subj_name);
        
        % Get filter coefficients for the numerator
        b =  get(Hd,'Numerator');   % N = 75
        a = 1;  % The denominator of FIR filters is, by definition, equal to 1
    else
        % Retrieve filter coefficients saved in file
        load(fullfile(fileparts(mfilename('fullpath')), 'resp_filt_coeff.mat'))
    end
    %% Apply filter to the signal
    respFilt = filtfilt(b, a, physio.Resp);
    
    %% Plot filter and unfiltered signals
    figure;
    plot(t,physio.Resp)
    hold on
    plot(t,respFilt)
    xlim([115 118]);
    legend({'Raw' 'Filtered'}, 'FontSize', 14)
    title(IOI.subj_name);
    
    %% Compute resampling ratio
    % Free up memory
    clear b physio
%     pack
    % resampling at p/q times the original sampling rate
    q = round(resp_freq);
    p = round(Fs);
    % Cut data to experiment length
    nT = min(size(y,2)*resp_freq, numel(t));
    respFilt = respFilt(1:nT);
    t = t(1:nT);
    
    %% Memory-mapping
    tic
    clear im_obj
    fmem_name=fullfile(pathNameParent, [IOI.subj_name '_resampled.dat']);
    fid = fopen(fmem_name,'w');
    nPix = size(y,1)+1;             % Number of pixel time-traces
    fwrite(fid, zeros([nPix nT]), 'double');
    fclose(fid);
    im_obj = memmapfile(fmem_name,...
        'Format',{'double' [nPix nT] 'Frames'}, ...
        'writable',true);
    toc
    
    %% Interpolate & Resample
    % [ECGresamp,b] = resample(ECGfilt,p,q);
    tic
    ioi_text_waitbar(0, 'Please wait...');
    % Concatenate data, first channel is respiration data
    im_obj.Data.Frames(1,:) = respFilt;
    % [yResamp,b] = resample(y,q,p);
    
    for iPix = 2:nPix
        im_obj.Data.Frames(iPix,:) = resample(y(iPix-1,:),q,p);
        % Update progress bar
        ioi_text_waitbar(iPix/nPix, sprintf('Processing event %d from %d', iPix, nPix));
    end
    ioi_text_waitbar('Clear');
    toc
    
    %% Write data to be imported in Analyzer2
    DataFileNIR         = fullfile(pathNameParent,[IOI.subj_name '_resp_noCB.nir']);
    fwrite_NIR(DataFileNIR, im_obj.Data.Frames);
    % Delete memory-mapped file
    clear im_obj
    delete(fmem_name);
    
    %% .vhdr parameters
    DataFile            = fullfile('D:\Edgar\IOI_raw\',[IOI.subj_name '_resp_noCB.nir']);
    OutputFile          = fullfile(pathNameParent,[IOI.subj_name '_resp_noCB.vhdr']);
    CreatingFunction    = mfilename;
    ChannelResolution   = '1';
    ChannelUnits        = 'µV';
    ChannelLabels       = cellstr(int2str((1:nPix)'));
    ChannelLabels{1}    = 'Resp';               % 1st channel is respiration
    SamplingInterval    = round(1e6/resp_freq); % in µs
    dataPoints          = nT;
    % Create .vhdr file
    ioi_boxy_write_header(OutputFile,DataFile,...
    CreatingFunction,ChannelResolution,ChannelUnits,...
    ChannelLabels,SamplingInterval,dataPoints)
    OK = true;
catch exception
    disp(exception.identifier)
    disp(exception.stack(1))
    OK = false;
end

% EOF