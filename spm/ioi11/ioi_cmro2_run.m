function out = ioi_cmro2_run(job)
% Computes CMRO2 from HbO, HbR and blood flow data.
% M. Jones, J. Berwick, D. Johnston, and J. Mayhew, “Concurrent optical imaging
% spectroscopy and laser-Doppler flowmetry: the relationship between blood flow,
% oxygenation, and volume in rodent barrel cortex,” Neuroimage, vol. 13, no. 6
% Pt 1, pp. 1002–1015, Jun. 2001.
%_______________________________________________________________________________
% Copyright (C) 2012 LIOM Laboratoire d'Imagerie Optique et Moleculaire
%                    Ecole Polytechnique de Montreal
%_______________________________________________________________________________

% ------------------------------------------------------------------------------
% REMOVE AFTER FINISHING THE FUNCTION //EGC
% ------------------------------------------------------------------------------
% fprintf('Work in progress...\nEGC\n')
% out.IOImat = job.IOImat;
% return
% ------------------------------------------------------------------------------

% Select a subset of sessions
[all_sessions selected_sessions] = ioi_get_sessions(job);
% String that identifies CMRO2
str_CMRO2 = 'M';
tmp_str_CMRO2 = ['_' str_CMRO2 '_'];
% String that identifies flow
str_flow = 'F';
tmp_str_flow = ['_' str_flow '_'];
% String that identifies oxy-hemoglobin
str_HbO = 'O';
% String that identifies deoxy-hemoglobin
str_HbR = 'D';
% Vascular weighting constants
gammaT = job.constants.gammaT;
gammaR = job.constants.gammaR;
% Filter order
filterOrder = job.bpf.bpf_On.bpf_order;
% Band-pass cut-off frequencies
BPFfreq = job.bpf.bpf_On.bpf_freq;
% Filter type
fType = job.bpf.bpf_On.bpf_type;
% Passband/Stopband ripple in dB
Rp_Rs = [job.bpf.bpf_On.bpf_Rp job.bpf.bpf_On.bpf_Rs];

% Loop over subjects
for SubjIdx = 1:length(job.IOImat)
    try
        tic
        % Load IOI.mat information
        [IOI IOImat dir_ioimat]= ioi_get_IOI(job,SubjIdx);
        % Check if CMRO2 has not been computed or the job is redone,
        if (~isfield(IOI.res,'cmro2OK') || job.force_redo)
            %Check if flow and concentrations have already been computed.
            if (isfield(IOI.res,'flowOK') && isfield(IOI.res,'concOK'))
                % Original Sampling Frequency (5 Hz per color, data is sampled at 20
                % Hz for 4 colors RGYL)
                fs = 1/IOI.dev.TR;
                % Update color structure in IOI matrix
                IOI.color.CMRO2 = str_CMRO2;
                % Update string of colors
                if ~(IOI.color.eng==str_CMRO2)
                    IOI.color.eng = [IOI.color.eng str_CMRO2];
                end
                % Band-pass filter configuration. Computes zpk only once.
                [z, p, k] = temporalBPFconfig(fType, fs, BPFfreq, filterOrder, Rp_Rs);
                if isfield(IOI,'fcIOS')
                    if isfield(IOI.fcIOS,'mask')
                        % Check if a brain mask has already been selected
                        if IOI.fcIOS.mask.maskOK
                            % Read brain mask file
                            vol = spm_vol(IOI.fcIOS.mask.fname);
                            brainMask = logical(spm_read_vols(vol));
                        else
                            brainMask = [];
                        end
                    else
                        brainMask = [];
                    end
                else
                    brainMask = [];
                end
                % Check that sessions have been processed
                if isfield(IOI,'sess_res')
                    % Loop over sessions
                    for s1 = 1:length(IOI.sess_res)
                        if all_sessions || sum(s1==selected_sessions)
                            %tic
                            % Check available colors and files
                            [filesOK fname_list_flow fname_list_HbO fname_list_HbR] = local_check_NIfTI_files(IOI, s1, str_flow, str_HbO, str_HbR);
                            if filesOK
                                % Initialize cell with CMRO2 file names
                                fname_new_list = {};
                                % Initialize progress bar
                                spm_progress_bar('Init', length(fname_list_flow), sprintf('CMRO2 computation, %s, S%02d\n',IOI.subj_name,s1), 'Files');
                                % Loop over data files
                                for f1 = 1:length(fname_list_flow)
                                    % Save - substitute 'F' for 'M' in file name
                                    fname_new = regexprep(fname_list_flow{f1}, tmp_str_flow , tmp_str_CMRO2);
                                    if job.MemoryManagementMenu
                                        if ~(job.keepFiles && exist(fname_new,'file'))
                                            % Load all images at once
                                            [imagesFlow nx(1) ny(1) nt(1)] = local_read_NIfTI(fname_list_flow{f1});
                                            [imagesHbO nx(2) ny(2) nt(2)]  = local_read_NIfTI(fname_list_HbO{f1});
                                            [imagesHbR nx(3) ny(3) nt(3)]  = local_read_NIfTI(fname_list_HbR{f1});
                                            % If HbO, HbR & flow data have the same size
                                            if all(nx == nx(1)) && all(ny == ny(1)) && all(nt == nt(1))
                                                if ~isempty(brainMask)
                                                    % Test if there was shrinkage
                                                    if size(brainMask,1)~= nx(1) || size(brainMask,2)~= ny(1)
                                                        brainMask = ioi_MYimresize(brainMask, [nx(1) ny(1)]);
                                                    end
                                                else
                                                    brainMask = [];
                                                end
                                                %% Data filtering
                                                imagesFlow = local_filter_time_course(imagesFlow, z, p, k, s1, brainMask);
                                                imagesHbO  = local_filter_time_course(imagesHbO, z, p, k, s1, brainMask);
                                                imagesHbR  = local_filter_time_course(imagesHbR, z, p, k, s1, brainMask);
                                                %% CMRO2 computation
                                                imagesCMRO2 = ioi_cmro2_compute(imagesFlow, imagesHbO, imagesHbR, gammaT, gammaR, IOI.conc.baseline_hbt, IOI.conc.baseline_hbr);
                                            else
                                                error('ioi_cmro2_run:different_Sizes_Flow_HbO_HbR', 'Flow and Hb data have different sizes.');
                                            end
                                        end
                                    else
                                        % Load images one by one
                                        % TO DO...
                                    end
                                    % Reshrinks image if necessary
                                    vx = local_check_shrinkage(IOI);
                                    % Append new file name to file name list
                                    fname_new_list = [fname_new_list; fname_new];
                                    if ~(job.keepFiles && exist(fname_new,'file'))
                                        % Save as NIfTI file
                                        ioi_save_nifti(imagesCMRO2, fname_new, vx);
                                    end
                                    fprintf('%d of %d files done. %s saved.\n',f1, length(fname_list_flow), fname_new);
                                    % Update progress bar
                                    spm_progress_bar('Set', f1);
                                end % Loop over data files (NIfTI)
                                % Clear progress bar
                                spm_progress_bar('Clear');
                                IOI.sess_res{s1}.fname{IOI.color.eng==str_CMRO2} = fname_new_list;
                                disp(['CMRO2 calculation for session ' int2str(s1) ' complete']);
                            else
                                fprintf('Skipped CMRO2 computation.\nVerify if flow and concentrations are computed for session %d\n',s1)
                                % error('ioi_cmro2_run:dataNotComputed', 'Needs flow & Hb concentrations computed.');
                            end % F,O,D files OK
                        end
                    end % Loop over sessions
                end
                IOI.res.cmro2OK = 1;
                save(IOImat,'IOI');
            else
                fprintf('Skipped CMRO2 computation.\nVerify if flow and concentrations are computed for %s\n',IOI.subj_name)
                % error('ioi_cmro2_run:dataNotComputed', 'Needs flow & Hb concentrations computed.');
            end
        end
        out.IOImat{SubjIdx} = IOImat;
        disp(['Elapsed time: ' datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS')]);
        disp(['Subject ' int2str(SubjIdx) ' (' IOI.subj_name ')' ' complete']);
    catch exception
        disp(exception.identifier)
        disp(exception.stack(1))
        out.IOImat{SubjIdx} = job.IOImat{SubjIdx};
    end % End try
end % Loop over subjects
end % Main function

function [filesOK fname_list_flow fname_list_HbO fname_list_HbR] = local_check_NIfTI_files(IOI, s1, str_flow, str_HbO, str_HbR)
% Checks that F,O,D are available and NIfTI files exist for the current session.
filesOK = false;
fname_list_flow = [];
fname_list_HbO = [];
fname_list_HbR = [];
% check that flow is available
Fidx = find(IOI.color.eng == str_flow); %index for flow
% check that HbO is available
Oidx = find(IOI.color.eng == str_HbO); %index for HbO
% check that HbR is available
Didx = find(IOI.color.eng == str_HbR); %index for HbR
% Check that flow, HbO & HbR are saved as NIfTI files
if length(IOI.sess_res{s1}.fname) >= Fidx &&...
        length(IOI.sess_res{s1}.fname) >= Oidx &&...
        length(IOI.sess_res{s1}.fname) >= Didx
    fname_list_flow     = IOI.sess_res{s1}.fname{Fidx};
    fname_list_HbO      = IOI.sess_res{s1}.fname{Oidx};
    fname_list_HbR      = IOI.sess_res{s1}.fname{Didx};
    % Check same number of flow, HbO & HbR files
    if numel(fname_list_flow) ==  numel(fname_list_HbO) &&...
            numel(fname_list_flow) ==  numel(fname_list_HbR) &&...
            ~isempty(fname_list_flow) && ~isempty(fname_list_HbO) && ~isempty(fname_list_HbR)
        filesOK = true;
    end
end
end % local_check_NIfTI_files

function [dataOut nx ny nt] = local_read_NIfTI(fname)
% Loads data from NIfTI files
vol = spm_vol(fname);
nx = vol(1).dim(1);
ny = vol(1).dim(2);
% NOTE: nt is not necessarily the largest dimension of vol //EGC
if length(vol) == 1
    nt = vol(1).dim(3);
else
    nt = length(vol);
end
dataOut = spm_read_vols(vol);
end % local_read_NIfTI

function [vx] = local_check_shrinkage(IOI)
% Reshrinks image if necessary
if IOI.res.shrinkageOn
    vx = [IOI.res.shrink_x IOI.res.shrink_y 1];
else
    vx = [1 1 1];
end
end % local_check_shrinkage

function filtY = local_filter_time_course(Y, z, p, k, s1, brainMask)
slice2D = false;
% Filters an images time-course along the time dimension.
if ndims(Y) == 4,
    % We assume 4th dimension is time
    filtY = zeros(size(Y));
elseif ndims(Y) == 3,
    % We assume 3rd dimension is time
    filtY = zeros([size(Y,1) size(Y,2) 1 size(Y,3)]);
else
    % Only one temporal point, we have a single slice
    filtY = zeros([size(Y,1) size(Y,2) 1 1]);
    slice2D = true;
end
if ~slice2D
    % Can these loops be vectorized? //EGC
    for iX = 1:size(filtY,1)
        for iY = 1:size(filtY,2)
            if ~isempty(brainMask)
                if brainMask(iX,iY) == 1
                    % Only non-masked pixels are band-passs filtered
                    filtY(iX,iY,1,:) = temporalBPFrun(squeeze(Y(iX,iY,:)), z, p, k);
                end
            else
                filtY(iX,iY,1,:) = temporalBPFrun(squeeze(Y(iX,iY,:)), z, p, k);
            end
        end
    end
else
    fprintf('2D flow matrix found at session %2d.\n',s1)
    filtY = Y;
end
end % local_filter_time_course

% EOF
