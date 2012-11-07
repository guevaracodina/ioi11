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
fprintf('Work in progress...\nEGC\n')
out.IOImat = job.IOImat;
return
% ------------------------------------------------------------------------------

% Select a subset of sessions
[all_sessions selected_sessions] = ioi_get_sessions(job);
% String that identifies CMRO2
str_CMRO2 = 'M';
tmp_str_CMRO2 = ['_' str_CMRO2 '_'];
str_flow = 'F';
tmp_str_flow = ['_' str_flow '_'];
str_HbO = 'O';
str_HbR = 'D';
% Vascular weighting constants
gammaT = job.constants.gammaT;
gammaR = job.constants.gammaR;

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
                % Update color structure in IOI matrix
                IOI.color.CMRO2 = str_CMRO2;
                % Update string of colors
                if ~(IOI.color.eng==str_CMRO2)
                    IOI.color.eng = [IOI.color.eng str_CMRO2];
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
                                spm_progress_bar('Init', length(fname_list_flow), sprintf('CMRO2 computation, session %d\n',s1), 'Files');
                                % Loop over data files
                                for f1 = 1:length(fname_list_flow)
                                    if job.MemoryManagementMenu
                                        % Load all images at once
                                        [dataFlow nx(1) ny(1) nt(1)] = local_read_NIfTI(fname_list_flow{f1});
                                        [dataHbO nx(2) ny(2) nt(2)] = local_read_NIfTI(fname_list_HbO{f1});
                                        [dataHbR nx(3) ny(3) nt(3)] = local_read_NIfTI(fname_list_HbR{f1});
                                        % If HbO, HbR & flow data have the same
                                        % size
                                        if all(nx == nx(1)) && all(ny == ny(1)) && all(nt == nt(1))
                                            %% CMRO2 computation
                                            image_CMRO2 = ioi_cmro2_compute(dataFlow, dataHbO, dataHbR, gammaT, gammaR);
                                        else
                                            error('ioi_cmro2_run:different_Sizes_Flow_HbO_HbR', 'Flow and Hb data have different sizes.');
                                        end
                                    else
                                        % Load images one by one
                                        % TO DO...
                                    end
                                    % Reshrinks image if necessary
                                    vx = local_check_shrinkage(IOI);
                                    % Save - substitute 'F' for 'M' in file name
                                    fname_new = regexprep(fname_list_flow{f1}, tmp_str_flow , tmp_str_CMRO2);
                                    % Append new file name to file name list
                                    fname_new_list = [fname_new_list; fname_new];
                                    % Save as NIfTI file
                                    ioi_save_nifti(image_CMRO2, fname_new, vx);
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
        disp(['Subject ' int2str(SubjIdx) ' complete']);
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
% if vx(1) > 1 || vx(2) > 1
%     nx = size(vx(1):vx(1):size(image_CMRO2_in,1),2);
%     ny = size(vx(2):vx(2):size(image_CMRO2_in,2),2);
%     image_CMRO2_out = ioi_imresize(image_CMRO2_in, 1, nx, ny, vx(1), vx(2));
% end
end % local_check_shrinkage

% EOF
