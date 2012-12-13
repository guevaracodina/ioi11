function out = ioi_fcIOS_maps_run(job)
% 
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

% Loop over subjects
for SubjIdx = 1:length(job.IOImat)
    try
        tic
        % Load IOI.mat information
        [IOI IOImat dir_ioimat]= ioi_get_IOI(job,SubjIdx);
        % Load brain mask
        brainMaskVol = spm_vol(IOI.fcIOS.mask.fname);
        brainMask = logical(spm_read_vols(brainMaskVol));
        % Read anatomical NIFTI file
        im_anat = ioi_read_vol(IOI.res.file_anat);
        % Loop over sessions
        for s1 = 1:length(IOI.sess_res)
            if all_sessions || sum(s1==selected_sessions)
                % Loop over available colors
                for c1=1:length(IOI.sess_res{s1}.fname)
                    doColor = ioi_doColor(IOI,c1,IC);
                    if doColor
                        % Load correlation data
                        if size(brainMask,1)~= size(y,1)|| size(brainMask,2)~= size(y,2)
                            brainMask = ioi_MYimresize(brainMask, [size(y,1) size(y,2)]);
                        end
                        % Mask data
                        
                        % display overlay data
                       
                    end
                end % colors loop
            end
        end % sessions loop
        out.IOImat{SubjIdx} = IOImat;
        disp(['Elapsed time: ' datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS')]);
        disp(['Subject ' int2str(SubjIdx) ' (' IOI.subj_name ')' ' complete']);
    catch exception
        disp(exception.identifier)
        disp(exception.stack(1))
        out.IOImat{SubjIdx} = job.IOImat{SubjIdx};
    end
end % loop over subjects
end % Main function

% EOF
