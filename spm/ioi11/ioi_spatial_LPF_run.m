function out = ioi_spatial_LPF_run(job)
% Low-pass filtering of 2-D images with a rotationally symmetric gaussian kernel
%_______________________________________________________________________________
% Copyright (C) 2012 LIOM Laboratoire d'Imagerie Optique et Mol�culaire
%                    �cole Polytechnique de Montr�al
%_______________________________________________________________________________

% Get sessions info
[all_sessions selected_sessions] = ioi_get_sessions(job);

%Big loop over subjects
for SubjIdx=1:length(job.IOImat)
    try
        tic
        clear IOI
        %Load IOI.mat information
        IOImat = job.IOImat{SubjIdx};
        [dir_ioimat dummy] = fileparts(job.IOImat{SubjIdx});
        if isfield(job.IOImatCopyChoice,'IOImatCopy')
            newDir = job.IOImatCopyChoice.IOImatCopy.NewIOIdir;
            newDir = fullfile(dir_ioimat,newDir);
            if ~exist(newDir,'dir'),mkdir(newDir); end
            IOImat = fullfile(newDir,'IOI.mat');
        else
            newDir = dir_ioimat;
        end
        try
            load(IOImat);
        catch
            load(job.IOImat{SubjIdx});
        end
        if ~isfield(IOI.res,'concOK') % Concentrations OK
            disp(['No concentrations available for subject ' int2str(SubjIdx) ' ... skipping low-pass filtering']);
        else
            if ~isfield(IOI.res,'flowOK') % Flow OK
                disp(['No flow available for subject ' int2str(SubjIdx) ' ... skipping low-pass filtering']);
            else
                if ~isfield(IOI.fcIOS.LPF,'LPFOK') || job.force_redo
                    % Get colors to include information
                    IC = job.IC;
                    colorNames = fieldnames(IOI.color);
                    % Radius of the gaussian kernel
                    K.radius = round(job.spatial_LPF.spatial_LPF_On.spatial_LPF_radius);
                    if IOI.res.shrinkageOn %ne pas oublier de decommenter
                        vx = [IOI.res.shrink_x IOI.res.shrink_y 1];
                    else
                        vx = [1 1 1];
                    end
                    % Loop over sessions
                    for s1=1:length(IOI.sess_res)
                        if all_sessions || sum(s1==selected_sessions)
                            % Loop over available colors
                            for c1=1:length(IOI.sess_res{s1}.fname)
                                doColor = ioi_doColor(IOI,c1,IC);
                                if doColor
                                    colorOK = true;
                                    %skip laser - only extract for flow
                                    if ~(IOI.color.eng(c1)==IOI.color.laser)
                                        %% Correlation map
                                        % Results filenames
                                        fname_list = IOI.sess_res{s1}.fname{c1};
                                        % Color names
                                        colorNames = fieldnames(IOI.color);
                                        % Initialize progress bar
                                        spm_progress_bar('Init', length(fname_list), sprintf('Spatial LPF session %d, color %d (%s)\n',s1,c1,colorNames{1+c1}), 'Files');
                                        % Loop over files
                                        for f1 = 1:length(fname_list)
                                            try
                                                fname = fname_list{f1};
                                                vols = spm_vol(fname);
                                                imagesTimeCourse = spm_read_vols(vols);
                                                [K.k1 K.k2 d3 d4] = size(imagesTimeCourse);
                                                imagesTimeCourseLPF = zeros([K.k1 K.k2 d3 d4]);
                                                if K.k1 <= 1 || K.k2 <= 1
                                                    colorOK = false;
                                                end
                                            catch
                                                colorOK = false;
                                            end
                                            %time dimension in 3rd dimension for colors
                                            %R, G, Y, but in 4th dimension for O, D, F
                                            if c1==1 || c1==2 || c1==3
                                                nt = d3;
                                            else
                                                nt = d4;
                                            end
                                            % Low-passfiltering and saving
                                            K = ioi_spatial_LPF('set', K);
                                            
                                            for t1 = 1:nt,
                                                imagesTimeCourseLPF(:,:,1,t1) = ioi_spatial_LPF('lpf', K, squeeze(imagesTimeCourse(:,:,1,t1)));
                                            end
                                            % Overwrite the image
                                            ioi_save_nifti(imagesTimeCourseLPF, fname, vx);
                                            % Update progress bar
                                            spm_progress_bar('Set', f1);
                                        end
                                        % Clear progress bar
                                        spm_progress_bar('Clear');
                                        if colorOK
                                            fprintf('Spatial low-pass filtering for session %d and color %d (%s) completed\n',s1,c1,colorNames{1+c1})
                                        end
                                    end
                                end
                            end % colors loop
                        end
                    end % Sessions Loop
                    % LPF succesful!
                    IOI.fcIOS.LPF(1).LPFOK = true;
                    % Save LPF settings
                    IOI.fcIOS.LPF(1).radius = K.radius;
                    IOI.fcIOS.LPF(1).sigma = K.radius/2;
                    save(IOImat,'IOI');
                end % LPF OK or redo job
            end % Flow OK
        end % Concentrations OK
        disp(['Elapsed time: ' datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS')]);
        disp(['Subject ' int2str(SubjIdx) ' (' IOI.subj_name ')' ' complete']);
        out.IOImat{SubjIdx} = IOImat;
    catch exception
        out.IOImat{SubjIdx} = IOImat;
        disp(exception.identifier)
        disp(exception.stack(1))
    end
end

% EOF
