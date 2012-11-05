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

% Loop over subjects
for SubjIdx = 1:length(job.IOImat)
    try
        tic
        %Load IOI.mat information
        [IOI IOImat dir_ioimat]= ioi_get_IOI(job,SubjIdx);
        
        if ~isfield(IOI.res,'cmro2OK') || job.force_redo
            wsize=job.configuration.window_size;
            if (mod(wsize,2)==0)
                disp(['Flow Computation: Need to specify an odd window size, changing to: ' num2str(wsize+1)]);
                wsize=wsize+1;
                job.configuration.window_size=wsize;  % For reference
            end
            
            str_flow = 'F';
            str_laser = IOI.color.laser;
            tmp_str_laser = ['_' str_laser '_'];
            tmp_str_flow = ['_' str_flow '_'];
            IOI.color.flow = str_flow;
            if ~(IOI.color.eng==str_flow)
                IOI.color.eng = [IOI.color.eng str_flow];
            end
            
            IOI.res.flow.window_size=wsize;
            
            % Get integration time directly from the info.txt file //EGC
            try
                % Get the raw data folder
                [dummy,dirs] = cfg_getfile('FPListRec',IOI.dir.dir_subj_raw,'dir');
                % Get the info.txt from the 1st session as a string
                [file_info,dummy] = cfg_getfile('FPList',dirs{1},'info.txt');
                t0 = fileread(file_info{1});
                % Find the exposure time
                [startIndex endIndex]= regexp(t0, 'Exposure time', 'once');
                tmpString = t0(endIndex:end);
                % Find the time units prefix
                tFactorStr = regexp(tmpString, '\(..\)', 'match', 'once');
                switch tFactorStr
                    case '(us)'
                        tFactor = 1e6;
                    case '(ms)'
                        tFactor = 1e3;
                    otherwise
                        tFactor = 1;
                end
                % Gets the first number in the string as integration time
                intTimeStr = regexp(tmpString,'\d+','match', 'once');
                % Save integration time in job (in seconds)
                if ~isempty(intTimeStr)
                    job.configuration.integ_time = str2double(intTimeStr)/tFactor;                    
                end
                % Update value in IOI matrix
                IOI.res.flow.T = job.configuration.integ_time;
            catch
                % Do as usual, use whatever value the user entered
                IOI.res.flow.T = job.configuration.integ_time;
            end
            
            %Loop over sessions
            if isfield(IOI,'sess_res')
                for s1=1:length(IOI.sess_res)
                    if all_sessions || sum(s1==selected_sessions)
                        %tic
                        %check that laser speckle is available
                        Li = find(IOI.color.eng==str_laser); %index for laser speckle
                        if length(IOI.sess_res{s1}.fname)>=Li
                            fname_list = IOI.sess_res{s1}.fname{Li};
                            fname_new_list = {};
                            if ~isempty(fname_list)
                                % Initialize progress bar
                                spm_progress_bar('Init', length(fname_list), sprintf('CMRO2 computation, session %d\n',s1), 'Files');
                                %loop over data files
                                for f1=1:length(fname_list)
                                    fname = fname_list{f1};
                                    %load data
                                    vol = spm_vol(fname);
                                    nx = vol(1).dim(1);
                                    ny = vol(1).dim(2);
                                    if length(vol) == 1
                                        nt = vol(1).dim(3);
                                    else
                                        nt = length(vol); 
                                    end
                                    % NOTE: nt is not necessarily the largest
                                    % dimension of vol //EGC
                                    % nt = length(vol); 
                                    win2 = ones(wsize,wsize);
                                    laser=spm_read_vols(vol);
                                    %mean_laser = zeros(size(laser));
                                    %std_laser = zeros(size(laser));
                                    image_flow=zeros(nx,ny,1,nt);
                                    OPTIONS.GPU = 0;
                                    OPTIONS.Power2Flag = 0;
                                    OPTIONS.Brep = 0;
                                    for i3=1:nt
                                        tmp_laser = squeeze(laser(:,:,i3));
                                        std_laser=stdfilt(tmp_laser,win2);
                                        %this is much faster (4 times) than conv2
                                        %tic
                                        %THIS IS SLOW:
                                        mean_laser = convnfft(tmp_laser,win2,'same',1:2,OPTIONS)/wsize^2;
                                        %toc
                                        contrast=std_laser./mean_laser;
                                        image_flow(:,:,1,i3)=private_flow_from_contrast(contrast,job.configuration.integ_time);
                                        % Here I compute the CMRO2
                                        cmro2 = ioi_cmro2_compute();
                                    end
                                    
                                    if IOI.res.shrinkageOn
                                        vx=[IOI.res.shrink_x IOI.res.shrink_y 1];
                                    else
                                        vx = [1 1 1];
                                    end
                                    if vx(1) > 1 || vx(2) > 1
                                        nx=size(vx(1):vx(1):size(image_flow,1),2);
                                        ny=size(vx(2):vx(2):size(image_flow,2),2);
                                        image_flow = ioi_imresize(image_flow,0,nx,ny,vx(1),vx(2));
                                    end
                                    
                                    %save - substitute 'L' for 'F' in file name
                                    fname_new = regexprep(fname, tmp_str_laser , tmp_str_flow);
                                    fname_new_list = [fname_new_list; fname_new];
                                    ioi_save_nifti(image_flow, fname_new, vx);
                                    % Update progress bar
                                    spm_progress_bar('Set', f1);
                                end % files loop
                                % Clear progress bar
                                spm_progress_bar('Clear');
                                IOI.sess_res{s1}.fname{IOI.color.eng==str_flow} = fname_new_list;
                            end
                        end
                        %toc
                        disp(['CMRO2 calculation for session ' int2str(s1) ' complete']);
                    end
                end
            end
            IOI.res.cmro2OK = 1;
            save(IOImat,'IOI');
            
            
        end
        out.IOImat{SubjIdx} = IOImat;
        disp(['Elapsed time: ' datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS')]);
        disp(['Subject ' int2str(SubjIdx) ' complete']);
    catch exception
        disp(exception.identifier)
        disp(exception.stack(1))
        out.IOImat{SubjIdx} = job.IOImat{SubjIdx};
    end
end

% EOF
