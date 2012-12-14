function out = ioi_flow_run(job)
% Computes decorrelation velocity from speckle contrast in IR laser images.
%_______________________________________________________________________________
% Copyright (C) 2011 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________

% select a subset of sessions
[all_sessions selected_sessions] = ioi_get_sessions(job);
RemoveLC = job.RemoveLC;
for SubjIdx=1:length(job.IOImat)
    try
        tic
        %Load IOI.mat information
        [IOI IOImat dir_ioimat]= ioi_get_IOI(job,SubjIdx);
        
        if ~isfield(IOI.res,'flowOK') || job.force_redo
            str_flow = 'F';
            str_laser = IOI.color.laser;
            tmp_str_contrast = ['_C_'];
            tmp_str_flow = ['_' str_flow '_'];
            IOI.color.flow = str_flow;
            if ~(IOI.color.eng==str_flow)
                IOI.color.eng = [IOI.color.eng str_flow];
            end
                        
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
                                spm_progress_bar('Init', length(fname_list), sprintf('Flow computation, session %d\n',s1), 'Files');
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
                                    contrast=spm_read_vols(vol);
                                    image_flow=zeros(nx,ny,1,nt);
                                    for i3=1:nt
                                        image_flow(:,:,1,i3)=private_flow_from_contrast(squeeze(contrast(:,:,1,i3)),job.configuration.integ_time);
                                        % Normalization of the decorrelation velocity
                                        if isfield(job, 'flow_lambda_norm')
                                            if isfield (job.flow_lambda_norm, 'flow_lambda_norm_On')
                                                % Establishes the relationship on the correlation time and the velocity of the scattering particles.
                                                image_flow(:,:,1,i3) = image_flow(:,:,1,i3).*job.flow_lambda_norm.flow_lambda_norm_On.IR_laser_lambda/(2*pi);
                                            else
                                                % Decorrelation velocity is inversely proportional to correlation time.
                                            end
                                        end
                                    end
                                    
                                    if IOI.res.shrinkageOn
                                        vx=[IOI.res.shrink_x IOI.res.shrink_y 1];
                                    else
                                        vx = [1 1 1];
                                    end
                                  
                                    %save - substitute 'C' for 'F' in file name
                                    fname_new = regexprep(fname, tmp_str_contrast , tmp_str_flow);
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
                        disp(['Flow calculation for session ' int2str(s1) ' complete']);
                    end
                end
            end
            IOI.res.flowOK = 1;
            save(IOImat,'IOI');
            
            %remove LC images
            if RemoveLC
                for s1=1:length(IOI.sess_res)  
                    if all_sessions || sum(s1==selected_sessions)
                         if length(IOI.sess_res{s1}.fname)>=Li
                            fname_list = IOI.sess_res{s1}.fname{Li};                        
                            if ~isempty(fname_list)
                                for f1=1:length(fname_list)
                                    remove_vols_each_color(IOI,str_laser,f1,s1);
                                    if isfield(IOI.color,'str_contrast')
                                        remove_vols_each_color(IOI,IOI.color.str_contrast,f1,s1);
                                    end
                                    
                                end
                            end
                        end
                    end
                end
            end
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
%%%END

function speed = private_flow_from_contrast(contrast,T)
% trouve la vitesse à partir des images de contraste
%[nx ny] = size(contrast);

% Correct for points that cannot be...
contrast(isnan(contrast)|contrast<0)=0;

% Not sure about this, verify
contrast2=contrast(3:end-2,3:end-2);
mmean=mean(contrast2(1:end));
sstd=std(contrast2(1:end));

% Build non-linear curve between contrast and correlation time (tau)
tau=(logspace(-11,-2,30).^.5); % Correlation time
K  = ((tau/(2*T)).*(1-exp(-2*T*ones(size(tau))./tau))).^(1/2);

% Find values for which the mean contrast is in the middle
[dummy index1]=find(K>(mmean-3*sstd),1);
[dummy index2]=find(K>(mmean+3*sstd),1);
if isempty(index1), index1=1; end
if isempty(index2)||index2==index1, index2=30; end

% For these values, build a log-linear vector on which contrast is computed
Tau2=(logspace(log10(tau(index1)),log10(tau(index2)),40));
K  = ((Tau2/(2*T)).*(1-exp(-2*T*ones(size(Tau2))./Tau2))).^(1/2);

% Add limit points for interpolation
Tau2=[Tau2(1) Tau2 Tau2(end)];
K= [ 0 K 1e30];
% Interpolate contrast image on these values to obtain correlation time
Tau3=interp1(K,Tau2,contrast); %This is SLOW
% Get speed from correlation time
speed=1./Tau3;
