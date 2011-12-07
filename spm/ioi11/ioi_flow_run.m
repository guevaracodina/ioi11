function out = ioi_flow_run(job)
%select a subset of sessions
if isfield(job.session_choice,'select_sessions')
    all_sessions = 0;
    selected_sessions = job.session_choice.select_sessions.selected_sessions;
else
    all_sessions = 1;
end
try
    RemoveLC = job.RemoveLC;
catch
    RemoveLC = 1;
end
for SubjIdx=1:length(job.IOImat)
    try
        tic
        clear IOI
        %Load IOI.mat information
        IOImat = job.IOImat{SubjIdx};
        load(IOImat);
        if ~isfield(IOI.res,'flowOK') || job.force_redo
            
            [dir_ioimat dummy] = fileparts(job.IOImat{SubjIdx});
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
            IOI.res.flow.T=job.configuration.integ_time;
            
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
                                %loop over data files
                                for f1=1:length(fname_list)
                                    fname = fname_list{f1};
                                    %load data
                                    vol=spm_vol(fname);
                                    nx=vol(1).dim(1);
                                    ny=vol(1).dim(2);
                                    nt = length(vol); 
                                    win2 = ones(wsize,wsize);
                                    laser=spm_read_vols(vol);
                                    %mean_laser = zeros(size(laser));
                                    %std_laser = zeros(size(laser));
                                    image_flow=zeros(nx,ny,1,nt);
                                    OPTIONS.GPU = 0;
                                    OPTIONS.Power2Flag = 0;
                                    OPTIONS.Brep = 0;
                                    for i3=1:nt
                                        tmp_laser = squeeze(laser(:,:,:,i3));
                                        std_laser=stdfilt(tmp_laser,win2);
                                        %this is much faster (4 times) than conv2
                                        %tic
                                        %THIS IS SLOW:
                                        mean_laser = convnfft(tmp_laser,win2,'same',1:2,OPTIONS)/wsize^2;
                                        %toc
                                        contrast=std_laser./mean_laser;
                                        image_flow(:,:,1,i3)=private_flow_from_contrast(contrast,job.configuration.integ_time);
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
                                end
                                IOI.sess_res{s1}.fname{IOI.color.eng==str_flow} = fname_new_list;
                            end
                        end
                        %toc
                        disp(['Flow calculation for session ' int2str(s1) ' complete']);
                    end
                end
            end
            IOI.res.flowOK = 1;
            if isfield(job.IOImatCopyChoice,'IOImatCopy')
                newDir = job.IOImatCopyChoice.IOImatCopy.NewIOIdir;
                newDir = fullfile(dir_ioimat,newDir);
                if ~exist(newDir,'dir'),mkdir(newDir); end
                IOImat = fullfile(newDir,'IOI.mat');
            end
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
        toc
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
