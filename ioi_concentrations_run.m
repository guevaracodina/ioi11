function out = ioi_concentrations_run(job)
rescaling_factor = 1e6;
lambda1=450;
lambda2=700;
npoints=512;
str_HbO = 'O'; %oxy
str_HbR = 'D'; %deoxy
tmp_str_HbO = ['_' str_HbO '_'];
tmp_str_HbR = ['_' str_HbR '_'];
try
    RemoveRGY = job.RemoveRGY;
catch
    RemoveRGY = 1;
end
for SubjIdx=1:length(job.IOImat)
    try
        tic
        clear IOI
        %Load IOI.mat information
        IOImat = job.IOImat{SubjIdx};
        load(IOImat);
        %select a subset of sessions
        if isfield(job.session_choice,'select_sessions')
            all_sessions = 0;
            selected_sessions = job.session_choice.select_sessions.selected_sessions;
        else
            all_sessions = 1;
        end
        if ~isfield(IOI.res,'concOK') || job.force_redo
            
            [dir_ioimat dummy] = fileparts(job.IOImat{SubjIdx});
            IOI.color.HbO = str_HbO;
            IOI.color.HbR = str_HbR;
            if ~(IOI.color.eng==str_HbO)
                IOI.color.eng = [IOI.color.eng str_HbO];
            end
            if ~(IOI.color.eng==str_HbR)
                IOI.color.eng = [IOI.color.eng str_HbR];
            end
            whichCurve=job.configuration.pathlength;
            if IOI.res.shrinkageOn %ne pas oublier de decommenter
                vx=[IOI.res.shrink_x IOI.res.shrink_y 1];
            else
                vx = [1 1 1];
            end
            %This function is hard-coded with wavelengths in the order Red,
            %Green, Yellow
            eps_pathlength = ioi_epsilon_pathlength(lambda1,lambda2,npoints,whichCurve);
            %Loop over sessions
            if isfield(IOI,'sess_res')
                for s1=1:length(IOI.sess_res)
                    if all_sessions || sum(s1==selected_sessions)
                        %tic
                        %check which colors are available
                        tmp_hasRGY = IOI.sess_res{s1}.hasRGY;
                        hasRGY = [];
                        for i1=1:length(tmp_hasRGY)
                            hasRGY = [hasRGY find(IOI.color.eng==tmp_hasRGY(i1))];
                        end
                        A=eps_pathlength(hasRGY,:);
                        if (size(A,1)<2)
                            disp('ioi_run_concentrations: Need at least two wavelenghts to compute concentrations')
                        else
                            Ainv=rescaling_factor*pinv(A);
                        end
                        fname_new_HbO_list = {};
                        fname_new_HbR_list = {};                                              
                        fname_list = IOI.sess_res{s1}.fname{hasRGY(1)};
                        str_avail_color = tmp_hasRGY(1);
                        tmp_str_avail_color = ['_' str_avail_color '_'];
                        if ~isempty(fname_list)
                            for f1=1:length(fname_list)
                                vols = {}; vi = 1;
                                %get volume headers for each color 
                                [vols vi] = get_vols_each_color(IOI,vols,vi,IOI.color.red,f1,s1);
                                [vols vi] = get_vols_each_color(IOI,vols,vi,IOI.color.green,f1,s1);
                                [vols vi] = get_vols_each_color(IOI,vols,vi,IOI.color.yellow,f1,s1);
                                
                                % Loop over volumes within each files
                                nx = vols{1}(1).dim(1);
                                ny = vols{1}(1).dim(2);
                                nt = length(vols{1});
                                image_hbo=zeros(nx,ny,1,nt);
                                image_hbr=zeros(nx,ny,1,nt);
                                clear slice
                                if job.MemoryManagementMenu %load all at once
                                    for c1 = 1:length(vols) %for each color
                                        slice(c1,:,:,:) = squeeze(spm_read_vols(vols{c1}));
                                    end
                                    slice=reshape(slice,[length(vols),nx*ny*nt]);
                                    slice=Ainv*slice;
                                    % A this point we have HbO and HbR, reform image
                                    slice=reshape(slice,[2,nx,ny,nt]);
                                    image_hbo(:,:,1,:)=slice(1,:,:,:);
                                    image_hbr(:,:,1,:)=slice(2,:,:,:);
                                else %load one image at a time
                                    for i1=1:nt
                                        for c1 = 1:length(vols)
                                            slice(c1,:,:)=ioi_read_time_vol(vols{c1},i1);
                                        end
                                        slice=reshape(slice,[length(vols),nx*ny]);
                                        slice=Ainv*slice;
                                        % A this point we have HbO and HbR, reform image
                                        slice=reshape(slice,[2,nx,ny]);
                                        image_hbo(:,:,1,i1)=slice(1,:,:);
                                        image_hbr(:,:,1,i1)=slice(2,:,:);
                                    end
                                end
                                %save - substitute 'O' and 'D' in file name
                                fname = IOI.sess_res{s1}.fname{hasRGY(1)}{f1};
                                fname_new_HbO = regexprep(fname,tmp_str_avail_color ,tmp_str_HbO);
                                fname_new_HbO_list = [fname_new_HbO_list; fname_new_HbO];
                                ioi_save_nifti(image_hbo, fname_new_HbO, vx);
                                fname_new_HbR = regexprep(fname,tmp_str_avail_color ,tmp_str_HbR);
                                fname_new_HbR_list = [fname_new_HbR_list; fname_new_HbR];
                                ioi_save_nifti(image_hbr, fname_new_HbR, vx);
                            end
                        end
                        IOI.sess_res{s1}.fname{IOI.color.eng==str_HbO} = fname_new_HbO_list;
                        IOI.sess_res{s1}.fname{IOI.color.eng==str_HbR} = fname_new_HbR_list;
                        %toc
                        disp(['Concentration calculation for session ' int2str(s1) ' complete']);
                    end
                end
            end
            
            IOI.res.concOK = 1;
            if isfield(job.IOImatCopyChoice,'IOImatCopy')
                newDir = job.IOImatCopyChoice.IOImatCopy.NewIOIdir;
                newDir = fullfile(dir_ioimat,newDir);
                if ~exist(newDir,'dir'),mkdir(newDir); end
                IOImat = fullfile(newDir,'IOI.mat');
            end
            save(IOImat,'IOI');
            
            %remove RGY images
            if RemoveRGY
                for s1=1:length(IOI.sess_res)
                    if all_sessions || sum(s1==selected_sessions)
                        fname_list = IOI.sess_res{s1}.fname{hasRGY(1)};
                        for c1 = 1:length(hasRGY)
                            if ~isempty(fname_list)
                                for f1=1:length(fname_list)
                                    remove_vols_each_color(IOI,tmp_hasRGY(c1),f1,s1);
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
    end
end

function [vols vi] = get_vols_each_color(IOI,vols,vi,str_color,f1,s1)
c1 = find(IOI.color.eng==str_color);
if length(IOI.sess_res{s1}.fname)>=c1
    fname_list = IOI.sess_res{s1}.fname{c1};
    if ~isempty(fname_list)
        fname = fname_list{f1};
        vols{vi}=spm_vol(fname);
        vi = vi+1;
    end
end
