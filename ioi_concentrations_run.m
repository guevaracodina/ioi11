function out = ioi_concentrations_run(job)
rescaling_factor = 1e6;
lambda1=450;
lambda2=700;
npoints=512;
str_HbO = 'O'; %oxy
str_HbR = 'D'; %deoxy
tmp_str_HbO = ['_' str_HbO '_'];
tmp_str_HbR = ['_' str_HbR '_'];
%basehbt1 = job.basehbt1;
baseline_hbt = job.configuration.HbT0;
baseline_hbo = job.configuration.HbO0;
baseline_hbr = job.configuration.HbR0;
RemoveRGY = job.RemoveRGY;
try
    if isfield(job.normalization_choice,'select_norm_session')
        normalization_choice = 1;
        selected_norm_session = job.normalization_choice.select_norm_session.selected_norm_session;
    else
        normalization_choice = 0;
    end
catch
    normalization_choice = 0;
end
for SubjIdx=1:length(job.IOImat)
    try
        tic
        clear IOI
        [all_sessions selected_sessions] = ioi_get_sessions(job);
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
        
        if ~isfield(IOI.res,'concOK') || job.force_redo
            IOI.conc.baseline_hbt = baseline_hbt;
            IOI.conc.baseline_hbo = baseline_hbo;
            IOI.conc.baseline_hbr = baseline_hbr;
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
            if isfield(IOI,'sess_raw')
                whichSystem = 0; %old
            else
                whichSystem = 1; %new
            end
            eps_pathlength = ioi_epsilon_pathlength(lambda1,lambda2,npoints,whichSystem,whichCurve,baseline_hbt,baseline_hbo,baseline_hbr);
            
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
                        if normalization_choice
                            %get median images
                            for i1=1:length(hasRGY)
                                try
                                    V0 = spm_vol(IOI.sess_res{selected_norm_session}.fname_median{hasRGY(i1)});
                                    V1 = spm_vol(IOI.sess_res{s1}.fname_median{hasRGY(i1)});
                                catch
                                    V0 = spm_vol([IOI.sess_res{selected_norm_session}.fname_median{hasRGY(i1)} '.nii']);
                                    V1 = spm_vol([IOI.sess_res{s1}.fname_median{hasRGY(i1)} '.nii']);
                                end
                                med0(i1,:,:) = spm_read_vols(V0);
                                med1(i1,:,:) = spm_read_vols(V1);
                            end
                        end
                        if ~isempty(fname_list)
                            % Initialize progress bar
                            spm_progress_bar('Init', length(fname_list), sprintf('Concentrations computation, session %d\n',s1), 'Files');
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
                                    %normalization_choice
                                    if normalization_choice
                                        slice = ioi_normalization_choice(slice,med0,med1);
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
                                        if normalization_choice
                                            slice = ioi_normalization_choice(slice,med0,med1);
                                        end
                                        slice=reshape(slice,[length(vols),nx*ny]);
                                        slice=Ainv*slice;
                                        % A this point we have HbO and HbR, reform image
                                        slice=reshape(slice,[2,nx,ny]);
                                        image_hbo(:,:,1,i1)=slice(1,:,:);
                                        image_hbr(:,:,1,i1)=slice(2,:,:);
                                    end
                                end
                                oNaN = sum(isnan(image_hbo(:)));
                                rNaN = sum(isnan(image_hbr(:)));
                                oInf = sum(isinf(image_hbo(:)));
                                rInf = sum(isinf(image_hbr(:)));
                                omax = max(image_hbo(:));
                                rmax = max(image_hbr(:));
                                if oNaN
                                    IOI = disp_msg(IOI,[int2str(oNaN) ' NaN in HbO in session ' int2str(s1)]);
                                    image_hbo(isnan(image_hbo(:))) = omax;
                                end
                                if rNaN
                                    IOI = disp_msg(IOI,[int2str(rNaN) ' NaN in HbR in session ' int2str(s1)]);
                                    image_hbr(isnan(image_hbr(:))) = rmax;
                                end
                                if oInf
                                    IOI = disp_msg(IOI,[int2str(oInf) ' Inf in HbO in session ' int2str(s1)]);
                                    image_hbo(isinf(image_hbo(:))) = omax;
                                end
                                if rInf
                                    IOI = disp_msg(IOI,[int2str(rInf) ' Inf in HbR in session ' int2str(s1)]);
                                    image_hbr(isinf(image_hbr(:))) = rmax;
                                end
                                %save - substitute 'O' and 'D' in file name
                                fname = IOI.sess_res{s1}.fname{hasRGY(1)}{f1};
                                if isfield(job.IOImatCopyChoice,'IOImatCopy')
                                    [dir0 fil0 ext0] = fileparts(fname);
                                    fdir = fullfile(newDir,['S' gen_num_str(s1,2)]);
                                    if ~exist(fdir,'dir'), mkdir(fdir); end
                                    fname = fullfile(fdir,[fil0 ext0]);
                                end
                                
                                fname_new_HbO = regexprep(fname,tmp_str_avail_color ,tmp_str_HbO);
                                fname_new_HbO_list = [fname_new_HbO_list; fname_new_HbO];
                                ioi_save_nifti(image_hbo, fname_new_HbO, vx);
                                fname_new_HbR = regexprep(fname,tmp_str_avail_color ,tmp_str_HbR);
                                fname_new_HbR_list = [fname_new_HbR_list; fname_new_HbR];
                                ioi_save_nifti(image_hbr, fname_new_HbR, vx);
                                % Update progress bar
                                spm_progress_bar('Set', f1);
                            end % files loop
                        end
                        % Clear progress bar
                        spm_progress_bar('Clear');
                        IOI.sess_res{s1}.fname{IOI.color.eng==str_HbO} = fname_new_HbO_list;
                        IOI.sess_res{s1}.fname{IOI.color.eng==str_HbR} = fname_new_HbR_list;
                        %toc
                        disp(['Concentration calculation for session ' int2str(s1) ' complete']);
                    end
                end
            end
            
            IOI.res.concOK = 1;
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
        disp(['Elapsed time: ' datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS')]);
        disp(['Subject ' int2str(SubjIdx) ' complete']);
    catch exception
        disp(exception.identifier)
        disp(exception.stack(1))
        out.IOImat{SubjIdx} = job.IOImat{SubjIdx};
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

function slice = ioi_normalization_choice(slice,med0,med1)
if length(size(slice)) == 4
    med0 = repmat(med0, [1 1 1 size(slice,4)]);
    med1 = repmat(med1, [1 1 1 size(slice,4)]);
end
slice = slice - log(med1./med0);
