function out = ioi_GLM_on_ROI_run(job)
%Volterra
volt = job.volt;
%bases
bases = job.bases;
%HPF filter
if isfield(job.hpf_butter,'hpf_butter_On')
    hpf_butter_On = 1;
    hpf_butter_freq = job.hpf_butter.hpf_butter_On.hpf_butter_freq;
    hpf_butter_order = job.hpf_butter.hpf_butter_On.hpf_butter_order;
else
    hpf_butter_On = 0;
end
%LPF filter
fwhm = job.lpf_gauss.fwhm1;

%select a subset of sessions
if isfield(job.session_choice,'select_sessions')
    all_sessions = 0;
    selected_sessions = job.session_choice.select_sessions.selected_sessions;
else
    all_sessions = 1;
end
%select a subset of ROIs
if isfield(job.ROI_choice,'select_ROIs')
    all_ROIs = 0;
    selected_ROIs = job.ROI_choice.select_ROIs.selected_ROIs;
else
    all_ROIs = 1;
end
%select onsets
if isfield(job.stim_choice,'electro_stims')
    electro_stims = 1;
    %number of standard deviations
    nSD = job.stim_choice.electro_stims.nSD;
    %minimal distance between peaks in seconds
    dP = job.stim_choice.electro_stims.dP/1000;
else
    %default stims
    electro_stims = 0;
end
%Big loop over subjects
for SubjIdx=1:length(job.IOImat)
    try
        tic
        clear IOI ROI SPM
        %Load IOI.mat information
        IOImat = job.IOImat{SubjIdx};
        load(IOImat);
        if ~isfield(IOI.res,'ROIOK')
            disp(['No ROI available for subject ' int2str(SubjIdx) ' ... skipping series extraction']);
        else
            if ~isfield(IOI.res,'seriesOK')
                disp(['Extracted series not available for subject ' int2str(SubjIdx) ' ... skipping GLM']);
            else
                if ~isfield(IOI.res,'GLMOK') || job.force_redo
                    if isfield(IOI,'X'), IOI = rmfield(IOI,'X'); end
                    [dir_ioimat dummy] = fileparts(job.IOImat{SubjIdx});
                    %loop over sessions
                    for s1=1:length(IOI.sess_res)
                        if all_sessions || sum(s1==selected_sessions)
                            %Electrophysiology, for each subject and session
                            if electro_stims
                                [pkh ons] = ioi_get_onsets(IOI,s1,nSD,dP); %pk in seconds; pkh in arbitrary units
                                dur = 1;
                                name = '';
                            else
                                ons = IOI.sess_res{s1}.onsets{1}; %already in seconds *IOI.dev.TR;
                                dur = IOI.sess_res{s1}.durations{1}; %*IOI.dev.TR;
                                name = IOI.sess_res{s1}.names{1};
                            end
                            %convolve with hemodynamic response function
                            [X U] = ioi_get_X(IOI,name,ons,dur,s1,bases,volt);
                            IOI.Sess(s1).U = U; %store onsets for each session
                            IOI.X{s1}.X0 = X;
                            %filter X - HPF
                            if hpf_butter_On
                                X = ButterHPF(1/IOI.dev.TR,hpf_butter_freq,hpf_butter_order,X);
                             end
                            %add a constant
                            X = [X ones(size(X,1),1)];
                            %get K for low pass filtering:
                            K = get_K(X,fwhm,IOI.dev.TR);
                            %filter X - LPF
                            %calculate X inverse
                            %Xm = pinv(X);
                            KX = spm_sp('Set', spm_filter_HPF_LPF_WMDL(K,X));
                            KX.X = full(KX.X); %Filtered X
                            Xm = spm_sp('x-',KX); % projector
                            %covariance
                            bcov = Xm * K.KL;
                            bcov = bcov * bcov';
                            IOI.X{s1}.X = X;
                            IOI.X{s1}.Xm = Xm;
                            IOI.X{s1}.bcov = bcov;
                            %load ROI
                            load(IOI.ROI.ROIfname);
                            %approximate calculation of effective degrees of freedom
                            [trRV trRVRV] = approx_trRV(KX.X,Xm,K.KL);
                            IOI.X{s1}.trRV = trRV;
                            IOI.X{s1}.trRVRV = trRVRV;
                            IOI.X{s1}.erdf = (trRV)^2/trRVRV;
                            %loop over available colors
                            for c1=1:length(IOI.sess_res{s1}.fname)
                                if ~(IOI.color.eng(c1)==IOI.color.laser)
                                    %Loop over ROIs
                                    for r1=1:length(ROI)
                                        if all_ROIs || sum(r1==selected_ROIs)
                                            y = ROI{r1}{s1,c1};
                                            if ~isempty(y)
                                                %filtering of the data: HPF
                                                if hpf_butter_On
                                                    y =  ButterHPF(1/IOI.dev.TR,hpf_butter_freq,hpf_butter_order,y);
                                                end
                                                %filtering of the data: LPF (Gaussian)
                                                y = spm_filter_HPF_LPF_WMDL(K,y')';
                                                %GLM inversion: calculating beta and residuals - would be SPM.Vbeta,
                                                b = Xm * y'; % beta : least square estimate
                                                IOI.X{1,s1}.b{r1,c1} = b;
                                                %Compute t stat
                                                res = y'-X*b;
                                                res2 = sum(res.^2);
                                                IOI.X{1,s1}.r(r1,c1) = res2; % Residuals
                                                IOI.X{1,s1}.t(r1,c1) = b(1)/(res2*bcov(1,1)/trRV)^0.5;
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                    
                    IOI.res.GLMOK = 1;
                    if isfield(job.IOImatCopyChoice,'IOImatCopy')
                        newDir = job.IOImatCopyChoice.IOImatCopy.NewIOIdir;
                        newDir = fullfile(dir_ioimat,newDir);
                        if ~exist(newDir,'dir'),mkdir(newDir); end
                        IOImat = fullfile(newDir,'IOI.mat');
                    end
                    save(IOImat,'IOI');
                end
            end
        end
        toc
        disp(['Subject ' int2str(SubjIdx) ' complete']);
        out.IOImat{SubjIdx} = IOImat;
        
    catch exception
        disp(exception.identifier)
        disp(exception.stack(1))
    end
end
end
function [pkh pk] = ioi_get_onsets(IOI,s1,nSD,nP)
%load raw electrophysiology vector
load(IOI.res.el{s1});
%remove time stamps for actual or spurious stimulations
ind = el>2; ind2 = [ind(7:end) false false false false false false];
el(ind)=el(ind2);
SD = std(el); %used 0.05 before
MN = mean(el);
sf = 10000;
%find onsets peaks: pkh: peak height; pk: onset time at sf sampling frequency
[pkh pk] = findpeaks(el,'MINPEAKHEIGHT',MN+nSD*SD,'MINPEAKDISTANCE',floor(sf*nP)); %25ms
pk = pk/sf;
%we2 = zeros(1,ns);
%we2(pk) = -0.05;
%inx = 0*1e6+(1:1e6);
%figure; plot(ve(inx),'k'); hold on; plot(we2(inx),'r'); hold off
end

function [X U] = ioi_get_X(IOI,name,ons,dur,s1,bases,volt)
SPM.xBF.dt = IOI.dev.TR;
SPM.xBF.T = 1;
SPM.xBF.T0 = 0;
SPM.xBF.UNITS = 'secs';
nambase = fieldnames(bases);
if ischar(nambase)
    nam=nambase;
else
    nam=nambase{1};
end
if strcmp(fieldnames(bases),'hrf') || strcmp(fieldnames(bases),'rat') ...
        || strcmp(fieldnames(bases),'mouse')
    %canonical HRF or rat or mouse
    if all(bases.(nam).derivs == [0 0])
        SPM.xBF.name = 'hrf';
    elseif all(bases.(nam).derivs == [1 0])
        SPM.xBF.name = 'hrf (with time derivative)';
    elseif all(bases.(nam).derivs == [1 1])
        SPM.xBF.name = 'hrf (with time and dispersion derivatives)';
    else
        disp('Unrecognized hrf derivative choices.')
    end
    SPM.xBF.nam = nam;
else
    switch nam,
        case 'fourier',
            SPM.xBF.name = 'Fourier set';
        case 'fourier_han',
            SPM.xBF.name = 'Fourier set (Hanning)';
        case 'gamma',
            SPM.xBF.name = 'Gamma functions';
        case 'fir',
            SPM.xBF.name = 'Finite Impulse Response';
        otherwise
            error('Unrecognized hrf derivative choices.')
    end
    SPM.xBF.length = bases.(nam).length;
    SPM.xBF.order  = bases.(nam).order;
end

SPM.xBF = spm_get_bf_rat_mouse(SPM.xBF);
if size(SPM.xBF.bf,1) == 1 || size(SPM.xBF.bf,2) == 1
    SPM.xBF.bf = SPM.xBF.bf/sum(SPM.xBF.bf); %normalize
end

% Get inputs, neuronal causes or stimulus functions U
%------------------------------------------------------------------
SPM.nscan = IOI.sess_res{s1}.n_frames;
P.name = 'none';
P.h    = 0;
if isempty(name)
    SPM.Sess.U.name = {'Spk'};
else
    SPM.Sess.U.name = {name};
end
SPM.Sess.U.ons = ons;
SPM.Sess.U.dur = dur;
SPM.Sess.U.P = P;
U = spm_get_ons(SPM,1);
% Convolve stimulus functions with basis functions
%------------------------------------------------------------------
[X,Xn,Fc] = spm_Volterra(U,SPM.xBF.bf,volt);

% Resample regressors at acquisition times (32 bin offset)
%-------------------------------------------------
X = X((0:(SPM.nscan - 1))*SPM.xBF.T + SPM.xBF.T0 + 32,:);

% and orthogonalise (within trial type)
%--------------------------------------
for i = 1:length(Fc)
    X(:,Fc(i).i) = spm_orth(X(:,Fc(i).i));
end
end

function K = get_K(X,fwhm,TR)
svec = 1:size(X,1);
K.HParam.type = 'none';
K.LParam.FWHM = fwhm;
K.LParam.type = 'Gaussian';
K = struct( 'HParam', K.HParam,...
    'row', svec ,...
    'RT', TR,...
    'LParam', K.LParam);
K = spm_filter_HPF_LPF_WMDL(K);
K.row = 1:length(svec);
end