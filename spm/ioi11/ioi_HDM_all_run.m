function out = ioi_HDM_all_run(job)
%Based on SPM HDM
%select a subset of sessions
if isfield(job.session_choice,'select_sessions')
    O.all_sessions = 0;
    O.selected_sessions = job.session_choice.select_sessions.selected_sessions;
else
    O.all_sessions = 1;
end
%find data selection mode
if isfield(job.data_selection_choice,'ROI_mode')
    O.image_mode = 0;
    O.ROImat = job.data_selection_choice.ROI_mode.ROImat;
    %select a subset of ROIs
    if isfield(job.data_selection_choice.ROI_mode.ROI_choice,'select_ROIs')
        O.all_ROIs = 0;
        O.selected_ROIs = job.data_selection_choice.ROI_mode.ROI_choice.select_ROIs.selected_ROIs;
    else
        O.all_ROIs = 1;
    end
    DO.show_mse = job.data_selection_choice.ROI_mode.show_mse;
else
    O.image_mode = 1;
end

%HPF
if isfield(job.hpf_butter,'hpf_butter_On')
    HPF.hpf_butter_On = 1;
    HPF.hpf_butter_freq = job.hpf_butter.hpf_butter_On.hpf_butter_freq;
    HPF.hpf_butter_order = job.hpf_butter.hpf_butter_On.hpf_butter_order;
else
    HPF.hpf_butter_On = 0;
end
%LPF
if isfield(job.lpf_choice,'lpf_gauss_On')
    LPF.lpf_gauss_On = 1;
    LPF.fwhm1 = job.lpf_choice.lpf_gauss_On.fwhm1;
else
    LPF.lpf_gauss_On = 0;
end
%Hemodynamic model choice
O.PhysioModel_Choice = job.PhysioModel_Choice;
%modalities to include
O.includeHbR = job.includeHbR;
O.includeHbT = job.includeHbT;
O.includeFlow = job.includeFlow;
%EM parameters
EM.spm_integrator = job.EM_parameters.spm_integrator;
EM.Niterations = job.EM_parameters.Niterations;
EM.dFcriterion = job.EM_parameters.dFcriterion;
EM.LogAscentRate = job.EM_parameters.LogAscentRate;
EM.Mstep_iterations = job.EM_parameters.Mstep_iterations;
%Onset amplitudes
O.use_onset_amplitudes = job.use_onset_amplitudes;
%Perform HDM (only_display = 0) or only display previous results (=1)
O.only_display = job.only_display;
%Baseline choice
if isfield(job.baseline_choice,'baseline_percentile_choice')
    O.baseline_choice = 1;
    O.baseline_correction{1} = job.baseline_choice.baseline_percentile_choice.baseline_percentile_HbR/100;
    O.baseline_correction{2} = job.baseline_choice.baseline_percentile_choice.baseline_percentile_HbT/100;
    O.baseline_correction{3} = job.baseline_choice.baseline_percentile_choice.baseline_percentile_flow/100;
else
    if isfield(job.baseline_choice,'baseline_offset_choice')
        O.baseline_choice = 2;
        O.baseline_correction{1} = job.baseline_choice.baseline_offset_choice.baseline_offset_HbR;
        O.baseline_correction{2} = job.baseline_choice.baseline_offset_choice.baseline_offset_HbT;
        O.baseline_correction{3} = job.baseline_choice.baseline_offset_choice.baseline_offset_flow;
    else
        O.baseline_choice = 0;
        O.baseline_correction = [];
    end
end
%Display options
DO.plot_algebraic_CMRO2 = job.plot_algebraic_CMRO2;
DO.save_figures = job.save_figures;
DO.generate_figures = job.generate_figures;
DO.show_normalized_parameters = job.show_normalized_parameters;
%Simulations
if isfield(job.simuOn,'simuYes')
    S.simuS     = job.simuOn.simuYes.simuS; %Stimuli types to include
    S.simuIt    = job.simuOn.simuYes.simuIt; %Number of random iterations
    S.simuA     = job.simuOn.simuYes.simuA; %Signal amplitude, as % of BOLD signal
    S.simuP     = job.simuOn.simuYes.simuP; %Parameters to vary
    S.simuPrior = job.simuOn.simuYes.simuPrior; %Priors to use
    S.simuR     = job.simuOn.simuYes.simuR; %Range to sample
    S.simuUpsample = job.simuOn.simuYes.simuUpsample; %Upsampling factor on data
    S.simuNoise = job.simuOn.simuYes.simuNoise; %Yes to include background noise based on restscans
    S.simuOn = 1;
else
    S.simuOn = 0;
end
%Big loop over subjects
for SubjIdx=1:length(job.IOImat)
    try
        tic
        clear IOI ROI HDM
        %Load IOI.mat information
        IOImat = job.IOImat{SubjIdx};
        load(IOImat);
        %Loop over sessions and ROIs
        
        if ~isfield(IOI.res,'ROIOK') && ~O.image_mode
            disp(['No ROI available for subject ' int2str(SubjIdx) ' ... skipping HDM']);
        else
            if ~isfield(IOI.res,'seriesOK') && ~O.image_mode
                disp(['Extracted series not available for subject ' int2str(SubjIdx) ' ... skipping HDM']);
            else
                %if ~isfield(IOI.res,'OnsetsOK')
                %    disp(['Onsets (from GLM module) not available for subject ' int2str(SubjIdx) ' ... skipping HDM']);
                %else
                if ~isfield(IOI.res,'HDMOK') || job.force_redo
                    [dir_ioimat dummy] = fileparts(job.IOImat{SubjIdx});
                    if ~O.image_mode
                        %load ROI
                        if ~isempty(O.ROImat)
                            load(O.ROImat{SubjIdx});
                        else
                            try
                                load(IOI.ROI.ROIfname);
                            catch
                                load(fullfile(dir_ioimat,'ROI.mat'));
                            end
                        end
                        if O.all_ROIs
                            O.selected_ROIs = 1:length(IOI.res.ROI);
                        end
                    end
                    if isfield(job.IOImatCopyChoice,'IOImatCopy')
                        newDir = job.IOImatCopyChoice.IOImatCopy.NewIOIdir;
                        newDir = fullfile(dir_ioimat,newDir);
                        if ~exist(newDir,'dir'),mkdir(newDir); end
                        IOImat = fullfile(newDir,'IOI.mat');
                    end
                    [dir1 dummy] = fileparts(IOImat);
                    
                    O = ioi_get_PS(IOI,O);
                    if O.all_sessions
                        O.selected_sessions = 1:length(IOI.sess_res);
                    end
                    
                    %loop over sessions
                    for s1=1:length(IOI.sess_res)
                        if O.all_sessions || sum(s1==O.selected_sessions)
                            if O.image_mode
                                HDM_str = ['S' gen_num_str(s1,2)];
                                HDMfname = fullfile(dir1,['HDM_' HDM_str '.mat']);
                                IOI.HDM{s1}.HDMfname = HDMfname;
                                if ~O.only_display
                                    HDM0 = [];
                                    %Various options
                                    HDM0.O = O;
                                    HDM0.DO = DO;
                                    HDM0.EM = EM;
                                    HDM0.dir1 = dir1;
                                    HDM0.HDM_str = HDM_str;
                                    %Simulation option
                                    HDM0.S = S;
                                    %onsets: they are now specified in the module create_onsets
                                    name = IOI.sess_res{s1}.names{1};
                                    ons = IOI.sess_res{s1}.onsets{1};
                                    dur = IOI.sess_res{s1}.durations{1};
                                    if O.use_onset_amplitudes
                                        amp = IOI.sess_res{s1}.parameters{1};
                                        %normalize -- but what if amp takes extreme values?
                                        amp = amp/mean(amp);
                                    else
                                        amp = [];
                                    end
                                    bases.hrf.derivs = [0 0]; %not used
                                    [dummyX U] = ioi_get_X(IOI,name,ons,dur,amp,s1,bases,1); %call only to get U
                                    U.u = U.u(33:end); %?
                                    HDM0.U = U;
                                    HDM0.dt = IOI.dev.TR;
                                    HDM0.N = 20/IOI.dev.TR; %too large?
                                    %Filters
                                    HDM0.HPF = HPF; %High pass filter on data
                                    HDM0.LPF = LPF;
                                    %data specification - which modalities to include:
                                    HDM0=ioi_get_data(ROI,HDM0,r1,s1);
                                    HDM0=ioi_set_physiomodel(HDM0);
                                    %choose priors
                                    HDM0=ioi_set_priors(HDM0);
                                    %scaling factors on covariance ... do we need that?
                                    %HDM0.pC = 10*HDM0.pC;
                                    %HDM0.pC(end,end) = 10*HDM0.pC(end,end);
                                    if S.simuOn
                                        HDM0 = ioi_simu_gen_parameters(HDM0);
                                        simuIt = HDM0.S.simuIt;
                                        HDM0.Y0 = HDM0.Y; %background
                                        HDM0.HDM_str0 = HDM0.HDM_str; %save each figure of fit
                                        SHDM = []; %simulation HDM structure
                                    else
                                        simuIt = 1;
                                    end
                                    warning('off','MATLAB:nearlySingularMatrix');
                                    warning('off','MATLAB:singularMatrix');
                                    % big loop over simulation iterations
                                    %--------------------------------------------------------------------------
                                    for it1=1:simuIt
                                        if S.simuOn
                                            HDM0 = ioi_set_simu(HDM0,it1);
                                        end
                                        % nonlinear system identification
                                        %--------------------------------------------------------------------------
                                        HDM0 = ioi_nlsi(HDM0);
                                        if S.simuOn
                                            SHDM{it1} = ioi_save_simu(HDM0,it1);
                                        end
                                    end
                                    warning('on','MATLAB:nearlySingularMatrix');
                                    warning('on','MATLAB:singularMatrix');
                                else
                                    try
                                        load(HDMfname);
                                        HDM0 = HDM{s1};
                                    catch
                                        disp(['HDM.mat not found for Session ' int2str(s1) ' and ROI ' int2str(r1)]);
                                    end
                                end
                                if ~S.simuOn
                                    HDM{s1} = HDM0;
                                    save(HDMfname,'HDM');
                                    ioi_HDM_display(HDM0);
                                    
                                    %Store some of the information in IOI structure
                                    IOI.HDM{s1}.Ep = HDM0.Ep;
                                    IOI.HDM{s1}.Cp = HDM0.Cp;
                                    IOI.HDM{s1}.K1 = HDM0.K1;
                                    IOI.HDM{s1}.H1 = HDM0.H1;
                                    IOI.HDM{s1}.F  = HDM0.F;
                                else
                                    save(HDMfname,'SHDM');
                                end
                            else
                                %On ROIs
                                HDM_str = ['S' gen_num_str(s1,2)];
                                HDMfname = fullfile(dir1,['HDM_' HDM_str '.mat']);
                                IOI.HDM{s1}.HDMfname = HDMfname;
                                if ~O.only_display
                                    HDM0 = [];
                                    %Various options
                                    HDM0.O = O;
                                    HDM0.DO = DO;
                                    HDM0.EM = EM;
                                    HDM0.dir1 = dir1;
                                    HDM0.HDM_str = HDM_str;
                                    %Simulation option
                                    HDM0.S = S;
                                    %onsets: they are now specified in the module create_onsets
                                    name = IOI.sess_res{s1}.names{1};
                                    ons = IOI.sess_res{s1}.onsets{1};
                                    dur = IOI.sess_res{s1}.durations{1};
                                    if O.use_onset_amplitudes
                                        amp = IOI.sess_res{s1}.parameters{1};
                                        %normalize -- but what if amp takes extreme values?
                                        amp = amp/mean(amp);
                                    else
                                        amp = [];
                                    end
                                    bases.hrf.derivs = [0 0]; %not used
                                    [dummyX U] = ioi_get_X(IOI,name,ons,dur,amp,s1,bases,1); %call only to get U
                                    U.u = U.u(33:end); %?
                                    HDM0.U = U;
                                    HDM0.dt = IOI.dev.TR;
                                    HDM0.N = 20/IOI.dev.TR; %too large?
                                    %Filters
                                    HDM0.HPF = HPF; %High pass filter on data
                                    HDM0.LPF = LPF;
                                    %data specification - which modalities to include:
                                    %loop over ROIs
                                    ct = 0;
                                    for r1=1:length(IOI.res.ROI)
                                        if O.all_ROIs || sum(r1==O.selected_ROIs)
                                            ct = ct+1;
                                            HDM0=ioi_get_data(ROI,HDM0,r1,s1);
                                            Yt = HDM0.Y.y;
                                            if ct == 1
                                                Y = zeros(size(Yt,1),size(Yt,2),length(O.selected_ROIs));
                                            end
                                            Y(:,:,ct) = Yt;
                                        end
                                    end
                                    Y = permute(Y,[1 3 2]);
                                    Y = reshape(Y,size(Yt,1),size(Yt,2)*length(O.selected_ROIs));
                                    HDM0.Y.y = Y;
                                    HDM0=ioi_set_physiomodel(HDM0);
                                    %include number of ROIs in number of outputs
                                    HDM0.l = HDM0.l*length(O.selected_ROIs);
                                    %choose priors
                                    HDM0=ioi_set_priors(HDM0);
                                    %scaling factors on covariance ... do we need that?
                                    %HDM0.pC = 10*HDM0.pC;
                                    %HDM0.pC(end,end) = 10*HDM0.pC(end,end);
                                    if S.simuOn
                                        HDM0 = ioi_simu_gen_parameters(HDM0);
                                        simuIt = HDM0.S.simuIt;
                                        HDM0.Y0 = HDM0.Y; %background
                                        HDM0.HDM_str0 = HDM0.HDM_str; %save each figure of fit
                                        SHDM = []; %simulation HDM structure
                                    else
                                        simuIt = 1;
                                    end
                                    warning('off','MATLAB:nearlySingularMatrix');
                                    warning('off','MATLAB:singularMatrix');
                                    % big loop over simulation iterations
                                    %--------------------------------------------------------------------------
                                    for it1=1:simuIt
                                        if S.simuOn
                                            HDM0 = ioi_set_simu(HDM0,it1);
                                        end
                                        % nonlinear system identification
                                        %--------------------------------------------------------------------------
                                        HDM0 = ioi_nlsi(HDM0);
                                        if S.simuOn
                                            SHDM{it1} = ioi_save_simu(HDM0,it1);
                                        end
                                    end
                                    warning('on','MATLAB:nearlySingularMatrix');
                                    warning('on','MATLAB:singularMatrix');
                                else
                                    try
                                        load(HDMfname);
                                        HDM0 = HDM{s1,r1};
                                    catch
                                        disp(['HDM.mat not found for Session ' int2str(s1) ' and ROI ' int2str(r1)]);
                                    end
                                end
                                if ~S.simuOn
                                    HDM{s1,r1} = HDM0;
                                    save(HDMfname,'HDM');
                                    ioi_HDM_display(HDM0);
                                    
                                    %Store some of the information in IOI structure
                                    IOI.HDM{s1,r1}.Ep = HDM0.Ep;
                                    IOI.HDM{s1,r1}.Cp = HDM0.Cp;
                                    IOI.HDM{s1,r1}.K1 = HDM0.K1;
                                    IOI.HDM{s1,r1}.H1 = HDM0.H1;
                                    IOI.HDM{s1,r1}.F  = HDM0.F;
                                else
                                    save(HDMfname,'SHDM');
                                end
                                
                            end
                            disp(['HDM for session ' int2str(s1) ' completed']);
                        end
                    end
                    IOI.res.HDMOK = 1;
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
        out.IOImat{SubjIdx} = job.IOImat{SubjIdx};
    end
end
end