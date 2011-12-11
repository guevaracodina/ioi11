function out = ioi_HDM_run(job)
%Based on SPM HDM

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
%select onsets, default is stimulation based
stim_choice=0;
if isfield(job.stim_choice,'electro_stims')
    stim_choice = 1;
end
if isfield(job.stim_choice,'manual_stims')
    stim_choice = 2;
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
%modalities to include
includeHbR = job.includeHbR;
includeHbT = job.includeHbT;
includeFlow = job.includeFlow;
%EM parameters
spm_integrator = job.EM_parameters.spm_integrator;
Niterations = job.EM_parameters.Niterations;
dFcriterion = job.EM_parameters.dFcriterion;
LogAscentRate = job.EM_parameters.LogAscentRate;
Mstep_iterations = job.EM_parameters.Mstep_iterations;
%Baseline choice
if isfield(job.baseline_choice,'baseline_percentile_choice')
    baseline_choice = 1;
    baseline_correction{1} = job.baseline_choice.baseline_percentile_choice.baseline_percentile_HbR/100;
    baseline_correction{2} = job.baseline_choice.baseline_percentile_choice.baseline_percentile_HbT/100;
    baseline_correction{3} = job.baseline_choice.baseline_percentile_choice.baseline_percentile_flow/100;
else
    if isfield(job.baseline_choice,'baseline_offset_choice')
        baseline_choice = 2;
        baseline_correction{1} = job.baseline_choice.baseline_offset_choice.baseline_offset_HbR;
        baseline_correction{2} = job.baseline_choice.baseline_offset_choice.baseline_offset_HbT;
        baseline_correction{3} = job.baseline_choice.baseline_offset_choice.baseline_offset_flow;
    else
        baseline_choice = 0;
    end
end
%save_figures
save_figures = job.save_figures;
generate_figures = job.generate_figures;
show_normalized_parameters = job.show_normalized_parameters;
%Big loop over subjects
for SubjIdx=1:length(job.IOImat)
    try
        tic
        clear IOI ROI HDM
        %Load IOI.mat information
        IOImat = job.IOImat{SubjIdx};
        load(IOImat);
        %Loop over sessions and ROIs
        
        if ~isfield(IOI.res,'ROIOK')
            disp(['No ROI available for subject ' int2str(SubjIdx) ' ... skipping HDM']);
        else
            if ~isfield(IOI.res,'seriesOK')
                disp(['Extracted series not available for subject ' int2str(SubjIdx) ' ... skipping HDM']);
            else
                if ~isfield(IOI.res,'GLMOK')
                    disp(['Onsets (from GLM module) not available for subject ' int2str(SubjIdx) ' ... skipping HDM']);
                else
                    if ~isfield(IOI.res,'HDMOK') || job.force_redo
                        [dir_ioimat dummy] = fileparts(job.IOImat{SubjIdx});
                        %load ROI
                        if ~isempty(job.ROImat)
                            load(job.ROImat{SubjIdx});
                        else
                            try
                                load(IOI.ROI.ROIfname);
                            catch
                                load(fullfile(dir_ioimat,'ROI.mat'));
                            end
                        end
                        if isfield(job.IOImatCopyChoice,'IOImatCopy')
                            newDir = job.IOImatCopyChoice.IOImatCopy.NewIOIdir;
                            newDir = fullfile(dir_ioimat,newDir);
                            if ~exist(newDir,'dir'),mkdir(newDir); end
                            IOImat = fullfile(newDir,'IOI.mat');
                        end
                        [dir1 dummy] = fileparts(IOImat);
                                               
                        PS0 = ioi_get_PS(IOI,includeHbR,includeHbT,includeFlow,job.PhysioModel_Choice);
                        if all_sessions
                            selected_sessions = 1:length(IOI.sess_res);
                        end
                        if all_ROIs
                            selected_ROIs = 1:length(IOI.res.ROI);
                        end
                        %loop over sessions
                        for s1=1:length(IOI.sess_res)
                            if all_sessions || sum(s1==selected_sessions)
                                %loop over ROIs
                                for r1=1:length(IOI.res.ROI)
                                    if all_ROIs || sum(r1==selected_ROIs)
                                        HDM_str = ['S' gen_num_str(s1,2) '_ROI' gen_num_str(r1,3)];
                                        HDMfname = fullfile(dir1,['HDM_' HDM_str '.mat']); 
                                        IOI.HDM{s1,r1}.HDMfname = HDMfname;
                                        HDM0 = [];
                                        Y = [];
                                        %Various options
                                        HDM0.save_figures = save_figures;
                                        HDM0.generate_figures = generate_figures;
                                        HDM0.show_normalized_parameters = show_normalized_parameters;
                                        HDM0.dir1 = dir1;
                                        HDM0.HDM_str = HDM_str;
                                        HDM0.Niterations = Niterations;
                                        HDM0.spm_integrator = spm_integrator;
                                        HDM0.dFcriterion = dFcriterion;
                                        HDM0.LogAscentRate = LogAscentRate;
                                        HDM0.Mstep_iterations = Mstep_iterations;
                                        HDM0.selected_sessions = selected_sessions;
                                        HDM0.selected_ROIs = selected_ROIs;
                                        HDM0.baseline_choice = baseline_choice;
                                        HDM0.baseline_correction = baseline_correction;
                                        %Choose onsets: stimulations or electrophysiology
                                        if electro_stims
                                            U = IOI.Sess(s1).U;
                                        else
                                            name = IOI.sess_res{s1}.names{1};
                                            ons = IOI.sess_res{s1}.onsets{1};
                                            dur = IOI.sess_res{s1}.durations{1};
                                            bases.hrf.derivs = [0 0]; %not used
                                            volt = 1; %not used
                                            [X U] = ioi_get_X(IOI,name,ons,dur,s1,bases,volt);
                                        end
                                        U.u = U.u(33:end); %?
                                        HDM0.U = U;
                                        HDM0.TR = IOI.dev.TR;
                                        HDM0.dt = IOI.dev.TR; %do we need both?!
                                        HDM0.N = 20/IOI.dev.TR; %too large?
                                        %data specification - which modalities to include:
                                        HDM0.PS = PS0;
                                        HDM0.HPF = HPF; %High pass filter on data
                                        HDM0.LPF = LPF;
                                        HDM0=ioi_get_data(ROI,HDM0,r1,s1);
                                        HDM0=ioi_set_physiomodel(HDM0);
                                        %choose priors
                                        HDM0=ioi_set_priors(HDM0);
                                        %setup for priors
                                        HDM0.pE = HDM0.PS.pE;
                                        HDM0.pC = 10*HDM0.PS.pC;
                                        HDM0.pC(end,end) = 10*HDM0.pC(end,end);
                                        warning('off','MATLAB:nearlySingularMatrix');
                                        warning('off','MATLAB:singularMatrix');
                                        % nonlinear system identification
                                        %--------------------------------------------------------------------------
                                        HDM0 = ioi_nlsi(HDM0);
                                        warning('on','MATLAB:nearlySingularMatrix');
                                        warning('on','MATLAB:singularMatrix');
                                        HDM{r1,s1} = HDM0;
                                        save(HDMfname,'HDM');  
                                        IOI.HDM{s1,r1}.Ep = HDM0.Ep;
                                        IOI.HDM{s1,r1}.Cp = HDM0.Cp;
                                        IOI.HDM{s1,r1}.K1 = HDM0.K1;
                                        IOI.HDM{s1,r1}.H1 = HDM0.H1;
                                        IOI.HDM{s1,r1}.PS = HDM0.PS;
                                        IOI.HDM{s1,r1}.F  = HDM0.F;
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