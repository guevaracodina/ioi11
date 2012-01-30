%-----------------------------------------------------------------------
% Job configuration created by cfg_util (rev $Rev: 4252 $)
%-----------------------------------------------------------------------
matlabbatch{1}.spm.tools.ioi11.msioi1.top_bin_dir = {'V:\TRI\epiRat\EP\11_11_15,TR03\'};
matlabbatch{1}.spm.tools.ioi11.msioi1.force_redo = 0;
matlabbatch{1}.spm.tools.ioi11.msioi1.configuration_choice.no_shrinkage = struct([]);
matlabbatch{1}.spm.tools.ioi11.msioi1.output_path_choice.output_path_default = struct([]);
matlabbatch{1}.spm.tools.ioi11.msioi1.session_choice.all_sessions = struct([]);
matlabbatch{1}.spm.tools.ioi11.msioi1.save_choice = 2;
matlabbatch{1}.spm.tools.ioi11.msioi1.memmapfileOn = 1;
matlabbatch{1}.spm.tools.ioi11.msioi1.sess_min_image_files = 60;
matlabbatch{1}.spm.tools.ioi11.msioi1.forceProcessingOn = 0;
matlabbatch{2}.spm.tools.ioi11.conc1.IOImat(1) = cfg_dep;
matlabbatch{2}.spm.tools.ioi11.conc1.IOImat(1).tname = 'Select IOI.mat';
matlabbatch{2}.spm.tools.ioi11.conc1.IOImat(1).tgt_spec{1}(1).name = 'filter';
matlabbatch{2}.spm.tools.ioi11.conc1.IOImat(1).tgt_spec{1}(1).value = 'mat';
matlabbatch{2}.spm.tools.ioi11.conc1.IOImat(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{2}.spm.tools.ioi11.conc1.IOImat(1).tgt_spec{1}(2).value = 'e';
matlabbatch{2}.spm.tools.ioi11.conc1.IOImat(1).sname = 'Read Multi-Spectral IOI: IOI.mat';
matlabbatch{2}.spm.tools.ioi11.conc1.IOImat(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{2}.spm.tools.ioi11.conc1.IOImat(1).src_output = substruct('.','IOImat');
matlabbatch{2}.spm.tools.ioi11.conc1.force_redo = 0;
matlabbatch{2}.spm.tools.ioi11.conc1.IOImatCopyChoice.IOImatOverwrite = struct([]);
matlabbatch{2}.spm.tools.ioi11.conc1.configuration.pathlength = 'Dunn';
matlabbatch{2}.spm.tools.ioi11.conc1.configuration.camera_corr = 1;
matlabbatch{2}.spm.tools.ioi11.conc1.configuration.leds_conv = 1;
matlabbatch{2}.spm.tools.ioi11.conc1.configuration.HbT0 = 100;
matlabbatch{2}.spm.tools.ioi11.conc1.MemoryManagementMenu = 1;
matlabbatch{2}.spm.tools.ioi11.conc1.session_choice.all_sessions = struct([]);
matlabbatch{2}.spm.tools.ioi11.conc1.RemoveRGY = 1;
matlabbatch{3}.spm.tools.ioi11.flow1.IOImat(1) = cfg_dep;
matlabbatch{3}.spm.tools.ioi11.flow1.IOImat(1).tname = 'Select IOI.mat';
matlabbatch{3}.spm.tools.ioi11.flow1.IOImat(1).tgt_spec{1}(1).name = 'filter';
matlabbatch{3}.spm.tools.ioi11.flow1.IOImat(1).tgt_spec{1}(1).value = 'mat';
matlabbatch{3}.spm.tools.ioi11.flow1.IOImat(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{3}.spm.tools.ioi11.flow1.IOImat(1).tgt_spec{1}(2).value = 'e';
matlabbatch{3}.spm.tools.ioi11.flow1.IOImat(1).sname = 'Compute Concentrations: IOI.mat';
matlabbatch{3}.spm.tools.ioi11.flow1.IOImat(1).src_exbranch = substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{3}.spm.tools.ioi11.flow1.IOImat(1).src_output = substruct('.','IOImat');
matlabbatch{3}.spm.tools.ioi11.flow1.force_redo = 0;
matlabbatch{3}.spm.tools.ioi11.flow1.IOImatCopyChoice.IOImatOverwrite = struct([]);
matlabbatch{3}.spm.tools.ioi11.flow1.configuration.integ_time = 0.2;
matlabbatch{3}.spm.tools.ioi11.flow1.configuration.window_size = 5;
matlabbatch{3}.spm.tools.ioi11.flow1.session_choice.all_sessions = struct([]);
matlabbatch{3}.spm.tools.ioi11.flow1.RemoveLC = 1;
matlabbatch{4}.spm.tools.ioi11.cine2D1.IOImat(1) = cfg_dep;
matlabbatch{4}.spm.tools.ioi11.cine2D1.IOImat(1).tname = 'Select IOI.mat';
matlabbatch{4}.spm.tools.ioi11.cine2D1.IOImat(1).tgt_spec{1}(1).name = 'filter';
matlabbatch{4}.spm.tools.ioi11.cine2D1.IOImat(1).tgt_spec{1}(1).value = 'mat';
matlabbatch{4}.spm.tools.ioi11.cine2D1.IOImat(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{4}.spm.tools.ioi11.cine2D1.IOImat(1).tgt_spec{1}(2).value = 'e';
matlabbatch{4}.spm.tools.ioi11.cine2D1.IOImat(1).sname = 'Compute Flow: IOI.mat';
matlabbatch{4}.spm.tools.ioi11.cine2D1.IOImat(1).src_exbranch = substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{4}.spm.tools.ioi11.cine2D1.IOImat(1).src_output = substruct('.','IOImat');
matlabbatch{4}.spm.tools.ioi11.cine2D1.force_redo = 0;
matlabbatch{4}.spm.tools.ioi11.cine2D1.IOImatCopyChoice.IOImatOverwrite = struct([]);
matlabbatch{4}.spm.tools.ioi11.cine2D1.session_choice.all_sessions = struct([]);
matlabbatch{4}.spm.tools.ioi11.cine2D1.window_after = 20;
matlabbatch{4}.spm.tools.ioi11.cine2D1.window_before = 2;
matlabbatch{4}.spm.tools.ioi11.cine2D1.normalize_choice = 1;
matlabbatch{4}.spm.tools.ioi11.cine2D1.which_onset_type = 1;
matlabbatch{4}.spm.tools.ioi11.cine2D1.include_flow = 1;
matlabbatch{4}.spm.tools.ioi11.cine2D1.include_OD = 0;
matlabbatch{4}.spm.tools.ioi11.cine2D1.include_HbT = 1;
matlabbatch{4}.spm.tools.ioi11.cine2D1.lpf_choice.lpf_Off = struct([]);
matlabbatch{4}.spm.tools.ioi11.cine2D1.show_movie = 0;
matlabbatch{5}.spm.tools.ioi11.create_onsets1.IOImat(1) = cfg_dep;
matlabbatch{5}.spm.tools.ioi11.create_onsets1.IOImat(1).tname = 'Select IOI.mat';
matlabbatch{5}.spm.tools.ioi11.create_onsets1.IOImat(1).tgt_spec{1}(1).name = 'filter';
matlabbatch{5}.spm.tools.ioi11.create_onsets1.IOImat(1).tgt_spec{1}(1).value = 'mat';
matlabbatch{5}.spm.tools.ioi11.create_onsets1.IOImat(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{5}.spm.tools.ioi11.create_onsets1.IOImat(1).tgt_spec{1}(2).value = 'e';
matlabbatch{5}.spm.tools.ioi11.create_onsets1.IOImat(1).sname = '2D Cine: IOI.mat';
matlabbatch{5}.spm.tools.ioi11.create_onsets1.IOImat(1).src_exbranch = substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{5}.spm.tools.ioi11.create_onsets1.IOImat(1).src_output = substruct('.','IOImat');
matlabbatch{5}.spm.tools.ioi11.create_onsets1.elDir = '';
matlabbatch{5}.spm.tools.ioi11.create_onsets1.force_redo = 0;
matlabbatch{5}.spm.tools.ioi11.create_onsets1.IOImatCopyChoice.IOImatCopy.NewIOIdir = 'Stims';
matlabbatch{5}.spm.tools.ioi11.create_onsets1.session_choice.all_sessions = struct([]);
matlabbatch{5}.spm.tools.ioi11.create_onsets1.stim_choice.default_stims = struct([]);
matlabbatch{6}.spm.tools.ioi11.create_roi1.IOImat(1) = cfg_dep;
matlabbatch{6}.spm.tools.ioi11.create_roi1.IOImat(1).tname = 'Select IOI.mat';
matlabbatch{6}.spm.tools.ioi11.create_roi1.IOImat(1).tgt_spec{1}(1).name = 'filter';
matlabbatch{6}.spm.tools.ioi11.create_roi1.IOImat(1).tgt_spec{1}(1).value = 'mat';
matlabbatch{6}.spm.tools.ioi11.create_roi1.IOImat(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{6}.spm.tools.ioi11.create_roi1.IOImat(1).tgt_spec{1}(2).value = 'e';
matlabbatch{6}.spm.tools.ioi11.create_roi1.IOImat(1).sname = 'Create onsets: IOI.mat';
matlabbatch{6}.spm.tools.ioi11.create_roi1.IOImat(1).src_exbranch = substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{6}.spm.tools.ioi11.create_roi1.IOImat(1).src_output = substruct('.','IOImat');
matlabbatch{6}.spm.tools.ioi11.create_roi1.force_redo = 0;
matlabbatch{6}.spm.tools.ioi11.create_roi1.RemovePreviousROI = 0;
matlabbatch{6}.spm.tools.ioi11.create_roi1.IOImatCopyChoice.IOImatOverwrite = struct([]);
matlabbatch{6}.spm.tools.ioi11.create_roi1.select_names = 1;
matlabbatch{6}.spm.tools.ioi11.create_roi1.AutoROIchoice.ManualROI = struct([]);
matlabbatch{7}.spm.tools.ioi11.extract_roi1.IOImat(1) = cfg_dep;
matlabbatch{7}.spm.tools.ioi11.extract_roi1.IOImat(1).tname = 'Select IOI.mat';
matlabbatch{7}.spm.tools.ioi11.extract_roi1.IOImat(1).tgt_spec{1}(1).name = 'filter';
matlabbatch{7}.spm.tools.ioi11.extract_roi1.IOImat(1).tgt_spec{1}(1).value = 'mat';
matlabbatch{7}.spm.tools.ioi11.extract_roi1.IOImat(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{7}.spm.tools.ioi11.extract_roi1.IOImat(1).tgt_spec{1}(2).value = 'e';
matlabbatch{7}.spm.tools.ioi11.extract_roi1.IOImat(1).sname = 'Create ROI: IOI.mat';
matlabbatch{7}.spm.tools.ioi11.extract_roi1.IOImat(1).src_exbranch = substruct('.','val', '{}',{6}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{7}.spm.tools.ioi11.extract_roi1.IOImat(1).src_output = substruct('.','IOImat');
matlabbatch{7}.spm.tools.ioi11.extract_roi1.force_redo = 0;
matlabbatch{7}.spm.tools.ioi11.extract_roi1.IOImatCopyChoice.IOImatOverwrite = struct([]);
matlabbatch{7}.spm.tools.ioi11.extract_roi1.session_choice.all_sessions = struct([]);
matlabbatch{7}.spm.tools.ioi11.extract_roi1.ROI_choice.all_ROIs = struct([]);
matlabbatch{8}.spm.tools.ioi11.stim_mean1.IOImat(1) = cfg_dep;
matlabbatch{8}.spm.tools.ioi11.stim_mean1.IOImat(1).tname = 'Select IOI.mat';
matlabbatch{8}.spm.tools.ioi11.stim_mean1.IOImat(1).tgt_spec{1}(1).name = 'filter';
matlabbatch{8}.spm.tools.ioi11.stim_mean1.IOImat(1).tgt_spec{1}(1).value = 'mat';
matlabbatch{8}.spm.tools.ioi11.stim_mean1.IOImat(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{8}.spm.tools.ioi11.stim_mean1.IOImat(1).tgt_spec{1}(2).value = 'e';
matlabbatch{8}.spm.tools.ioi11.stim_mean1.IOImat(1).sname = 'Extract ROI: IOI.mat';
matlabbatch{8}.spm.tools.ioi11.stim_mean1.IOImat(1).src_exbranch = substruct('.','val', '{}',{7}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{8}.spm.tools.ioi11.stim_mean1.IOImat(1).src_output = substruct('.','IOImat');
matlabbatch{8}.spm.tools.ioi11.stim_mean1.force_redo = 0;
matlabbatch{8}.spm.tools.ioi11.stim_mean1.IOImatCopyChoice.IOImatOverwrite = struct([]);
matlabbatch{8}.spm.tools.ioi11.stim_mean1.stim_choice.default_stims = struct([]);
matlabbatch{8}.spm.tools.ioi11.stim_mean1.session_choice.all_sessions = struct([]);
matlabbatch{8}.spm.tools.ioi11.stim_mean1.ROI_choice.all_ROIs = struct([]);
matlabbatch{8}.spm.tools.ioi11.stim_mean1.window_after = 20;
matlabbatch{8}.spm.tools.ioi11.stim_mean1.window_before = 2;
matlabbatch{8}.spm.tools.ioi11.stim_mean1.normalize_choice = 1;
matlabbatch{8}.spm.tools.ioi11.stim_mean1.hpf_butter.hpf_butter_On.hpf_butter_freq = 0.01;
matlabbatch{8}.spm.tools.ioi11.stim_mean1.hpf_butter.hpf_butter_On.hpf_butter_order = 3;
matlabbatch{8}.spm.tools.ioi11.stim_mean1.lpf_choice.lpf_Off = struct([]);
matlabbatch{8}.spm.tools.ioi11.stim_mean1.include_flow = 1;
matlabbatch{8}.spm.tools.ioi11.stim_mean1.include_HbT = 1;
matlabbatch{8}.spm.tools.ioi11.stim_mean1.include_OD = 0;
matlabbatch{8}.spm.tools.ioi11.stim_mean1.extract_HRF = 1;
matlabbatch{8}.spm.tools.ioi11.stim_mean1.fit_3_gamma = 0;
matlabbatch{8}.spm.tools.ioi11.stim_mean1.generate_global = 0;
matlabbatch{8}.spm.tools.ioi11.stim_mean1.generate_figures = 0;
matlabbatch{8}.spm.tools.ioi11.stim_mean1.save_figures = 1;
matlabbatch{8}.spm.tools.ioi11.stim_mean1.add_error_bars = 0;
matlabbatch{8}.spm.tools.ioi11.stim_mean1.remove_segment_drift = 0;
matlabbatch{9}.spm.tools.ioi11.glm_roi1.IOImat(1) = cfg_dep;
matlabbatch{9}.spm.tools.ioi11.glm_roi1.IOImat(1).tname = 'Select IOI.mat';
matlabbatch{9}.spm.tools.ioi11.glm_roi1.IOImat(1).tgt_spec{1}(1).name = 'filter';
matlabbatch{9}.spm.tools.ioi11.glm_roi1.IOImat(1).tgt_spec{1}(1).value = 'mat';
matlabbatch{9}.spm.tools.ioi11.glm_roi1.IOImat(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{9}.spm.tools.ioi11.glm_roi1.IOImat(1).tgt_spec{1}(2).value = 'e';
matlabbatch{9}.spm.tools.ioi11.glm_roi1.IOImat(1).sname = 'Average stimulations: IOI.mat';
matlabbatch{9}.spm.tools.ioi11.glm_roi1.IOImat(1).src_exbranch = substruct('.','val', '{}',{8}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{9}.spm.tools.ioi11.glm_roi1.IOImat(1).src_output = substruct('.','IOImat');
matlabbatch{9}.spm.tools.ioi11.glm_roi1.ROImat = '';
matlabbatch{9}.spm.tools.ioi11.glm_roi1.elDir = '';
matlabbatch{9}.spm.tools.ioi11.glm_roi1.force_redo = 0;
matlabbatch{9}.spm.tools.ioi11.glm_roi1.IOImatCopyChoice.IOImatOverwrite = struct([]);
matlabbatch{9}.spm.tools.ioi11.glm_roi1.session_choice.all_sessions = struct([]);
matlabbatch{9}.spm.tools.ioi11.glm_roi1.ROI_choice.all_ROIs = struct([]);
matlabbatch{9}.spm.tools.ioi11.glm_roi1.bases.hrf.derivs = [0 0];
matlabbatch{9}.spm.tools.ioi11.glm_roi1.volt = 1;
matlabbatch{9}.spm.tools.ioi11.glm_roi1.hpf_butter.hpf_butter_On.hpf_butter_freq = 0.01;
matlabbatch{9}.spm.tools.ioi11.glm_roi1.hpf_butter.hpf_butter_On.hpf_butter_order = 3;
matlabbatch{9}.spm.tools.ioi11.glm_roi1.lpf_gauss.fwhm1 = 1;
matlabbatch{9}.spm.tools.ioi11.glm_roi1.stim_choice.default_stims = struct([]);
matlabbatch{9}.spm.tools.ioi11.glm_roi1.generate_figures = 0;
matlabbatch{9}.spm.tools.ioi11.glm_roi1.save_figures = 1;
matlabbatch{9}.spm.tools.ioi11.glm_roi1.figure_show_stim = 1;
matlabbatch{9}.spm.tools.ioi11.glm_roi1.figure_rebase_to_zero_at_stim = 0;
matlabbatch{10}.spm.tools.ioi11.hdm1.IOImat(1) = cfg_dep;
matlabbatch{10}.spm.tools.ioi11.hdm1.IOImat(1).tname = 'Select IOI.mat';
matlabbatch{10}.spm.tools.ioi11.hdm1.IOImat(1).tgt_spec{1}(1).name = 'filter';
matlabbatch{10}.spm.tools.ioi11.hdm1.IOImat(1).tgt_spec{1}(1).value = 'mat';
matlabbatch{10}.spm.tools.ioi11.hdm1.IOImat(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{10}.spm.tools.ioi11.hdm1.IOImat(1).tgt_spec{1}(2).value = 'e';
matlabbatch{10}.spm.tools.ioi11.hdm1.IOImat(1).sname = 'GLM on ROI: IOI.mat';
matlabbatch{10}.spm.tools.ioi11.hdm1.IOImat(1).src_exbranch = substruct('.','val', '{}',{9}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{10}.spm.tools.ioi11.hdm1.IOImat(1).src_output = substruct('.','IOImat');
matlabbatch{10}.spm.tools.ioi11.hdm1.ROImat = '';
matlabbatch{10}.spm.tools.ioi11.hdm1.force_redo = 0;
matlabbatch{10}.spm.tools.ioi11.hdm1.IOImatCopyChoice.IOImatOverwrite = struct([]);
matlabbatch{10}.spm.tools.ioi11.hdm1.session_choice.all_sessions = struct([]);
matlabbatch{10}.spm.tools.ioi11.hdm1.ROI_choice.all_ROIs = struct([]);
matlabbatch{10}.spm.tools.ioi11.hdm1.PhysioModel_Choice = 0;
matlabbatch{10}.spm.tools.ioi11.hdm1.includeHbR = 1;
matlabbatch{10}.spm.tools.ioi11.hdm1.includeHbT = 1;
matlabbatch{10}.spm.tools.ioi11.hdm1.includeFlow = 1;
matlabbatch{10}.spm.tools.ioi11.hdm1.stim_choice.default_stims = struct([]);
matlabbatch{10}.spm.tools.ioi11.hdm1.hpf_butter.hpf_butter_On.hpf_butter_freq = 0.01;
matlabbatch{10}.spm.tools.ioi11.hdm1.hpf_butter.hpf_butter_On.hpf_butter_order = 3;
matlabbatch{10}.spm.tools.ioi11.hdm1.lpf_choice.lpf_gauss_On.fwhm1 = 1;
matlabbatch{10}.spm.tools.ioi11.hdm1.baseline_choice.baseline_percentile_choice.baseline_percentile_HbR = 90;
matlabbatch{10}.spm.tools.ioi11.hdm1.baseline_choice.baseline_percentile_choice.baseline_percentile_HbT = 10;
matlabbatch{10}.spm.tools.ioi11.hdm1.baseline_choice.baseline_percentile_choice.baseline_percentile_flow = 10;
matlabbatch{10}.spm.tools.ioi11.hdm1.EM_parameters.Niterations = 128;
matlabbatch{10}.spm.tools.ioi11.hdm1.EM_parameters.spm_integrator = 'spm_int';
matlabbatch{10}.spm.tools.ioi11.hdm1.EM_parameters.dFcriterion = 1;
matlabbatch{10}.spm.tools.ioi11.hdm1.EM_parameters.LogAscentRate = -2;
matlabbatch{10}.spm.tools.ioi11.hdm1.EM_parameters.Mstep_iterations = 8;
matlabbatch{10}.spm.tools.ioi11.hdm1.show_normalized_parameters = 0;
matlabbatch{10}.spm.tools.ioi11.hdm1.generate_figures = 0;
matlabbatch{10}.spm.tools.ioi11.hdm1.save_figures = 0;
matlabbatch{11}.spm.tools.ioi11.SCKS1.IOImat(1) = cfg_dep;
matlabbatch{11}.spm.tools.ioi11.SCKS1.IOImat(1).tname = 'Select IOI.mat';
matlabbatch{11}.spm.tools.ioi11.SCKS1.IOImat(1).tgt_spec{1}(1).name = 'filter';
matlabbatch{11}.spm.tools.ioi11.SCKS1.IOImat(1).tgt_spec{1}(1).value = 'mat';
matlabbatch{11}.spm.tools.ioi11.SCKS1.IOImat(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{11}.spm.tools.ioi11.SCKS1.IOImat(1).tgt_spec{1}(2).value = 'e';
matlabbatch{11}.spm.tools.ioi11.SCKS1.IOImat(1).sname = 'HDM on ROI: IOI.mat';
matlabbatch{11}.spm.tools.ioi11.SCKS1.IOImat(1).src_exbranch = substruct('.','val', '{}',{10}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{11}.spm.tools.ioi11.SCKS1.IOImat(1).src_output = substruct('.','IOImat');
matlabbatch{11}.spm.tools.ioi11.SCKS1.ROImat = '';
matlabbatch{11}.spm.tools.ioi11.SCKS1.force_redo = 0;
matlabbatch{11}.spm.tools.ioi11.SCKS1.IOImatCopyChoice.IOImatOverwrite = struct([]);
matlabbatch{11}.spm.tools.ioi11.SCKS1.session_choice.all_sessions = struct([]);
matlabbatch{11}.spm.tools.ioi11.SCKS1.ROI_choice.all_ROIs = struct([]);
matlabbatch{11}.spm.tools.ioi11.SCKS1.PhysioModel_Choice = 0;
matlabbatch{11}.spm.tools.ioi11.SCKS1.includeHbR = 1;
matlabbatch{11}.spm.tools.ioi11.SCKS1.includeHbT = 1;
matlabbatch{11}.spm.tools.ioi11.SCKS1.includeFlow = 1;
matlabbatch{11}.spm.tools.ioi11.SCKS1.hpf_butter.hpf_butter_On.hpf_butter_freq = 0.01;
matlabbatch{11}.spm.tools.ioi11.SCKS1.hpf_butter.hpf_butter_On.hpf_butter_order = 3;
matlabbatch{11}.spm.tools.ioi11.SCKS1.SCKSparams.SCKSnoise = 0;
matlabbatch{11}.spm.tools.ioi11.SCKS1.SCKSparams.State_annealing = 0.9995;
matlabbatch{11}.spm.tools.ioi11.SCKS1.SCKSparams.Parameter_annealing = 0.9995;
matlabbatch{11}.spm.tools.ioi11.SCKS1.generate_figures = 0;
matlabbatch{11}.spm.tools.ioi11.SCKS1.save_figures = 0;
matlabbatch{12}.spm.tools.ioi11.ROC1.IOImat(1) = cfg_dep;
matlabbatch{12}.spm.tools.ioi11.ROC1.IOImat(1).tname = 'Select IOI.mat';
matlabbatch{12}.spm.tools.ioi11.ROC1.IOImat(1).tgt_spec{1}(1).name = 'filter';
matlabbatch{12}.spm.tools.ioi11.ROC1.IOImat(1).tgt_spec{1}(1).value = 'mat';
matlabbatch{12}.spm.tools.ioi11.ROC1.IOImat(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{12}.spm.tools.ioi11.ROC1.IOImat(1).tgt_spec{1}(2).value = 'e';
matlabbatch{12}.spm.tools.ioi11.ROC1.IOImat(1).sname = 'SCKS Deconvolution: IOI.mat';
matlabbatch{12}.spm.tools.ioi11.ROC1.IOImat(1).src_exbranch = substruct('.','val', '{}',{11}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{12}.spm.tools.ioi11.ROC1.IOImat(1).src_output = substruct('.','IOImat');
matlabbatch{12}.spm.tools.ioi11.ROC1.ROImat = '';
matlabbatch{12}.spm.tools.ioi11.ROC1.force_redo = 0;
matlabbatch{12}.spm.tools.ioi11.ROC1.IOImatCopyChoice.IOImatOverwrite = struct([]);
matlabbatch{12}.spm.tools.ioi11.ROC1.session_choice.all_sessions = struct([]);
matlabbatch{12}.spm.tools.ioi11.ROC1.ROI_choice.all_ROIs = struct([]);
matlabbatch{12}.spm.tools.ioi11.ROC1.PhysioModel_Choice = 0;
matlabbatch{12}.spm.tools.ioi11.ROC1.includeHbR = 1;
matlabbatch{12}.spm.tools.ioi11.ROC1.includeHbT = 1;
matlabbatch{12}.spm.tools.ioi11.ROC1.includeFlow = 1;
matlabbatch{12}.spm.tools.ioi11.ROC1.hpf_butter.hpf_butter_On.hpf_butter_freq = 0.01;
matlabbatch{12}.spm.tools.ioi11.ROC1.hpf_butter.hpf_butter_On.hpf_butter_order = 3;
matlabbatch{12}.spm.tools.ioi11.ROC1.SCKSparams.SCKSnoise = 0;
matlabbatch{12}.spm.tools.ioi11.ROC1.SCKSparams.State_annealing = 0.9995;
matlabbatch{12}.spm.tools.ioi11.ROC1.SCKSparams.Parameter_annealing = 0.9995;
matlabbatch{12}.spm.tools.ioi11.ROC1.generate_figures = 0;
matlabbatch{12}.spm.tools.ioi11.ROC1.save_figures = 0;