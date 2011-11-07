% -----------------------------------------------------------------------
% Job configuration created by cfg_util (rev $Rev: 4252 $)
% -----------------------------------------------------------------------
matlabbatch{1}.spm.tools.ioi11.msioi1.top_bin_dir = '';
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
matlabbatch{2}.spm.tools.ioi11.conc1.configuration.HbT0 = 0.02;
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
matlabbatch{4}.spm.tools.ioi11.create_roi1.IOImat(1) = cfg_dep;
matlabbatch{4}.spm.tools.ioi11.create_roi1.IOImat(1).tname = 'Select IOI.mat';
matlabbatch{4}.spm.tools.ioi11.create_roi1.IOImat(1).tgt_spec{1}(1).name = 'filter';
matlabbatch{4}.spm.tools.ioi11.create_roi1.IOImat(1).tgt_spec{1}(1).value = 'mat';
matlabbatch{4}.spm.tools.ioi11.create_roi1.IOImat(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{4}.spm.tools.ioi11.create_roi1.IOImat(1).tgt_spec{1}(2).value = 'e';
matlabbatch{4}.spm.tools.ioi11.create_roi1.IOImat(1).sname = 'Compute Flow: IOI.mat';
matlabbatch{4}.spm.tools.ioi11.create_roi1.IOImat(1).src_exbranch = substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{4}.spm.tools.ioi11.create_roi1.IOImat(1).src_output = substruct('.','IOImat');
matlabbatch{4}.spm.tools.ioi11.create_roi1.force_redo = 0;
matlabbatch{4}.spm.tools.ioi11.create_roi1.RemovePreviousROI = 0;
matlabbatch{4}.spm.tools.ioi11.create_roi1.IOImatCopyChoice.IOImatOverwrite = struct([]);
matlabbatch{4}.spm.tools.ioi11.create_roi1.select_names = 1;
matlabbatch{4}.spm.tools.ioi11.create_roi1.AutoROIchoice.ManualROI = struct([]);
matlabbatch{5}.spm.tools.ioi11.extract_roi1.IOImat(1) = cfg_dep;
matlabbatch{5}.spm.tools.ioi11.extract_roi1.IOImat(1).tname = 'Select IOI.mat';
matlabbatch{5}.spm.tools.ioi11.extract_roi1.IOImat(1).tgt_spec{1}(1).name = 'filter';
matlabbatch{5}.spm.tools.ioi11.extract_roi1.IOImat(1).tgt_spec{1}(1).value = 'mat';
matlabbatch{5}.spm.tools.ioi11.extract_roi1.IOImat(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{5}.spm.tools.ioi11.extract_roi1.IOImat(1).tgt_spec{1}(2).value = 'e';
matlabbatch{5}.spm.tools.ioi11.extract_roi1.IOImat(1).sname = 'Create ROI: IOI.mat';
matlabbatch{5}.spm.tools.ioi11.extract_roi1.IOImat(1).src_exbranch = substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{5}.spm.tools.ioi11.extract_roi1.IOImat(1).src_output = substruct('.','IOImat');
matlabbatch{5}.spm.tools.ioi11.extract_roi1.force_redo = 0;
matlabbatch{5}.spm.tools.ioi11.extract_roi1.IOImatCopyChoice.IOImatOverwrite = struct([]);
matlabbatch{5}.spm.tools.ioi11.extract_roi1.session_choice.all_sessions = struct([]);
matlabbatch{5}.spm.tools.ioi11.extract_roi1.ROI_choice.all_ROIs = struct([]);
matlabbatch{6}.spm.tools.ioi11.stim_mean1.IOImat(1) = cfg_dep;
matlabbatch{6}.spm.tools.ioi11.stim_mean1.IOImat(1).tname = 'Select IOI.mat';
matlabbatch{6}.spm.tools.ioi11.stim_mean1.IOImat(1).tgt_spec{1}(1).name = 'filter';
matlabbatch{6}.spm.tools.ioi11.stim_mean1.IOImat(1).tgt_spec{1}(1).value = 'mat';
matlabbatch{6}.spm.tools.ioi11.stim_mean1.IOImat(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{6}.spm.tools.ioi11.stim_mean1.IOImat(1).tgt_spec{1}(2).value = 'e';
matlabbatch{6}.spm.tools.ioi11.stim_mean1.IOImat(1).sname = 'Extract ROI: IOI.mat';
matlabbatch{6}.spm.tools.ioi11.stim_mean1.IOImat(1).src_exbranch = substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{6}.spm.tools.ioi11.stim_mean1.IOImat(1).src_output = substruct('.','IOImat');
matlabbatch{6}.spm.tools.ioi11.stim_mean1.force_redo = 0;
matlabbatch{6}.spm.tools.ioi11.stim_mean1.IOImatCopyChoice.IOImatOverwrite = struct([]);
matlabbatch{6}.spm.tools.ioi11.stim_mean1.stim_choice.default_stims = struct([]);
matlabbatch{6}.spm.tools.ioi11.stim_mean1.session_choice.all_sessions = struct([]);
matlabbatch{6}.spm.tools.ioi11.stim_mean1.ROI_choice.all_ROIs = struct([]);
matlabbatch{6}.spm.tools.ioi11.stim_mean1.window_after = 15;
matlabbatch{6}.spm.tools.ioi11.stim_mean1.window_before = 3;
matlabbatch{6}.spm.tools.ioi11.stim_mean1.normalize_choice = 2;
matlabbatch{6}.spm.tools.ioi11.stim_mean1.generate_figures = 1;
matlabbatch{6}.spm.tools.ioi11.stim_mean1.save_figures = 0;
matlabbatch{6}.spm.tools.ioi11.stim_mean1.add_error_bars = 0;
matlabbatch{7}.spm.tools.ioi11.glm_roi1.IOImat(1) = cfg_dep;
matlabbatch{7}.spm.tools.ioi11.glm_roi1.IOImat(1).tname = 'Select IOI.mat';
matlabbatch{7}.spm.tools.ioi11.glm_roi1.IOImat(1).tgt_spec{1}(1).name = 'filter';
matlabbatch{7}.spm.tools.ioi11.glm_roi1.IOImat(1).tgt_spec{1}(1).value = 'mat';
matlabbatch{7}.spm.tools.ioi11.glm_roi1.IOImat(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{7}.spm.tools.ioi11.glm_roi1.IOImat(1).tgt_spec{1}(2).value = 'e';
matlabbatch{7}.spm.tools.ioi11.glm_roi1.IOImat(1).sname = 'Average stimulations: IOI.mat';
matlabbatch{7}.spm.tools.ioi11.glm_roi1.IOImat(1).src_exbranch = substruct('.','val', '{}',{6}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{7}.spm.tools.ioi11.glm_roi1.IOImat(1).src_output = substruct('.','IOImat');
matlabbatch{7}.spm.tools.ioi11.glm_roi1.force_redo = 0;
matlabbatch{7}.spm.tools.ioi11.glm_roi1.IOImatCopyChoice.IOImatOverwrite = struct([]);
matlabbatch{7}.spm.tools.ioi11.glm_roi1.session_choice.all_sessions = struct([]);
matlabbatch{7}.spm.tools.ioi11.glm_roi1.ROI_choice.all_ROIs = struct([]);
matlabbatch{7}.spm.tools.ioi11.glm_roi1.bases.hrf.derivs = [0 0];
matlabbatch{7}.spm.tools.ioi11.glm_roi1.volt = 1;
matlabbatch{7}.spm.tools.ioi11.glm_roi1.hpf_butter.hpf_butter_On.hpf_butter_freq = 0.05;
matlabbatch{7}.spm.tools.ioi11.glm_roi1.hpf_butter.hpf_butter_On.hpf_butter_order = 3;
matlabbatch{7}.spm.tools.ioi11.glm_roi1.lpf_gauss.fwhm1 = 1;
matlabbatch{7}.spm.tools.ioi11.glm_roi1.stim_choice.default_stims = struct([]);
