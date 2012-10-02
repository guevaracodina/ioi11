function seizure_mean1 = ioi_seizure_mean_roi_cfg
%_______________________________________________________________________
% Copyright (C) 2011 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%****by Cong on 2012/09/14
IOImat = ioi_dfg_IOImat(1);
redo1 = ioi_dfg_redo(0);
ROImat = ioi_dfg_ROImat(1);
IOImatCopyChoice = ioi_dfg_IOImatCopyChoice('Mean');
ROI_choice = ioi_dfg_ROI_choice;
session_choice = ioi_dfg_session_choice;
[window_before window_after window_offset] = ioi_dfg_window(30,50,60);

% normalize_choice      = cfg_menu;
% normalize_choice.tag  = 'normalize_choice';
% normalize_choice.name = 'Additive normalization choice';
% normalize_choice.labels = {'Median over window before','Time zero','Mean'};
% normalize_choice.values = {1,2,3};
% normalize_choice.val  = {1};
% normalize_choice.help = {'Normalization choice. In one test,'
%     'The mean standard deviation was higher by 10% or more'
%     'when using time zero as the baseline, compared to taking '
%     'an average (median) over the window before.'}';

% hpf_butter = ioi_dfg_hpf_butter(1,0.01,3);
lpf_choice = ioi_dfg_lpf_choice(0,0.5); % mice the FWHM is 0.5 but for rat the FWHM is 0.67
% 
% mult_normalize_choice      = cfg_menu;
% mult_normalize_choice.tag  = 'mult_normalize_choice';
% mult_normalize_choice.name = 'Multiplicative normalization choice';
% mult_normalize_choice.labels = {'Divide by constant','Divide by baseline of individual seizure'};
% mult_normalize_choice.values = {1,0};
% mult_normalize_choice.val  = {0};
% mult_normalize_choice.help = {'For oxy, deoxy and total hemoglobin only, '
%     'choose to divide by a constant (e.g. 40, 60 or 100 microMolar'
%     'or to divide by a different value for each seizure, '
%     'namely its baseline value in microMolar,'
%     'however, this value may be too volatile and give nonsensical results.'};
% 
% generate_global      = cfg_menu;
% generate_global.tag  = 'generate_global';
% generate_global.name = 'Generate global data';
% generate_global.labels = {'Yes','No'};
% generate_global.values = {1,0};
% generate_global.val  = {0};
% generate_global.help = {'Generate data averaged over all sessions.'}';

add_error_bars      = cfg_menu;
add_error_bars.tag  = 'add_error_bars';
add_error_bars.name = 'Add error bars';
add_error_bars.labels = {'Yes','No'};
add_error_bars.values = {1,0};
add_error_bars.val  = {0};
add_error_bars.help = {'Add error bars.'}';


std_choice         = cfg_entry; 
std_choice.name    = 'Choice of standard deviation threshold';
std_choice.tag     = 'std_choice';       
std_choice.strtype = 'r';
std_choice.num     = [1 1];     
std_choice.val     = {1};
std_choice.help    = {'Choice of standard deviation threshold.'}';

[generate_figures save_figures] = ioi_dfg_generate_figures;
IC = ioi_dfg_include_colors(0,1,1,1,0);

% Executable Branch
seizure_mean1      = cfg_exbranch;       % This is the branch that has information about how to run this module
seizure_mean1.name = 'Average seizures (on ROIs)';             % The display name
seizure_mean1.tag  = 'seizure_mean1'; %Very important: tag is used when calling for execution
seizure_mean1.val  = {IOImat ROImat redo1 IOImatCopyChoice session_choice ...
    ROI_choice window_after window_before window_offset   ...
    lpf_choice IC  ...
    generate_figures save_figures ...
    add_error_bars};    % The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs
seizure_mean1.prog = @ioi_seizure_mean_roi_run;  % A function handle that will be called with the harvested job to run the computation
seizure_mean1.vout = @ioi_cfg_vout_seizure_mean; % A function handle that will be called with the harvested job to determine virtual outputs
seizure_mean1.help = {'Calculate average over seizures.'};
return

%make IOI.mat available as a dependency
function vout = ioi_cfg_vout_seizure_mean(job)
vout = cfg_dep;                     % The dependency object
vout.sname      = 'IOI.mat';       % Displayed dependency name
vout.src_output = substruct('.','IOImat'); %{1}); %,'IOImat');
%substruct('()',{1}); % The output subscript reference. This could be any reference into the output variable created during computation
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});