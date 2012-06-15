function glm_roi1 = ioi_GLM_on_ROI_cfg
%_______________________________________________________________________
% Copyright (C) 2011 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal

IOImat = ioi_dfg_IOImat(1);
redo1 = ioi_dfg_redo(0);
ROImat = ioi_dfg_ROImat(1);
IOImatCopyChoice = ioi_dfg_IOImatCopyChoice('GLM');

%%%%%%%%%%%%%%%%%%%%%%
ROI_choice = ioi_dfg_ROI_choice;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
session_choice = ioi_dfg_session_choice;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------------------------------------------------------------------
% derivs Model derivatives
% ---------------------------------------------------------------------
derivs         = cfg_menu;
derivs.tag     = 'derivs';
derivs.name    = 'Model derivatives';
derivs.help    = {'Model HRF Derivatives. The canonical HRF combined with time and dispersion derivatives comprise an ''informed'' basis set, as the shape of the canonical response conforms to the hemodynamic response that is commonly observed. The incorporation of the derivate terms allow for variations in subject-to-subject and voxel-to-voxel responses. The time derivative allows the peak response to vary by plus or minus a second and the dispersion derivative allows the width of the response to vary. The informed basis set requires an SPM{F} for inference. T-contrasts over just the canonical are perfectly valid but assume constant delay/dispersion. The informed basis set compares favourably with eg. FIR bases on many data sets. '};
derivs.labels = {
                 'No derivatives'
                 'Time derivatives'
                 'Time and Dispersion derivatives'
}';
derivs.values = {[0 0] [1 0] [1 1]};
derivs.val    = {[0 0]};

%Additional HRFs

HRF_ROI      = cfg_entry;
HRF_ROI.tag  = 'HRF_ROI';
HRF_ROI.name = 'Enter ROI for desired HRF';
HRF_ROI.strtype  = 'r';
HRF_ROI.num = [1 1];
HRF_ROI.val{1} = 1;
HRF_ROI.help = {'Enter ROI for desired HRF.'};

HRF_global         = cfg_branch;
HRF_global.tag     = 'HRF_global';
HRF_global.name    = 'HRF calculated from all sessions';
HRF_global.val     = {};
HRF_global.help    = {'Use global calculation of HRF obtained from all '
    'previously specified sessions in averaging module.'}';

HRF_selected_session      = cfg_entry;
HRF_selected_session.tag  = 'HRF_selected_session';
HRF_selected_session.name = 'Enter session where HRF was calculated';
HRF_selected_session.strtype  = 'r';
HRF_selected_session.num = [1 1];
HRF_selected_session.val{1} = 1;
HRF_selected_session.help = {'Enter session where HRF was calculated.'};

HRF_select_session         = cfg_branch;
HRF_select_session.tag     = 'HRF_select_session';
HRF_select_session.name    = 'Select session where HRF was calculated';
HRF_select_session.val     = {HRF_selected_session};
HRF_select_session.help    = {'Select session where HRF was calculated.'};

HRF_session_choice        = cfg_choice;
HRF_session_choice.name   = 'Choose session selection method';
HRF_session_choice.tag    = 'HRF_session_choice';
HRF_session_choice.values = {HRF_global,HRF_select_session};
HRF_session_choice.val    = {HRF_select_session};
HRF_session_choice.help   = {'Choose session selection method'}';

%%%%%%%%%%%%%

HRF_respective         = cfg_branch;
HRF_respective.tag     = 'HRF_respective';
HRF_respective.name    = 'Use HRF calculated over same chromophore as the data to do the GLM over';
HRF_respective.val     = {};
HRF_respective.help    = {'Use HRF calculated over same chromophore as the data to do the GLM over'}';

HRF_selected_chromophore         = cfg_menu;
HRF_selected_chromophore.tag     = 'HRF_selected_chromophore';
HRF_selected_chromophore.name    = 'Chromophore of HRF to select';
HRF_selected_chromophore.help    = {'Select chromophore for the HRF to be used on all data, irrespective of their chromophore'};
HRF_selected_chromophore.labels = {
                 'HbR'
                 'HbO'
                 'HbT'
                 'Flow'
}';
HRF_selected_chromophore.values = {0 1 2 3};
HRF_selected_chromophore.val    = {0};

HRF_select_chromophore         = cfg_branch;
HRF_select_chromophore.tag     = 'HRF_select_chromophore';
HRF_select_chromophore.name    = 'Select chromophore corresponding to desired HRF';
HRF_select_chromophore.val     = {HRF_selected_chromophore};
HRF_select_chromophore.help    = {'Select chromophore corresponding to desired HRF.'};

HRF_chromophore_choice        = cfg_choice;
HRF_chromophore_choice.name   = 'Choose chromophore selection method';
HRF_chromophore_choice.tag    = 'HRF_chromophore_choice';
HRF_chromophore_choice.values = {HRF_respective,HRF_select_chromophore};
HRF_chromophore_choice.val    = {HRF_respective};
HRF_chromophore_choice.help   = {'Choose session selection method'}';

HRF_selected_stimulus      = cfg_entry;
HRF_selected_stimulus.tag  = 'HRF_selected_stimulus';
HRF_selected_stimulus.name = 'Enter stimulus number (onset type)  where HRF was calculated';
HRF_selected_stimulus.strtype  = 'r';
HRF_selected_stimulus.num = [1 1];
HRF_selected_stimulus.val{1} = 1;
HRF_selected_stimulus.help = {'Enter stimulus number (onset type) where HRF was calculated.'};

use_onset_amplitudes      = cfg_menu;
use_onset_amplitudes.tag  = 'use_onset_amplitudes';
use_onset_amplitudes.name = 'Use onset amplitudes to weigh the hemodynamic response';
use_onset_amplitudes.labels = {'Yes','No'};
use_onset_amplitudes.values = {1,0};
use_onset_amplitudes.val  = {0};
use_onset_amplitudes.help = {'Use onset amplitudes as parameters to weigh the hemodynamic response.'}';

% ---------------------------------------------------------------------
% hrf Canonical HRF
% ---------------------------------------------------------------------
specific_nlinfit         = cfg_branch;
specific_nlinfit.tag     = 'specific_nlinfit'; %note tag "rat" is historical and no longer applies
specific_nlinfit.name    = 'animal-specific: nlinfit';
specific_nlinfit.val     = {HRF_ROI HRF_session_choice HRF_chromophore_choice ...
    HRF_selected_stimulus};
specific_nlinfit.help    = {'Animal specific Hemodynamic Response Function.'
    'Calculated with Matlab nlinfit'}';

specific_EM         = cfg_branch;
specific_EM.tag     = 'specific_EM'; %note tag "mouse" is historical and no longer applies
specific_EM.name    = 'animal-specific: EM';
specific_EM.val     = {HRF_ROI HRF_session_choice HRF_chromophore_choice ...
    HRF_selected_stimulus};
specific_EM.help    = {'Animal specific Hemodynamic Response Function.'
    'Calculated with EM (Expectation Maximization) - The preferred choice'}';

%%%%%%%%%%%%%%% HRF %%%%%%%%%%%%%%%%%%%%%%%%

% ---------------------------------------------------------------------
% hrf Canonical HRF
% ---------------------------------------------------------------------
hrf         = cfg_branch;
hrf.tag     = 'hrf';
hrf.name    = 'Canonical HRF';
hrf.val     = {derivs };
hrf.help    = {'Canonical Hemodynamic Response Function. This is the default option. Contrasts of these effects have a physical interpretation and represent a parsimonious way of characterising event-related responses. This option is also useful if you wish to look separately at activations and deactivations (this is implemented using a t-contrast with a +1 or -1 entry over the canonical regressor). '};
% ---------------------------------------------------------------------
% length Window length
% ---------------------------------------------------------------------
length         = cfg_entry;
length.tag     = 'length';
length.name    = 'Window length';
length.help    = {'Post-stimulus window length (in seconds)'};
length.strtype = 'e';
length.num     = [1 1];
% ---------------------------------------------------------------------
% order Order
% ---------------------------------------------------------------------
order         = cfg_entry;
order.tag     = 'order';
order.name    = 'Order';
order.help    = {'Number of basis functions'};
order.strtype = 'e';
order.num     = [1 1];
% ---------------------------------------------------------------------
% fourier Fourier Set
% ---------------------------------------------------------------------
fourier         = cfg_branch;
fourier.tag     = 'fourier';
fourier.name    = 'Fourier Set';
fourier.val     = {length order };
fourier.help    = {'Fourier basis functions. This option requires an SPM{F} for inference.'};
% ---------------------------------------------------------------------
% length Window length
% ---------------------------------------------------------------------
length         = cfg_entry;
length.tag     = 'length';
length.name    = 'Window length';
length.help    = {'Post-stimulus window length (in seconds)'};
length.strtype = 'e';
length.num     = [1 1];
% ---------------------------------------------------------------------
% order Order
% ---------------------------------------------------------------------
order         = cfg_entry;
order.tag     = 'order';
order.name    = 'Order';
order.help    = {'Number of basis functions'};
order.strtype = 'e';
order.num     = [1 1];
% ---------------------------------------------------------------------
% fourier_han Fourier Set (Hanning)
% ---------------------------------------------------------------------
fourier_han         = cfg_branch;
fourier_han.tag     = 'fourier_han';
fourier_han.name    = 'Fourier Set (Hanning)';
fourier_han.val     = {length order };
fourier_han.help    = {'Fourier basis functions with Hanning Window - requires SPM{F} for inference.'};
% ---------------------------------------------------------------------
% length Window length
% ---------------------------------------------------------------------
length         = cfg_entry;
length.tag     = 'length';
length.name    = 'Window length';
length.help    = {'Post-stimulus window length (in seconds)'};
length.strtype = 'e';
length.num     = [1 1];
% ---------------------------------------------------------------------
% order Order
% ---------------------------------------------------------------------
order         = cfg_entry;
order.tag     = 'order';
order.name    = 'Order';
order.help    = {'Number of basis functions'};
order.strtype = 'e';
order.num     = [1 1];
% ---------------------------------------------------------------------
% gamma Gamma Functions
% ---------------------------------------------------------------------
gamma         = cfg_branch;
gamma.tag     = 'gamma';
gamma.name    = 'Gamma Functions';
gamma.val     = {length order };
gamma.help    = {'Gamma basis functions - requires SPM{F} for inference.'};
% ---------------------------------------------------------------------
% length Window length
% ---------------------------------------------------------------------
length         = cfg_entry;
length.tag     = 'length';
length.name    = 'Window length';
length.help    = {'Post-stimulus window length (in seconds)'};
length.strtype = 'e';
length.num     = [1 1];
% ---------------------------------------------------------------------
% order Order
% ---------------------------------------------------------------------
order         = cfg_entry;
order.tag     = 'order';
order.name    = 'Order';
order.help    = {'Number of basis functions'};
order.strtype = 'e';
order.num     = [1 1];
% ---------------------------------------------------------------------
% fir Finite Impulse Response
% ---------------------------------------------------------------------
fir         = cfg_branch;
fir.tag     = 'fir';
fir.name    = 'Finite Impulse Response';
fir.val     = {length order };
fir.help    = {'Finite impulse response - requires SPM{F} for inference.'};
% ---------------------------------------------------------------------
% bases Basis Functions
% ---------------------------------------------------------------------
bases         = cfg_choice;
bases.tag     = 'bases';
bases.name    = 'Basis Functions';
bases.val     = {hrf };
bases.help    = {'The most common choice of basis function is the Canonical HRF with or without time and dispersion derivatives. '};
bases.values  = {hrf specific_nlinfit specific_EM fourier fourier_han gamma fir };
% ---------------------------------------------------------------------
% volt Model Interactions (Volterra)
% ---------------------------------------------------------------------
volt         = cfg_menu;
volt.tag     = 'volt';
volt.name    = 'Model Interactions (Volterra)';
volt.help    = {
                'Generalized convolution of inputs (U) with basis set (bf).'
                ''
                'For first order expansions the causes are simply convolved (e.g. stick functions) in U.u by the basis functions in bf to create a design matrix X.  For second order expansions new entries appear in ind, bf and name that correspond to the interaction among the orginal causes. The basis functions for these efects are two dimensional and are used to assemble the second order kernel. Second order effects are computed for only the first column of U.u.'
                'Interactions or response modulations can enter at two levels.  Firstly the stick function itself can be modulated by some parametric variate (this can be time or some trial-specific variate like reaction time) modeling the interaction between the trial and the variate or, secondly interactions among the trials themselves can be modeled using a Volterra series formulation that accommodates interactions over time (and therefore within and between trial types).'
}';
volt.labels = {
               'Do not model Interactions'
               'Model Interactions'
}';
volt.values = {1 2};
volt.val    = {1};

hpf_butter = ioi_dfg_hpf_butter(1,0.01,3);

% ---------------------------------------------------------------------
% lpf Low-pass filter
% ---------------------------------------------------------------------
fwhm1      = cfg_entry;
fwhm1.tag  = 'fwhm1';
fwhm1.name = 'FWHM in seconds';
fwhm1.val = {0.67};
fwhm1.strtype = 'r';  
fwhm1.num     = [1 1]; 
fwhm1.help    = {'FWHM in seconds.'}; 

lpf_gauss         = cfg_branch;
lpf_gauss.tag     = 'lpf_gauss';
lpf_gauss.name    = 'Gaussian Filter';
lpf_gauss.val     = {fwhm1}; 
lpf_gauss.help    = {'Specify properties of Gaussian filter'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IC = ioi_dfg_include_colors(0,1,1,1,0);
[generate_figures save_figures] = ioi_dfg_generate_figures;

figure_show_stim      = cfg_menu;
figure_show_stim.tag  = 'figure_show_stim';
figure_show_stim.name = 'Show stim timings';
figure_show_stim.labels = {'Yes','No'};
figure_show_stim.values = {1,0};
figure_show_stim.val  = {1};
figure_show_stim.help = {'Show timings of stimuli on figures.'}';

figure_rebase_to_zero_at_stim      = cfg_menu;
figure_rebase_to_zero_at_stim.tag  = 'figure_rebase_to_zero_at_stim';
figure_rebase_to_zero_at_stim.name = 'Rebase series to zero at stimuli';
figure_rebase_to_zero_at_stim.labels = {'Yes','No'};
figure_rebase_to_zero_at_stim.values = {1,0};
figure_rebase_to_zero_at_stim.val  = {0};
figure_rebase_to_zero_at_stim.help = {'(Affects only the display of figures: rebase'
    'the series to zero after each stimuli.'}';

show_mse      = cfg_menu;
show_mse.tag  = 'show_mse';
show_mse.name = 'Show MSE';
show_mse.labels = {'Yes','No'};
show_mse.values = {1,0};
show_mse.val  = {1};
show_mse.help = {'Show mean square error on figures.'}';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Executable Branch
glm_roi1      = cfg_exbranch;       % This is the branch that has information about how to run this module
glm_roi1.name = 'GLM on ROI';             % The display name
glm_roi1.tag  = 'glm_roi1'; %Very important: tag is used when calling for execution
glm_roi1.val  = {IOImat ROImat redo1 IOImatCopyChoice session_choice ROI_choice...
     bases volt use_onset_amplitudes hpf_butter lpf_gauss IC ...
     generate_figures save_figures ...
     figure_show_stim figure_rebase_to_zero_at_stim show_mse};    % The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs
glm_roi1.prog = @ioi_GLM_on_ROI_run;  % A function handle that will be called with the harvested job to run the computation
glm_roi1.vout = @ioi_cfg_vout_glm_roi; % A function handle that will be called with the harvested job to determine virtual outputs
glm_roi1.help = {'Run GLMs on ROIs.'};

return

%make IOI.mat available as a dependency
function vout = ioi_cfg_vout_glm_roi(job)
vout = cfg_dep;                     % The dependency object
vout.sname      = 'IOI.mat';       % Displayed dependency name
vout.src_output = substruct('.','IOImat'); %{1}); %,'IOImat');
%substruct('()',{1}); % The output subscript reference. This could be any reference into the output variable created during computation
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
