function glm1 = ioi_GLM_cfg
%_______________________________________________________________________
% Copyright (C) 2011 LIOM Laboratoire d'Imagerie Optique et Mol�culaire
%                    �cole Polytechnique de Montr�al

IOImat         = cfg_files; %Select NIRS.mat for this subject 
IOImat.name    = 'Select IOI.mat'; % The displayed name
IOImat.tag     = 'IOImat';       %file names
IOImat.filter = 'mat';
IOImat.ufilter = '^IOI.mat$';    
IOImat.num     = [1 Inf];     % Number of inputs required 
IOImat.help    = {'Select IOImat dependency if available. '
    'Otherwise, for each subject, select IOI.mat.'}'; % help text displayed

ROImat         = cfg_files; 
ROImat.name    = 'Select ROI.mat'; 
ROImat.tag     = 'ROImat';      
ROImat.filter = 'mat';
ROImat.ufilter = '^ROI.mat$';
ROImat.val     = {''};
ROImat.num     = [1 Inf];    
ROImat.help    = {'Optional: Select ROImat. This allows working on ROI data even'
    'if the paths are not correct in IOI.mat. If not specified, the ROI.mat '
    'specified in IOI.mat will be used.'
    'If several subjects are run, and ROImat is explicitly specified, then'
    'it should be specified for all subjects, in the same order as the list of IOI.mat'}'; 

redo1      = cfg_menu;
redo1.tag  = 'force_redo';
redo1.name = 'Force processing';
redo1.labels = {'False','True'};
redo1.values = {0,1};
redo1.val  = {0};
redo1.help = {'Force redoing this processing even when it has been done already.'
    'Use option below for treatment of previous ROIs.'}';

IOImatOverwrite         = cfg_branch;
IOImatOverwrite.tag     = 'IOImatOverwrite';
IOImatOverwrite.name    = 'Overwrite IOI.mat structure'; 
IOImatOverwrite.help    = {'Will not copy IOI structure.'
            'This will write over the previous NIRS.mat'}';

NewIOIdir         = cfg_entry;
NewIOIdir.name    = 'New directory for IOI.mat';
NewIOIdir.tag     = 'NewIOIdir';       
NewIOIdir.strtype = 's';
NewIOIdir.val{1}    = 'GLM';
NewIOIdir.num     = [1 Inf];     
NewIOIdir.help    = {'Directory for IOI.mat.'}'; 

IOImatCopy         = cfg_branch;
IOImatCopy.tag     = 'IOImatCopy';
IOImatCopy.name    = 'Create new directory and copy IOI structure there'; 
IOImatCopy.val     = {NewIOIdir};
IOImatCopy.help    = {'Create new directory and copy IOI structure there.'}';
        
%Common to most modules: for creating a new directory and copying IOI.mat
IOImatCopyChoice           = cfg_choice;
IOImatCopyChoice.name      = 'Choose IOI copy method';
IOImatCopyChoice.tag       = 'IOImatCopyChoice';
IOImatCopyChoice.values    = {IOImatOverwrite IOImatCopy}; 
IOImatCopyChoice.val       = {IOImatOverwrite}; 
IOImatCopyChoice.help      = {'Choose whether to overwrite the IOI.mat structure'
            'or to create a new directory'
            'and copy the IOI.mat structure there'}'; 
        
%%%%%%%%%%%%%%%%%%%%%%

shrink_x      = cfg_entry;
shrink_x.tag  = 'shrink_x';
shrink_x.name = 'Shrink factor for x dimension';
shrink_x.strtype  = 'i';
shrink_x.num = [1 1];
shrink_x.val  = {2};
shrink_x.help = {'Data reduction factor in x.'};

shrink_y      = cfg_entry;
shrink_y.tag  = 'shrink_y';
shrink_y.name = 'Shrink factor for y dimension';
shrink_y.strtype  = 'i';
shrink_y.num = [1 1];
shrink_y.val = {2};
shrink_y.help = {'Data reduction factor in y.'};

force_shrink_recompute      = cfg_menu;
force_shrink_recompute.tag  = 'force_shrink_recompute';
force_shrink_recompute.name = 'Force shrink recompute';
force_shrink_recompute.labels = {'Yes','No'};
force_shrink_recompute.values = {1,0};
force_shrink_recompute.val  = {0};
force_shrink_recompute.help = {'This option is used when images of one shrunk size are present'
    'But a different size is desired. This forces generating new shrunk images. The old ones'
    'will be kept too, but will not be available from the new IOI.mat. Therefore, one should'
    'usually select the option to place the new IOI.mat in a new folder, so that the old IOI.mat'
    'can still be used to access the first set of shrunk images.'}';


configuration_shrink         = cfg_branch;
configuration_shrink.tag     = 'configuration_shrink';
configuration_shrink.name    = 'Configuration shrinkage';
configuration_shrink.val     = {shrink_x shrink_y force_shrink_recompute};
configuration_shrink.help    = {'Select values.'};

no_shrinkage         = cfg_branch;
no_shrinkage.tag     = 'no_shrinkage';
no_shrinkage.name    = 'No shrinkage';
no_shrinkage.val     = {};
no_shrinkage.help    = {};

shrinkage_choice        = cfg_choice;
shrinkage_choice.name   = 'Choose shrinkage method';
shrinkage_choice.tag    = 'shrinkage_choice';
shrinkage_choice.values = {no_shrinkage,configuration_shrink};
shrinkage_choice.val    = {configuration_shrink};
shrinkage_choice.help   = {'Choose whether to shrink the data. Images will then be stored. And could be reused later.'}';

%%%%%%%%%%%%%%%%%%%%%%%%%%%

vasomotion_on         = cfg_branch;
vasomotion_on.tag     = 'vasomotion_on';
vasomotion_on.name    = 'Include vasomotion regressor';
vasomotion_on.val     = {};
vasomotion_on.help    = {'Include vasomotion regressor.'};

no_vasomotion         = cfg_branch;
no_vasomotion.tag     = 'no_vasomotion';
no_vasomotion.name    = 'No vasomotion regressor';
no_vasomotion.val     = {};
no_vasomotion.help    = {};

vasomotion_choice        = cfg_choice;
vasomotion_choice.name   = 'Choose vasomotion method';
vasomotion_choice.tag    = 'vasomotion_choice';
vasomotion_choice.values = {no_vasomotion,vasomotion_on};
vasomotion_choice.val    = {vasomotion_on};
vasomotion_choice.help   = {'Choose whether to include a vasomotion regressor in the GLM.'}';

%%%%%%%%%%%%%%%%%%%%%%%%%
spatial_LPF_radius         = cfg_entry;
spatial_LPF_radius.name    = 'Spatial LPF radius';
spatial_LPF_radius.tag     = 'spatial_LPF_radius';
spatial_LPF_radius.strtype = 'r';
spatial_LPF_radius.num     = [1 1];
spatial_LPF_radius.val     = {1};
spatial_LPF_radius.help    = {'Enter radius of spatial low pass filter in pixels.'
    'In practice, a radius of 1 gives a weight of 0.4 to the central pixel, of 0.1 to the 4 nearest'
    'and the remainder to the next 8 pixels.'}';

spatial_LPF_On         = cfg_branch;
spatial_LPF_On.tag     = 'spatial_LPF_On';
spatial_LPF_On.name    = 'Spatial LP filter';
spatial_LPF_On.val     = {spatial_LPF_radius};
spatial_LPF_On.help    = {'Spatial low-pass filter.'};

spatial_LPF_Off         = cfg_branch;
spatial_LPF_Off.tag     = 'spatial_LPF_Off';
spatial_LPF_Off.name    = 'Spatial filter off';
spatial_LPF_Off.val     = {};
spatial_LPF_Off.help    = {'Spatial low pass filter turned off.'};

spatial_LPF      = cfg_choice;
spatial_LPF.tag  = 'spatial_LPF';
spatial_LPF.name = 'Spatial Low Pass Filter';
spatial_LPF.values = {spatial_LPF_On spatial_LPF_Off};
spatial_LPF.val = {spatial_LPF_Off};
spatial_LPF.help = {'Choose whether to include a spatial Low Pass Filter'
    'on the data prior to running the GLM.'}';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%

all_ROIs         = cfg_branch;
all_ROIs.tag     = 'all_ROIs';
all_ROIs.name    = 'All ROIs';
all_ROIs.val     = {};
all_ROIs.help    = {'All ROIs will be processed'};

selected_ROIs      = cfg_entry;
selected_ROIs.tag  = 'selected_ROIs';
selected_ROIs.name = 'Enter list of ROIs';
selected_ROIs.strtype  = 'r';
selected_ROIs.num = [1 Inf];
selected_ROIs.val{1} = 1;
selected_ROIs.help = {'Enter list of ROIs to process.'};

select_ROIs         = cfg_branch;
select_ROIs.tag     = 'select_ROIs';
select_ROIs.name    = 'Select ROIs';
select_ROIs.val     = {selected_ROIs};
select_ROIs.help    = {'Choose some ROIs to be processed'};

ROI_choice        = cfg_choice;
ROI_choice.name   = 'Choose ROI selection method';
ROI_choice.tag    = 'ROI_choice';
ROI_choice.values = {all_ROIs,select_ROIs};
ROI_choice.val    = {all_ROIs};
ROI_choice.help   = {'Choose ROI selection method'}';

image_mode         = cfg_branch;
image_mode.tag     = 'image_mode';
image_mode.name    = 'Image mode';
image_mode.val     = {spatial_LPF shrinkage_choice};
image_mode.help    = {'Data will be images'};

ROI_mode         = cfg_branch;
ROI_mode.tag     = 'ROI_mode';
ROI_mode.name    = 'ROI mode';
ROI_mode.val     = {ROImat ROI_choice figure_show_stim figure_rebase_to_zero_at_stim show_mse};
ROI_mode.help    = {'Data will be ROIs'};

data_selection_choice        = cfg_choice;
data_selection_choice.name   = 'Choose ROI selection method';
data_selection_choice.tag    = 'data_selection_choice';
data_selection_choice.values = {image_mode,ROI_mode};
data_selection_choice.val    = {image_mode};
data_selection_choice.help   = {'Choose data selection method'}';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

all_sessions         = cfg_branch;
all_sessions.tag     = 'all_sessions';
all_sessions.name    = 'All sessions';
all_sessions.val     = {};
all_sessions.help    = {'All good enough sessions will be processed'};

selected_sessions      = cfg_entry;
selected_sessions.tag  = 'selected_sessions';
selected_sessions.name = 'Enter list of sessions';
selected_sessions.strtype  = 'r';
selected_sessions.num = [1 Inf];
selected_sessions.val{1} = 1;
selected_sessions.help = {'Enter list of sessions to process.'};

select_sessions         = cfg_branch;
select_sessions.tag     = 'select_sessions';
select_sessions.name    = 'Select sessions';
select_sessions.val     = {selected_sessions};
select_sessions.help    = {'Choose some sessions to be processed'};

session_choice        = cfg_choice;
session_choice.name   = 'Choose session selection method';
session_choice.tag    = 'session_choice';
session_choice.values = {all_sessions,select_sessions};
session_choice.val    = {all_sessions};
session_choice.help   = {'Choose session selection method'}';

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

hpf_butter_freq         = cfg_entry; 
hpf_butter_freq.name    = 'Cutoff frequency for HPF';
hpf_butter_freq.tag     = 'hpf_butter_freq';       
hpf_butter_freq.strtype = 'r';
hpf_butter_freq.num     = [1 1];     
hpf_butter_freq.val     = {0.01};
hpf_butter_freq.help    = {'Enter cutoff frequency in Hz for Butterworth HPF.'};

hpf_butter_order         = cfg_entry; 
hpf_butter_order.name    = 'Order of Butterworth HPF';
hpf_butter_order.tag     = 'hpf_butter_order';       
hpf_butter_order.strtype = 'r';
hpf_butter_order.num     = [1 1];     
hpf_butter_order.val     = {3};
hpf_butter_order.help    = {'Enter order of Butterworth HPF (preferred value = 3).'};

hpf_butter_On         = cfg_branch;
hpf_butter_On.tag     = 'hpf_butter_On';
hpf_butter_On.name    = 'Butterworth HP filter';
hpf_butter_On.val     = {hpf_butter_freq hpf_butter_order}; 
hpf_butter_On.help    = {'Butterworth high-pass filter.'};

hpf_butter_Off         = cfg_branch;
hpf_butter_Off.tag     = 'hpf_butter_Off';
hpf_butter_Off.name    = 'HP filter off';
hpf_butter_Off.val     = {}; 
hpf_butter_Off.help    = {'High pass filter turned off.'};

hpf_butter      = cfg_choice;
hpf_butter.tag  = 'hpf_butter';
hpf_butter.name = 'Butterworth High Pass Filter';
hpf_butter.values = {hpf_butter_On hpf_butter_Off};
hpf_butter.val = {hpf_butter_On};
hpf_butter.help = {'Choose whether to include a Butterworth High Pass Filter.'
        'Parameters are: order (e.g. 3) and frequency (e.g. 0.01 Hz)'}';
 
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

include_flow      = cfg_menu;
include_flow.tag  = 'include_flow';
include_flow.name = 'Include flow';
include_flow.labels = {'Yes','No'};
include_flow.values = {1,0};
include_flow.val  = {0};
include_flow.help = {'Include flow.'}';

include_HbT      = cfg_menu;
include_HbT.tag  = 'include_HbT';
include_HbT.name = 'Include HbT';
include_HbT.labels = {'Yes','No'};
include_HbT.values = {1,0};
include_HbT.val  = {1};
include_HbT.help = {'Include HbT.'}';

include_OD      = cfg_menu;
include_OD.tag  = 'include_OD';
include_OD.name = 'Include optical intensity';
include_OD.labels = {'Yes','No'};
include_OD.values = {1,0};
include_OD.val  = {0};
include_OD.help = {'If the optical intensity images (Green, Red, Yellow) have not been deleted'
    'previously, choose whether to generate movies for these colors.'}';

include_HbO      = cfg_menu;
include_HbO.tag  = 'include_HbO';
include_HbO.name = 'Include HbO';
include_HbO.labels = {'Yes','No'};
include_HbO.values = {1,0};
include_HbO.val  = {1};
include_HbO.help = {'Include HbO.'}';

include_HbR      = cfg_menu;
include_HbR.tag  = 'include_HbR';
include_HbR.name = 'Include HbR';
include_HbR.labels = {'Yes','No'};
include_HbR.values = {1,0};
include_HbR.val  = {1};
include_HbR.help = {'Include HbR.'}';

generate_figures      = cfg_menu;
generate_figures.tag  = 'generate_figures';
generate_figures.name = 'Show figures';
generate_figures.labels = {'Yes','No'};
generate_figures.values = {1,0};
generate_figures.val  = {0};
generate_figures.help = {'Show figures. When selecting this option, the figures will stay opened after the code has completed.'}';

save_figures      = cfg_menu;
save_figures.tag  = 'save_figures';
save_figures.name = 'Save figures';
save_figures.labels = {'Yes','No'};
save_figures.values = {1,0};
save_figures.val  = {1};
save_figures.help = {'Save figures.'}';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Executable Branch
glm1      = cfg_exbranch;       % This is the branch that has information about how to run this module
glm1.name = 'GLM (images or ROIs)';             % The display name
glm1.tag  = 'glm1'; %Very important: tag is used when calling for execution
glm1.val  = {IOImat data_selection_choice redo1 IOImatCopyChoice ...
     session_choice ...
     bases volt use_onset_amplitudes hpf_butter lpf_gauss vasomotion_choice include_flow ...
     include_HbT include_HbO include_HbR include_OD ...
     generate_figures save_figures};    % The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs
glm1.prog = @ioi_GLM_run;  % A function handle that will be called with the harvested job to run the computation
glm1.vout = @ioi_cfg_vout_glm; % A function handle that will be called with the harvested job to determine virtual outputs
glm1.help = {'Run GLMs.'};

return

%make IOI.mat available as a dependency
function vout = ioi_cfg_vout_glm(job)
vout = cfg_dep;                     % The dependency object
vout.sname      = 'IOI.mat';       % Displayed dependency name
vout.src_output = substruct('.','IOImat'); %{1}); %,'IOImat');
%substruct('()',{1}); % The output subscript reference. This could be any reference into the output variable created during computation
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});