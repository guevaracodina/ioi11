function hdm1 = ioi_HDM_cfg
%_______________________________________________________________________
% Copyright (C) 2011 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal

IOImat         = cfg_files; %Select NIRS.mat for this subject 
IOImat.name    = 'Select IOI.mat'; % The displayed name
IOImat.tag     = 'IOImat';       %file names
IOImat.filter = 'mat';
IOImat.ufilter = '^IOI.mat$';    
IOImat.num     = [1 Inf];     % Number of inputs required 
IOImat.help    = {'Select IOImat dependency if available. '
    'Otherwise, for each subject, select IOI.mat.'}'; % help text displayed

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
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hemodynamic Model choice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% BF         = cfg_branch;
% BF.tag     = 'BF';
% BF.name    = 'Buxton-Friston';
% BF.val     = {};
% BF.help    = {'Buxton-Friston'};
% 
% ZM         = cfg_branch;
% ZM.tag     = 'ZM';
% ZM.name    = 'Zheng-Mayhew';
% ZM.val     = {};
% ZM.help    = {'Zheng-Mayhew'};
% 
% BH         = cfg_branch;
% BH.tag     = 'BH';
% BH.name    = 'Boas-Huppert';
% BH.val     = {};
% BH.help    = {'Boas-Huppert'};
% 
% Model_Choice        = cfg_choice;
% Model_Choice.name   = 'Choose hemodynamic model';
% Model_Choice.tag    = 'Model_Choice';
% Model_Choice.values = {BF,ZM,BH};
% Model_Choice.val    = {BF};
% Model_Choice.help   = {'Choose hemodynamic model'}';

Model_Choice        = cfg_menu;
Model_Choice.name   = 'Choose hemodynamic model';
Model_Choice.tag    = 'Model_Choice';
Model_Choice.labels = {'Buxton-Friston','Zheng-Mayhew','Boas-Huppert'};
Model_Choice.values = {0,1,2};
Model_Choice.val    = {0};
Model_Choice.help   = {'Choose hemodynamic model'}';


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
% Electrophysiology
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sf      = cfg_entry;
sf.tag  = 'sf';
sf.name = 'Enter electrophysiology sampling frequency';
sf.strtype  = 'r';
sf.num = [1 1];
sf.val{1} = 10000;
sf.help = {'Enter electrophysiology sampling frequency.'}';

nSD      = cfg_entry;
nSD.tag  = 'nSD';
nSD.name = 'Enter number of standard deviations';
nSD.strtype  = 'r';
nSD.num = [1 1];
nSD.val{1} = 2; %1.8; %3;
nSD.help = {'Enter number of standard deviations above mean of electrophysiology'
    'signal to use for detection of peaks.'}';

mbSD      = cfg_entry;
mbSD.tag  = 'mbSD';
mbSD.name = 'Enter minimum standard deviation';
mbSD.strtype  = 'r';
mbSD.num = [1 1];
mbSD.val{1} = 0.02;
mbSD.help = {'Enter minimum standard deviation to be used for detection.'
    'The detection threshold will be the mean signal plus the specified'
    'multiple of standard deviations, times the maximum of the calculated'
    'standard deviation of the electrophysiological signal after filtering'
    'and the specified minimum standard deviation.'}';

dP      = cfg_entry;
dP.tag  = 'dP';
dP.name = 'Enter minimal peak distance';
dP.strtype  = 'r';
dP.num = [1 1];
dP.val{1} = 25;
dP.help = {'Enter minimal distance allowed between peaks in milliseconds'}';

%HPF
electro_hpf_butter_freq         = cfg_entry; 
electro_hpf_butter_freq.name    = 'Cutoff frequency for HPF';
electro_hpf_butter_freq.tag     = 'electro_hpf_butter_freq';       
electro_hpf_butter_freq.strtype = 'r';
electro_hpf_butter_freq.num     = [1 1];     
electro_hpf_butter_freq.val     = {0.2};
electro_hpf_butter_freq.help    = {'Enter cutoff frequency in Hz for Butterworth HPF.'};

electro_hpf_butter_order         = cfg_entry; 
electro_hpf_butter_order.name    = 'Order of Butterworth HPF';
electro_hpf_butter_order.tag     = 'electro_hpf_butter_order';       
electro_hpf_butter_order.strtype = 'r';
electro_hpf_butter_order.num     = [1 1];     
electro_hpf_butter_order.val     = {3};
electro_hpf_butter_order.help    = {'Enter order of Butterworth HPF (preferred value = 3).'};

electro_hpf_butter_On         = cfg_branch;
electro_hpf_butter_On.tag     = 'electro_hpf_butter_On';
electro_hpf_butter_On.name    = 'Butterworth HP filter';
electro_hpf_butter_On.val     = {electro_hpf_butter_freq electro_hpf_butter_order}; 
electro_hpf_butter_On.help    = {'Butterworth high-pass filter.'};

electro_hpf_butter_Off         = cfg_branch;
electro_hpf_butter_Off.tag     = 'electro_hpf_butter_Off';
electro_hpf_butter_Off.name    = 'HP filter off';
electro_hpf_butter_Off.val     = {}; 
electro_hpf_butter_Off.help    = {'High pass filter turned off.'};

electro_hpf_butter      = cfg_choice;
electro_hpf_butter.tag  = 'electro_hpf_butter';
electro_hpf_butter.name = 'Butterworth High Pass Filter';
electro_hpf_butter.values = {electro_hpf_butter_On electro_hpf_butter_Off};
electro_hpf_butter.val = {electro_hpf_butter_On};
electro_hpf_butter.help = {'Choose whether to include a Butterworth High Pass Filter.'
        'Parameters are: order (e.g. 3) and frequency (e.g. 0.01 Hz)'}';

%LPF
electro_lpf_butter_freq         = cfg_entry; 
electro_lpf_butter_freq.name    = 'Cutoff frequency for LPF';
electro_lpf_butter_freq.tag     = 'electro_lpf_butter_freq';       
electro_lpf_butter_freq.strtype = 'r';
electro_lpf_butter_freq.num     = [1 1];     
electro_lpf_butter_freq.val     = {130};
electro_lpf_butter_freq.help    = {'Enter cutoff frequency in Hz for Butterworth LPF.'};

electro_lpf_butter_order         = cfg_entry; 
electro_lpf_butter_order.name    = 'Order of Butterworth LPF';
electro_lpf_butter_order.tag     = 'electro_lpf_butter_order';       
electro_lpf_butter_order.strtype = 'r';
electro_lpf_butter_order.num     = [1 1];     
electro_lpf_butter_order.val     = {3};
electro_lpf_butter_order.help    = {'Enter order of Butterworth LPF (preferred value = 3).'};

electro_lpf_butter_On         = cfg_branch;
electro_lpf_butter_On.tag     = 'electro_lpf_butter_On';
electro_lpf_butter_On.name    = 'Butterworth LP filter';
electro_lpf_butter_On.val     = {electro_lpf_butter_freq electro_lpf_butter_order}; 
electro_lpf_butter_On.help    = {'Butterworth low-pass filter.'};

electro_lpf_butter_Off         = cfg_branch;
electro_lpf_butter_Off.tag     = 'electro_lpf_butter_Off';
electro_lpf_butter_Off.name    = 'LP filter off';
electro_lpf_butter_Off.val     = {}; 
electro_lpf_butter_Off.help    = {'Low pass filter turned off.'};

electro_lpf_butter      = cfg_choice;
electro_lpf_butter.tag  = 'electro_lpf_butter';
electro_lpf_butter.name = 'Butterworth Low Pass Filter';
electro_lpf_butter.values = {electro_lpf_butter_On electro_lpf_butter_Off};
electro_lpf_butter.val = {electro_lpf_butter_On};
electro_lpf_butter.help = {'Choose whether to include a Butterworth Low Pass Filter.'
        'Parameters are: order (e.g. 3) and frequency (e.g. 130 Hz)'}';

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

% ---------------------------------------------------------------------
% hrf Canonical HRF
% ---------------------------------------------------------------------
rat         = cfg_branch;
rat.tag     = 'rat';
rat.name    = 'HRF for rats';
rat.val     = {derivs };
rat.help    = {'Rat Hemodynamic Response Function.'};

mouse         = cfg_branch;
mouse.tag     = 'mouse';
mouse.name    = 'HRF for mice';
mouse.val     = {derivs };
mouse.help    = {'Mouse Hemodynamic Response Function.'};

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
bases.values  = {hrf rat mouse fourier fourier_han gamma fir };
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
hpf_butter_freq.val     = {0.05};
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
fwhm1.val = {1};
fwhm1.strtype = 'r';  
fwhm1.num     = [1 1]; 
fwhm1.help    = {'FWHM in seconds.'}; 

lpf_gauss         = cfg_branch;
lpf_gauss.tag     = 'lpf_gauss';
lpf_gauss.name    = 'Gaussian Filter';
lpf_gauss.val     = {fwhm1}; 
lpf_gauss.help    = {'Specify properties of Gaussian filter'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

use_epilepsy_convention      = cfg_menu;
use_epilepsy_convention.tag  = 'use_epilepsy_convention';
use_epilepsy_convention.name = 'Use epilepsy convention for plots';
use_epilepsy_convention.labels = {'False','True'};
use_epilepsy_convention.values = {0,1};
use_epilepsy_convention.val  = {1};
use_epilepsy_convention.help = {'Use epilepsy/neurology convention for plots '
    'of electrophysiology: vertical axis is inverted.'}';

write_pictures      = cfg_menu;
write_pictures.tag  = 'write_pictures';
write_pictures.name = 'Write plots of electrophysiology';
write_pictures.labels = {'False','True'};
write_pictures.values = {0,1};
write_pictures.val  = {1};
write_pictures.help = {'Generate plots of electrophysiology.'
    'Currently only applies to choice Stimulations from electrophysiology.'}';

electro_stims         = cfg_branch;
electro_stims.tag     = 'electro_stims';
electro_stims.name    = 'Stimulations from electrophysiology';
electro_stims.val     = {sf nSD mbSD dP electro_hpf_butter electro_lpf_butter ...
        write_pictures use_epilepsy_convention};
electro_stims.help    = {'Electrophysiology information'
    'Stimulations are assumed to last one data point.'
    'Information stored in IOI.Sess (not to be confused with protocol info in IOI.sess_res).'}';

default_stims         = cfg_branch;
default_stims.tag     = 'default_stims';
default_stims.name    = 'Default stimulations';
default_stims.val     = {};
default_stims.help    = {'The names, onsets, durations fields in IOI.sess_res will be used'};

stim_choice        = cfg_choice;
stim_choice.name   = 'Choose stimulation selection method';
stim_choice.tag    = 'stim_choice';
stim_choice.values = {default_stims,electro_stims};
stim_choice.val    = {electro_stims};
stim_choice.help   = {'Choose stimulation selection method'}';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Executable Branch
hdm1      = cfg_exbranch;       % This is the branch that has information about how to run this module
hdm1.name = 'HDM on ROI';             % The display name
hdm1.tag  = 'hdm1'; %Very important: tag is used when calling for execution
hdm1.val  = {IOImat redo1 IOImatCopyChoice Model_Choice session_choice ROI_choice...
     bases volt hpf_butter lpf_gauss stim_choice};    % The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs
hdm1.prog = @ioi_HDM_run;  % A function handle that will be called with the harvested job to run the computation
hdm1.vout = @ioi_cfg_vout_HDM; % A function handle that will be called with the harvested job to determine virtual outputs
hdm1.help = {'Run hemodynamic modelling (HDM) on average ROI time courses.'};

return

%make IOI.mat available as a dependency
function vout = ioi_cfg_vout_HDM(job)
vout = cfg_dep;                     % The dependency object
vout.sname      = 'IOI.mat';       % Displayed dependency name
vout.src_output = substruct('.','IOImat'); %{1}); %,'IOImat');
%substruct('()',{1}); % The output subscript reference. This could be any reference into the output variable created during computation
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
