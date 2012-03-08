function [IOImat ROImat redo1 IOImatCopyChoice includeHbR includeHbO ...
    includeHbT includeFlow PhysioModel_Choice ROI_choice session_choice baseline_choice ...
    lpf_choice hpf_butter] = ioi_common_fields_SCKS_HDM

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
NewIOIdir.val{1}  = 'HDM';
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

includeHbR      = cfg_menu;
includeHbR.tag  = 'includeHbR';
includeHbR.name = 'Include HbR';
includeHbR.labels = {'Yes','No'};
includeHbR.values = {1,0};
includeHbR.val  = {1};
includeHbR.help = {'Include modality HbR.'}';

includeHbO      = cfg_menu;
includeHbO.tag  = 'includeHbO';
includeHbO.name = 'Include HbO';
includeHbO.labels = {'Yes','No'};
includeHbO.values = {1,0};
includeHbO.val  = {1};
includeHbO.help = {'Include modality HbO.'}';

includeHbT      = cfg_menu;
includeHbT.tag  = 'includeHbT';
includeHbT.name = 'Include HbT';
includeHbT.labels = {'Yes','No'};
includeHbT.values = {1,0};
includeHbT.val  = {1};
includeHbT.help = {'Include modality HbT.'}';

includeFlow      = cfg_menu;
includeFlow.tag  = 'includeFlow';
includeFlow.name = 'Include Flow';
includeFlow.labels = {'Yes','No'};
includeFlow.values = {1,0};
includeFlow.val  = {1};
includeFlow.help = {'Include modality Flow.'}';

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

PhysioModel_Choice        = cfg_menu;
PhysioModel_Choice.name   = 'Choose hemodynamic model';
PhysioModel_Choice.tag    = 'PhysioModel_Choice';
PhysioModel_Choice.labels = {'Buxton-Friston','Zheng-Mayhew','Boas-Huppert'};
PhysioModel_Choice.values = {0,1,2};
PhysioModel_Choice.val    = {0};
PhysioModel_Choice.help   = {'Choose hemodynamic model'}';

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

lpf_gauss_On         = cfg_branch;
lpf_gauss_On.tag     = 'lpf_gauss_On';
lpf_gauss_On.name    = 'Gaussian LP filter';
lpf_gauss_On.val     = {fwhm1}; 
lpf_gauss_On.help    = {'Gaussian low-pass filter '}';
   % '(applied forward then backward so that it does not create a time shift)'}';

lpf_Off         = cfg_branch;
lpf_Off.tag     = 'lpf_Off';
lpf_Off.name    = 'LP filter off';
lpf_Off.val     = {}; 
lpf_Off.help    = {'Low pass filter turned off.'};

lpf_choice      = cfg_choice;
lpf_choice.tag  = 'lpf_choice';
lpf_choice.name = 'Choose Low Pass Filter';
lpf_choice.values = {lpf_gauss_On lpf_Off};
lpf_choice.val = {lpf_gauss_On};
lpf_choice.help = {'Choose whether to include a Low Pass Filter.'
        'Parameters'}';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
no_baseline_correction         = cfg_branch;
no_baseline_correction.tag     = 'no_baseline_correction';
no_baseline_correction.name    = 'No baseline correction';
no_baseline_correction.val     = {};
no_baseline_correction.help    = {'No correction will be made to baseline'};

baseline_percentile_HbR      = cfg_entry;
baseline_percentile_HbR.tag  = 'baseline_percentile_HbR';
baseline_percentile_HbR.name = 'Choose data percentile to set baseline for HbR';
baseline_percentile_HbR.strtype  = 'r';
baseline_percentile_HbR.num = [Inf Inf];
baseline_percentile_HbR.val{1} = 90;
baseline_percentile_HbR.help = {'Enter percentile of data (after filtering) to set baseline at, for HbR.'
    'Enter either a single number,'
    'to be applied to each selected session and ROI, or '
    'enter an array of offsets, by selected session and by selected ROI.'}';

baseline_percentile_HbT      = cfg_entry;
baseline_percentile_HbT.tag  = 'baseline_percentile_HbT';
baseline_percentile_HbT.name = 'Choose data percentile to set baseline for HbT';
baseline_percentile_HbT.strtype  = 'r';
baseline_percentile_HbT.num = [Inf Inf];
baseline_percentile_HbT.val{1} = 10;
baseline_percentile_HbT.help = {'Enter percentile of data (after filtering) to set baseline at, for HbT.'
    'Enter either a single number,'
    'to be applied to each selected session and ROI, or '
    'enter an array of offsets, by selected session and by selected ROI.'}';

baseline_percentile_flow      = cfg_entry;
baseline_percentile_flow.tag  = 'baseline_percentile_flow';
baseline_percentile_flow.name = 'Choose data percentile to set baseline for flow';
baseline_percentile_flow.strtype  = 'r';
baseline_percentile_flow.num = [Inf Inf];
baseline_percentile_flow.val{1} = 10;
baseline_percentile_flow.help = {'Enter percentile of data (after filtering) to set baseline at, for flow.'
    'Enter either a single number,'
    'to be applied to each selected session and ROI, or '
    'enter an array of offsets, by selected session and by selected ROI.'}';

baseline_percentile_choice         = cfg_branch;
baseline_percentile_choice.tag     = 'baseline_percentile_choice';
baseline_percentile_choice.name    = 'Baseline percentile choice';
baseline_percentile_choice.val     = {baseline_percentile_HbR baseline_percentile_HbT baseline_percentile_flow};
baseline_percentile_choice.help    = {'Set baseline to a chosen percentile'}';

baseline_offset_HbR      = cfg_entry;
baseline_offset_HbR.tag  = 'baseline_offset_HbR';
baseline_offset_HbR.name = 'Baseline offset for HbR';
baseline_offset_HbR.strtype  = 'r';
baseline_offset_HbR.num = [Inf Inf];
baseline_offset_HbR.val{1} = 0;
baseline_offset_HbR.help = {'Enter baseline offset for HbR. Enter either a single number,'
    'to be applied to each selected session and ROI, or  '
    'enter an array of offsets, by selected session and by selected ROI.'}';

baseline_offset_HbT      = cfg_entry;
baseline_offset_HbT.tag  = 'baseline_offset_HbT';
baseline_offset_HbT.name = 'Baseline offset for HbT';
baseline_offset_HbT.strtype  = 'r';
baseline_offset_HbT.num = [Inf Inf];
baseline_offset_HbT.val{1} = 0;
baseline_offset_HbT.help = {'Enter baseline offset for HbT. Enter either a single number,'
    'to be applied to each selected session and ROI, or  '
    'enter an array of offsets, by selected session and by selected ROI.'}';

baseline_offset_flow      = cfg_entry;
baseline_offset_flow.tag  = 'baseline_offset_flow';
baseline_offset_flow.name = 'Baseline offset for Flow';
baseline_offset_flow.strtype  = 'r';
baseline_offset_flow.num = [Inf Inf];
baseline_offset_flow.val{1} = 0;
baseline_offset_flow.help = {'Enter baseline offset for flow. Enter either a single number,'
    'to be applied to each selected session and ROI, or  '
    'enter an array of offsets, by selected session and by selected ROI.'}';

baseline_offset_choice         = cfg_branch;
baseline_offset_choice.tag     = 'baseline_offset_choice';
baseline_offset_choice.name    = 'Baseline offset choice';
baseline_offset_choice.val     = {baseline_offset_HbR baseline_offset_HbT baseline_offset_flow};
baseline_offset_choice.help    = {'Offset baseline by a chosen amount'};

baseline_choice        = cfg_choice;
baseline_choice.name   = 'Choose baseline selection method';
baseline_choice.tag    = 'baseline_choice';
baseline_choice.values = {no_baseline_correction,baseline_percentile_choice,baseline_offset_choice};
baseline_choice.val    = {baseline_percentile_choice};
baseline_choice.help   = {'Choose baseline selection method'}';
