function ROC1 = ioi_ROC_cfg
%_______________________________________________________________________
% Copyright (C) 2011 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal

IOImat         = cfg_files; 
IOImat.name    = 'Select IOI.mat'; 
IOImat.tag     = 'IOImat';       
IOImat.filter = 'mat';
IOImat.ufilter = '^IOI.mat$';    
IOImat.num     = [1 Inf];     
IOImat.help    = {'Select IOImat dependency if available. '
    'Otherwise, for each subject, select IOI.mat.'}';

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
NewIOIdir.val{1}  = 'SCKS';
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

electro_stims         = cfg_branch;
electro_stims.tag     = 'electro_stims';
electro_stims.name    = 'Stimulations from electrophysiology';
electro_stims.val     = {};
electro_stims.help    = {'Electrophysiology information'
    'Taken from IOI.Sess().U, which is generated by GLM module.'
    'Thus to use this option, GLM module must be run first.'}';

manual_stims         = cfg_branch;
manual_stims.tag     = 'manual_stims';
manual_stims.name    = 'Stimulations manually defined';
manual_stims.val     = {};
manual_stims.help    = {'Stimulations are defined by user inputs'
    'This defines IOI.Sess().U, which is generated by user inputs.'}';

default_stims         = cfg_branch;
default_stims.tag     = 'default_stims';
default_stims.name    = 'Default stimulations';
default_stims.val     = {};
default_stims.help    = {'The names, onsets, durations fields in IOI.sess_res will be used'};

stim_choice        = cfg_choice;
stim_choice.name   = 'Choose stimulation selection method';
stim_choice.tag    = 'stim_choice';
stim_choice.values = {default_stims,electro_stims,manual_stims};
stim_choice.val    = {default_stims};
stim_choice.help   = {'Choose stimulation selection method'
       'NOTE: this electrophysiology option only applies when electrophysiological events '
    'have been detected in an earlier module. Currently only the GLM module performs this detection.'}';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

includeHbR      = cfg_menu;
includeHbR.tag  = 'includeHbR';
includeHbR.name = 'Include HbR';
includeHbR.labels = {'Yes','No'};
includeHbR.values = {1,0};
includeHbR.val  = {1};
includeHbR.help = {'Include modality HbR.'}';

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

PhysioModel_Choice        = cfg_menu;
PhysioModel_Choice.name   = 'Choose hemodynamic model';
PhysioModel_Choice.tag    = 'PhysioModel_Choice';
PhysioModel_Choice.labels = {'Buxton-Friston','Zheng-Mayhew','Boas-Huppert'};
PhysioModel_Choice.values = {0,1,2};
PhysioModel_Choice.val    = {0};
PhysioModel_Choice.help   = {'Choose hemodynamic model'}';

SCKSnoise      = cfg_menu;
SCKSnoise.tag  = 'SCKSnoise';
SCKSnoise.name = 'Promote parameters and noise to states';
SCKSnoise.labels = {'Yes','No'};
SCKSnoise.values = {1,0};
SCKSnoise.val  = {0};
SCKSnoise.help = {'Promote parameters and noise to time-dependent states.'}';

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

State_annealing         = cfg_entry; 
State_annealing.name    = 'Annealing factor on states';
State_annealing.tag     = 'State_annealing';       
State_annealing.strtype = 'r';
State_annealing.num     = [1 1];     
State_annealing.val     = {0.9995};
State_annealing.help    = {'Enter annealing factor on states. This is a number between 0 and 1.'
    'A value closer to 1 corresponds to more annealing. A value of 0.5 corresponds to very little annealing.'}';

Parameter_annealing         = cfg_entry; 
Parameter_annealing.name    = 'Annealing factor on parameters';
Parameter_annealing.tag     = 'Parameter_annealing';       
Parameter_annealing.strtype = 'r';
Parameter_annealing.num     = [1 1];     
Parameter_annealing.val     = {0.9995};
Parameter_annealing.help    = {'Enter annealing factor on parameters.' 
    'This is a number between 0 and 1.'
    'A value closer to 1 corresponds to more annealing. A value of 0.5 corresponds to very little annealing.'}';

SCKSparams         = cfg_branch;
SCKSparams.tag     = 'SCKSparams';
SCKSparams.name    = 'SCKS parameters';
SCKSparams.val     = {SCKSnoise State_annealing Parameter_annealing}; 
SCKSparams.help    = {'User-controlled SCKS parameters.'};

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
save_figures.val  = {0};
save_figures.help = {'Save figures.'}';

% Executable Branch
ROC1      = cfg_exbranch;       % This is the branch that has information about how to run this module
ROC1.name = 'ROC curves';             % The display name
ROC1.tag  = 'ROC1'; %Very important: tag is used when calling for execution
ROC1.val  = {IOImat ROImat redo1 IOImatCopyChoice session_choice ROI_choice ...
    PhysioModel_Choice includeHbR includeHbT includeFlow hpf_butter ...
    SCKSparams generate_figures save_figures};    % The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs
ROC1.prog = @ioi_ROC_run;  % A function handle that will be called with the harvested job to run the computation
ROC1.vout = @ioi_cfg_vout_ROC; % A function handle that will be called with the harvested job to determine virtual outputs
ROC1.help = {'Create regions of interest.'};

return

%make IOI.mat available as a dependency
function vout = ioi_cfg_vout_ROC(job)
vout = cfg_dep;                     % The dependency object
vout.sname      = 'IOI.mat';       % Displayed dependency name
vout.src_output = substruct('.','IOImat'); %{1}); %,'IOImat');
%substruct('()',{1}); % The output subscript reference. This could be any reference into the output variable created during computation
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
