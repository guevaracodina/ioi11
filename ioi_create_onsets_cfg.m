function create_onsets1 = ioi_create_onsets_cfg
%_______________________________________________________________________
% Copyright (C) 2011 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal

IOImat = ioi_dfg_IOImat(1);
redo1 = ioi_dfg_redo(0);

elDir         = cfg_files;
elDir.tag     = 'elDir';
elDir.name    = 'Subject directory for electrophysiology';
elDir.filter = 'dir'; 
elDir.num     = [1 Inf];
elDir.val     = {''};
elDir.help    = {'Optional: Select directory where processed electrophysiology files are located.'
    'This allows working on electrophysiology data even'
    'if the paths are not correct in IOI.mat. If not specified, the el '
    'specified in IOI.mat will be used.'
    'If several subjects are run, and elDir is explicitly specified, then'
    'it should be specified for all subjects, in the same order as the list of IOI.mat'}';

IOImatCopyChoice = ioi_dfg_IOImatCopyChoice('Onsets');        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
session_choice = ioi_dfg_session_choice;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Electrophysiology detection options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

electrophysiology_onset_name         = cfg_entry; 
electrophysiology_onset_name.name    = 'Electrophysiology onset name';
electrophysiology_onset_name.tag     = 'electrophysiology_onset_name';       
electrophysiology_onset_name.strtype = 's';
electrophysiology_onset_name.num     = [1 Inf];     
electrophysiology_onset_name.val{1}  = 'Spk';
electrophysiology_onset_name.help    = {'Specify the onset name.'};

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

tb         = cfg_entry; 
tb.name    = 'Time window before';
tb.tag     = 'tb';       
tb.strtype = 'r';
tb.num     = [1 1];     
tb.val{1}  = 0.2;
tb.help    = {'Specify the time window before onset to use for display (in seconds).'}';

ta         = cfg_entry; 
ta.name    = 'Time window after';
ta.tag     = 'ta';       
ta.strtype = 'r';
ta.num     = [1 1];     
ta.val{1}  = 0.5;
ta.help    = {'Specify the time window after onset to use for display (in seconds).'}';

%HPF
electro_hpf_butter = ioi_dfg_hpf_butter(1,0.2,3);

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
write_pictures.name = 'Save figures of electrophysiology';
write_pictures.labels = {'False','True'};
write_pictures.values = {0,1};
write_pictures.val  = {1};
write_pictures.help = {'Generate plots of electrophysiology.'}';

%%%%%%%%%%%%%%%%%%%%

electro_stims         = cfg_branch;
electro_stims.tag     = 'electro_stims';
electro_stims.name    = 'Onsets from electrophysiology';
electro_stims.val     = {electrophysiology_onset_name ...
    sf nSD mbSD dP tb ta electro_hpf_butter electro_lpf_butter ...
        write_pictures use_epilepsy_convention};
electro_stims.help    = {    'Electrophysiology information'
    'Stimulations are assumed to last one data point.'
    'Information stored in IOI.Sess (not to be confused with protocol info in IOI.sess_res).'}';

onset_name         = cfg_entry; 
onset_name.name    = 'Onset name';
onset_name.tag     = 'onset_name';       
onset_name.strtype = 's';
onset_name.num     = [1 Inf];     
onset_name.val{1}  = 'Stim1';
onset_name.help    = {'Specify the onset name.'};

onset_times         = cfg_entry; 
onset_times.name    = 'Onset times';
onset_times.tag     = 'onset_times';       
onset_times.strtype = 'r';
onset_times.num     = [0 Inf];     
onset_times.val{1}  = '';
onset_times.help    = {'Specify the onset times in seconds.'};

onset_durations         = cfg_entry; 
onset_durations.name    = 'Onset durations';
onset_durations.tag     = 'onset_durations';       
onset_durations.strtype = 'r';
onset_durations.num     = [1 Inf];     
onset_durations.val{1}  = 0;
onset_durations.help    = {'Specify the onset durations in seconds.'
    'Specify either one number, applicable to all onsets,'
    'or a list of durations for each onset'}';

onset_amplitude         = cfg_entry; 
onset_amplitude.name    = 'Onset amplitude';
onset_amplitude.tag     = 'onset_amplitude';       
onset_amplitude.strtype = 'r';
onset_amplitude.num     = [1 Inf];     
onset_amplitude.val{1}  = 1;
onset_amplitude.help    = {'Specify the onset amplitudes in seconds.'
    'Specify either one number, applicable to all onsets,'
    'or a list of amplitudes for each onset'}';

onset_type_info         = cfg_branch;
onset_type_info.tag     = 'onset_type_info';
onset_type_info.name    = 'New onset type (applies to this session)';
onset_type_info.val     = {onset_name onset_times onset_durations onset_amplitude};
onset_type_info.help    = {'Specify information for this onset type.'}';

onset_type_list         = cfg_repeat;
onset_type_list.tag     = 'onset_type_list';
onset_type_list.name    = 'New session';
onset_type_list.help    = {'Note: the use of this function always breaks the dependency to the next module.'
    'Therefore, after all the onset types have been specified, one needs to always relink the IOI.mat dependency'
    'in the next module (if there is one)'}';
onset_type_list.values  = {onset_type_info};
onset_type_list.num     = [0 Inf];

onsets_by_session         = cfg_repeat;
onsets_by_session.tag     = 'onsets_by_session';
onsets_by_session.name    = 'New subject';
onsets_by_session.help    = {'For each session selected above, specify the onsets for this session'
    'They will be applied to each subject.'
    'If there are several selected sessions and onsets are specified only for one session, '
    'they will be applied to all selected sessions.'
    'Note: the use of this function always breaks the dependency to the next module.'
    'Therefore, after all the onset types have been specified, one needs to always relink the IOI.mat dependency'
    'in the next module (if there is one)'}';
onsets_by_session.values  = {onset_type_list};
onsets_by_session.num     = [0 Inf];

onsets_by_subject         = cfg_repeat;
onsets_by_subject.tag     = 'onsets_by_subject';
onsets_by_subject.name    = 'Subjects';
onsets_by_subject.help    = {'For each subject, specify the onsets to be used'
      'If there are several subjects and onsets are specified only for one subject, '
    'they will be applied to all subjects.'
    'Note: the use of this function always breaks the dependency to the next module.'
    'Therefore, after all the onset types have been specified, one needs to always relink the IOI.mat dependency'
    'in the next module (if there is one)'}';
onsets_by_subject.values  = {onsets_by_session};
onsets_by_subject.num     = [0 Inf];

manual_stims         = cfg_branch;
manual_stims.tag     = 'manual_stims';
manual_stims.name    = 'Stimulations manually defined';
manual_stims.val     = {onsets_by_subject};
manual_stims.help    = {'Stimulations are defined by user inputs'
    'This defines IOI.Sess().U, which is generated by user inputs.'
    'If nothing is specified at this stage, the user will be queried during code execution.'
    'This is the most generic mechanism, to use if onsets are different for different subjects.'}';

default_stims         = cfg_branch;
default_stims.tag     = 'default_stims';
default_stims.name    = 'Default stimulations';
default_stims.val     = {};
default_stims.help    = {'The names, onsets, durations fields in IOI.sess_res will be used'};

stim_choice        = cfg_choice;
stim_choice.name   = 'Choose onset selection method';
stim_choice.tag    = 'stim_choice';
stim_choice.values = {default_stims,electro_stims,manual_stims};
stim_choice.val    = {default_stims};
stim_choice.help   = {'Choose onset selection method'
       'Manual: the user will be queried to enter onset times and durations.'
       'With default stims, nothing will be done, as stimulation times should '
       'have been found during the initial treatment of the images;'
       'However a check will be made that onsets are available.'
       'For electrophysiology, events will be detected on the electrophysiology file.'}';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Executable Branch
create_onsets1      = cfg_exbranch;       % This is the branch that has information about how to run this module
create_onsets1.name = 'Create onsets';             % The display name
create_onsets1.tag  = 'create_onsets1'; %Very important: tag is used when calling for execution
create_onsets1.val  = {IOImat elDir redo1 IOImatCopyChoice session_choice ...
     stim_choice};    % The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs
create_onsets1.prog = @ioi_create_onsets_run;  % A function handle that will be called with the harvested job to run the computation
create_onsets1.vout = @ioi_cfg_vout_create_onsets; % A function handle that will be called with the harvested job to determine virtual outputs
create_onsets1.help = {'Create onsets.'};

return

%make IOI.mat available as a dependency
function vout = ioi_cfg_vout_create_onsets(job)
vout = cfg_dep;                     % The dependency object
vout.sname      = 'IOI.mat';       % Displayed dependency name
vout.src_output = substruct('.','IOImat'); %{1}); %,'IOImat');
%substruct('()',{1}); % The output subscript reference. This could be any reference into the output variable created during computation
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
