function conc1 = ioi_concentrations_cfg
%_______________________________________________________________________________
% Copyright (C) 2011 LIOM Laboratoire d'Imagerie Optique et Mol�culaire
%                    �cole Polytechnique de Montr�al
%_______________________________________________________________________________
% Read IOI Multispectral files
% ---------------------------------------------------------------------

IOImat = ioi_dfg_IOImat(1);
redo1 = ioi_dfg_redo(0);

pathlength1      = cfg_menu;
pathlength1.tag  = 'pathlength';
pathlength1.name = 'Pathlength Factor';
pathlength1.labels = {'Kohl','Dunn'};
pathlength1.values = {'Kohl','Dunn'};
pathlength1.def  = @(val)ioi_get_defaults('conc1.pathlength', val{:});
pathlength1.help = {'Which pathlength factor to use from literature'};

camera1      = cfg_menu;
camera1.tag  = 'camera_corr';
camera1.name = 'Camera correction';
camera1.labels = {'False','True'};
camera1.values = {0,1};
camera1.def  = @(val)ioi_get_defaults('conc1.camera_correction', val{:});
camera1.help = {'Correct for camera wavelength sensitivity'};

led1      = cfg_menu;
led1.tag  = 'leds_conv';
led1.name = 'Leds Convolution';
led1.labels = {'False','True'};
led1.values = {0,1};
led1.def  = @(val)ioi_get_defaults('conc1.leds_spectra', val{:});
led1.help = {'Do we convolve by LEDs spectra curves to get effective epsilons?'};

basehbt1      = cfg_entry;
basehbt1.tag  = 'HbT0';
basehbt1.name = 'HbT0';
basehbt1.strtype  = 'r';
basehbt1.num = [1 1];
basehbt1.def  = @(val)ioi_get_defaults('conc1.baseline_hbt', val{:});
basehbt1.help = {'Baseline total hemoglobin, in microM, default 100 microM.'};

basehbr1      = cfg_entry;
basehbr1.tag  = 'HbR0';
basehbr1.name = 'HbR0';
basehbr1.strtype  = 'r';
basehbr1.num = [1 1];
basehbr1.def  = @(val)ioi_get_defaults('conc1.baseline_hbr', val{:});
basehbr1.help = {'Baseline deoxyhemoglobin, in microM, default 40 microM.'};

basehbo1      = cfg_entry;
basehbo1.tag  = 'HbO0';
basehbo1.name = 'HbO0';
basehbo1.strtype  = 'r';
basehbo1.num = [1 1];
basehbo1.def  = @(val)ioi_get_defaults('conc1.baseline_hbo', val{:});
basehbo1.help = {'Baseline oxyhemoglobin, in microM, default 60 microM.'};


configuration         = cfg_branch;
configuration.tag     = 'configuration';
configuration.name    = 'Configuration options';
configuration.val     = {pathlength1 camera1 led1 basehbt1 basehbo1 basehbr1};
configuration.help    = {'Select values.'};

IOImatCopyChoice        = ioi_dfg_IOImatCopyChoice('Conc');

MemoryManagementMenu    = ioi_dfg_MemoryManagement;

session_choice          = ioi_dfg_session_choice;

%%%%%%%%%%%%%%%%%%%%
% Normalization choice
%%%%%%%%%%%%%%%%%%%%%%%%%%

current_sessions         = cfg_branch;
current_sessions.tag     = 'current_sessions';
current_sessions.name    = 'Current session';
current_sessions.val     = {};
current_sessions.help    = {'Use normalization calculated previously' 
    'for this session'}';

selected_norm_session      = cfg_entry;
selected_norm_session.tag  = 'selected_norm_session';
selected_norm_session.name = 'Enter session to use for normalization';
selected_norm_session.strtype  = 'r';
selected_norm_session.num = [1 1];
selected_norm_session.val{1} = 1;
selected_norm_session.help = {'Enter session to use for normalization.'
    'The median for this session will be used for all the other sessions.'}';

select_norm_session         = cfg_branch;
select_norm_session.tag     = 'select_norm_session';
select_norm_session.name    = 'Select session median to use for normalization';
select_norm_session.val     = {selected_norm_session};
select_norm_session.help    = {'Select session median to use for normalization'};

normalization_choice        = cfg_choice;
normalization_choice.name   = 'Choose session normalization method';
normalization_choice.tag    = 'normalization_choice';
normalization_choice.values = {current_sessions,select_norm_session};
normalization_choice.val    = {current_sessions};
normalization_choice.help   = {'Choose session selection method'}';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RemoveRGY      = cfg_menu;
RemoveRGY.tag  = 'RemoveRGY';
RemoveRGY.name = 'Remove RGY nifti images';
RemoveRGY.labels = {'Yes','No'};
RemoveRGY.values = {1,0};
RemoveRGY.val  = {1};
RemoveRGY.help = {'After concentration concentrations are obtained'
    'Optical intensity images (RGY) are usually no longer required.'}';

% Executable Branch
conc1      = cfg_exbranch;       % This is the branch that has information about how to run this module
conc1.name = 'Compute Concentrations';             % The display name
conc1.tag  = 'conc1'; %Very important: tag is used when calling for execution
conc1.val  = {IOImat redo1 IOImatCopyChoice normalization_choice configuration ...
    MemoryManagementMenu session_choice RemoveRGY};    % The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs
conc1.prog = @ioi_concentrations_run;  % A function handle that will be called with the harvested job to run the computation
conc1.vout = @ioi_cfg_vout_concentrations; % A function handle that will be called with the harvested job to determine virtual outputs
conc1.help = {'Concentration computations.'};

return

%make IOI.mat available as a dependency
function vout = ioi_cfg_vout_concentrations(job)
vout = cfg_dep;                     % The dependency object
vout.sname      = 'IOI.mat';       % Displayed dependency name
vout.src_output = substruct('.','IOImat'); %{1}); %,'IOImat');
%substruct('()',{1}); % The output subscript reference. This could be any reference into the output variable created during computation
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
