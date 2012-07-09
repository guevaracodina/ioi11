function filtdown1 = ioi_filtdown_cfg
% Graphical interface configuration function for ioi_filtdown_cfg
%_______________________________________________________________________
% Copyright (C) 2010 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%______________________________________________________________________

% Select IOI.mat
IOImat              = ioi_dfg_IOImat(1);
% Force processing
redo1               = ioi_dfg_redo(0);
% IOI copy/overwrite method
IOImatCopyChoice    = ioi_dfg_IOImatCopyChoice('FiltNDown');
% Choose ROI selection method (all/selected)
ROI_choice          = ioi_dfg_ROI_choice;
% Choose session selection method (all/selected)
session_choice      = ioi_dfg_session_choice;
% Colors to include (OD,HbO,HbR,HbT,Flow)
IC                  = ioi_dfg_include_colors(0,1,1,1,1);

% Bandpass filtering
BPFfreq             = cfg_entry;
BPFfreq.name        = 'Band-pass filter cutoff frequencies';
BPFfreq.tag         = 'BPFfreq';       
BPFfreq.strtype     = 'r';
BPFfreq.val         = {[0.009 0.08]};
BPFfreq.num         = [1 2];     
BPFfreq.help        = {'Enter Wn in Hz, a two-element vector, Wn = [W_1 W_2] for the bandpass filter with passband  W_1 < W < W_2'}';

% Downsampling frequency
downFreq            = cfg_entry;
downFreq.name       = 'Downsampling frequency in Hz';   % The displayed name
downFreq.tag        = 'downFreq';                       % file names
downFreq.strtype    = 'r';                              % Real numbers
downFreq.num        = [1 1];                            % Number of inputs required
downFreq.val        = {1};                              % Default value
downFreq.help       = {'Enter downsampling frequency in Hz. A target sampling frequency will be generated, which may however be only approximately equal to the specified downsampling frequency, but it will correspond to the actual frequency of selecting every Nth point'};

% Filter and downsample whole images
wholeImage          = cfg_menu;
wholeImage.tag      = 'wholeImage';
wholeImage.name     = 'Filter and downsample whole image time-series';
wholeImage.labels   = {'Yes','No'};
wholeImage.values   = {1,0};
wholeImage.val      = {1};
wholeImage.help     = {'Filter and downsample whole image time-series. It creates a new sub-folder for each session'};

% Executable Branch
filtdown1           = cfg_exbranch; % This is the branch that has information about how to run this module
filtdown1.name      = 'Temporal filtering and downsampling of seeds & whole image series';             % The display name
filtdown1.tag       = 'filtdown1'; %Very important: tag is used when calling for execution
filtdown1.val       = {IOImat redo1 IOImatCopyChoice ROI_choice session_choice IC BPFfreq downFreq wholeImage};    % The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs
filtdown1.prog      = @ioi_filtdown_run;  % A function handle that will be called with the harvested job to run the computation
filtdown1.vout      = @ioi_cfg_vout_filtdown; % A function handle that will be called with the harvested job to determine virtual outputs
filtdown1.help      = {'Temporal band-pass filtering and downsampling of a given time trace [HbO/HbR/Flow], either on a seed or on the whole image series.'};

return

% Make IOI.mat available as a dependency
function vout = ioi_cfg_vout_filtdown(job)
vout                = cfg_dep; % The dependency object
vout.sname          = 'IOI.mat'; % Displayed dependency name
vout.src_output     = substruct('.','IOImat'); %{1}); %,'IOImat');
%substruct('()',{1}); % The output subscript reference. This could be any reference into the output variable created during computation
vout.tgt_spec       = cfg_findspec({{'filter','mat','strtype','e'}});

% EOF

