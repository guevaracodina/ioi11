function filtdown1 = ioi_filtdown_cfg
% Graphical interface configuration function for ioi_filtdown_cfg
%_______________________________________________________________________
% Copyright (C) 2010 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%______________________________________________________________________


IOImat          = cfg_files; %Select IOI.mat for this subject 
IOImat.name     = 'Select IOI.mat'; % The displayed name
IOImat.tag      = 'IOImat';       %file names
IOImat.filter   = 'mat';
IOImat.ufilter  = '^IOI.mat$';    
IOImat.num      = [1 Inf];     % Number of inputs required 
IOImat.help     = {'Select IOI.mat for this subject.'}; % help text displayed

redo1           = cfg_menu;
redo1.tag       = 'force_redo';
redo1.name      = 'Force processing';
redo1.labels    = {'False','True'};
redo1.values    = {0,1};
redo1.val       = {0};
redo1.help      = {'Force redoing this processing even when it has been done already'};

% Select mask file
maskFile            = cfg_files;
maskFile.tag        = 'maskFile';
maskFile.name       = 'Mask File';
maskFile.filter     = 'nifti';
maskFile.ufilter    = '.*_brainmask.nii$';
maskFile.num        = [1 1];
maskFile.help       = {'Select file containing the brain mask.'};

% Choose type of time trace to be used as a contrast [HbO/HbR/Flow]
HbOcontrast         = cfg_branch;
HbOcontrast.tag     = 'HbOcontrast';
HbOcontrast.name    = 'HbO time trace'; 
HbOcontrast.help    = {'HbO time trace will be used as a contrast'};

HbRcontrast         = cfg_branch;
HbRcontrast.tag     = 'HbRcontrast';
HbRcontrast.name    = 'HbR time trace'; 
HbRcontrast.help    = {'HbR time trace will be used as a contrast'};

Flowcontrast        = cfg_branch;
Flowcontrast.tag    = 'Flowcontrast';
Flowcontrast.name   = 'Blood flow time trace'; 
Flowcontrast.help   = {'Blood flow time trace will be used as a contrast'};
        
% Choose type of time trace to be used as a contrast [HbO/HbR/Flow]
contrastChoice          = cfg_choice;
contrastChoice.name     = 'Choose time trace to be used as a contrast';
contrastChoice.tag      = 'contrastChoice';
contrastChoice.values   = {HbOcontrast HbRcontrast Flowcontrast}; 
contrastChoice.val      = {HbOcontrast}; 
contrastChoice.help     = {'Choose what type of time trace [HbO/HbR/Flow]'
        ' will be used as a contrast.'}';

% Bandpass filtering
BPFfreq         = cfg_entry;
BPFfreq.name    = 'Band-pass filter cutoff frequencies';
BPFfreq.tag     = 'BPFfreq';       
BPFfreq.strtype = 'r';
BPFfreq.val{1}  = [0.009 0.08];
BPFfreq.num     = [1 2];     
BPFfreq.help    = {'Enter Wn, a two-element vector, Wn = [W1 W2] for the '
    'bandpass filter with passband  W1 < W < W2'}';

% Remove global mean signal
removeMean           = cfg_menu;
removeMean.tag       = 'removeMean';
removeMean.name      = 'Remove global mean signal';
removeMean.labels    = {'False','True'};
removeMean.values    = {0,1};
removeMean.val       = {1};
removeMean.help      = {'Remove global mean signal from the non-masked brain pixels'};

% Executable Branch
filtdown1      = cfg_exbranch;       % This is the branch that has information about how to run this module
filtdown1.name = 'Generate Network Mask';             % The display name
filtdown1.tag  = 'filtdown1'; %Very important: tag is used when calling for execution
filtdown1.val  = {IOImat redo1 maskFile contrastChoice BPFfreq removeMean};    % The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs
filtdown1.prog = @ioi_filtdown_run;  % A function handle that will be called with the harvested job to run the computation
filtdown1.vout = @ioi_cfg_vout_filtdown; % A function handle that will be called with the harvested job to determine virtual outputs
filtdown1.help = {'Temporal band-pass filtering and downsampling of a given time trace [HbO/HbR/Flow] to be used as a contrast.'};

return

% make IOI.mat available as a dependency
function vout = ioi_cfg_vout_filtdown(~)
vout = cfg_dep;                     % The dependency object
vout.sname      = 'IOI.mat';       % Displayed dependency name
vout.src_output = substruct('.','IOImat'); %{1}); %,'IOImat');
%substruct('()',{1}); % The output subscript reference. This could be any reference into the output variable created during computation
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
