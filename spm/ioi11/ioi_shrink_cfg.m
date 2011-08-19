function shrink1 = ioi_shrink_cfg
%_______________________________________________________________________
% Copyright (C) 2010 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%______________________________________________________________________
% Shrink specific .bin files
% ---------------------------------------------------------------------

%Select top level directory with .bin files
files1         = cfg_files;
files1.tag     = 'files1';
files1.name    = 'Select .bin files to shrink';
%files1.filter = ''; 
files1.ufilter = '.bin';
files1.num     = [1 Inf];
files1.help    = {'Select .bin files to shrink.'}';

shrink_x      = cfg_entry;
shrink_x.tag  = 'shrink_x';
shrink_x.name = 'Shrink factor for x dimension';
shrink_x.strtype  = 'i';
shrink_x.num = [1 1];
shrink_x.val{1}  = 2;
shrink_x.help = {'Data reduction factor in x.'};

shrink_y      = cfg_entry;
shrink_y.tag  = 'shrink_y';
shrink_y.name = 'Shrink factor for y dimension';
shrink_y.strtype  = 'i';
shrink_y.num = [1 1];
shrink_y.val{1}  = 2;
shrink_y.help = {'Data reduction factor in y.'};

% Executable Branch
shrink1      = cfg_exbranch;       % This is the branch that has information about how to run this module
shrink1.name = 'Shrink specific .bin files';             % The display name
shrink1.tag  = 'shrink1'; %Very important: tag is used when calling for execution
shrink1.val  = {files1 shrink_x shrink_y};    % The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs
shrink1.prog = @ioi_shrink_run;  % A function handle that will be called with the harvested job to run the computation
%shrink1.vout = @ioi_cfg_vout_msioi; % A function handle that will be called with the harvested job to determine virtual outputs
shrink1.help = {'Utility to shrink specific .bin files and replace '
    'them in the same spot '
    'Use only to handle unusual situations.'}';
return
