function create_roi1 = ioi_create_roi_cfg
%_______________________________________________________________________
% Copyright (C) 2010 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%______________________________________________________________________
% Read IOI Multispectral files
% ---------------------------------------------------------------------

IOImat = ioi_cfg_IOImat(1);
redo1 = ioi_dfg_redo(0);
ROImat = ioi_cfg_ROImat(1); %oops needs to be done for this unit
%(specifically: allow user to specify a new location for ROIs)

RemovePreviousROI      = cfg_menu;
RemovePreviousROI.tag  = 'RemovePreviousROI';
RemovePreviousROI.name = 'Treatment of previous ROIs';
RemovePreviousROI.labels = {'Keep','Remove'};
RemovePreviousROI.values = {0,1};
RemovePreviousROI.val  = {0};
RemovePreviousROI.help = {'If previous ROIs had been selected, choose whether to keep'
    'or to remove this information.'}';

select_names      = cfg_menu;
select_names.tag  = 'select_names';
select_names.name = 'Select names for ROIs';
select_names.labels = {'Yes','No'};
select_names.values = {1,0};
select_names.val  = {1};
select_names.help = {'Option for user to manually enter names of ROIs.'
    'If No is selected, ROI names will be a number (enumeration).'}';

IOImatCopyChoice = ioi_cfg_IOImatCopyChoice('ROI');
        
ArrayROI         = cfg_entry;
ArrayROI.name    = 'Size of array of ROIs';
ArrayROI.tag     = 'ArrayROI';       
ArrayROI.strtype = 'r';
ArrayROI.val{1}    = [2 3];
ArrayROI.num     = [1 2];     
ArrayROI.help    = {'Enter array (number of rows by number of columns) '
    'for automatically generated array of ROIs.'
    'Make sure the product of rows by columns (total number of ROIs)'
    'is less than 100.'
    'Names of ROIs will include array positions, increasing from down to up'
    'vertically and from left to right horizontally.'}'; 

AutoROI         = cfg_branch;
AutoROI.tag     = 'AutoROI';
AutoROI.name    = 'Automatic ROI selection'; 
AutoROI.val     = {ArrayROI};
AutoROI.help    = {'Automatic ROI selection.'}';
        
ManualROI         = cfg_branch;
ManualROI.tag     = 'ManualROI';
ManualROI.name    = 'Manual ROI selection: graphical tool'; 
ManualROI.val     = {};
ManualROI.help    = {'Manual ROI selection: graphical tool.'}';

ManualEnterROI         = cfg_branch;
ManualEnterROI.tag     = 'ManualEnterROI';
ManualEnterROI.name    = 'Manual ROI selection: coordinate entry'; 
ManualEnterROI.val     = {};
ManualEnterROI.help    = {'Manual ROI selection: coordinate entry.'}';

AutoROIchoice           = cfg_choice;
AutoROIchoice.name      = 'Choose ROI generation method';
AutoROIchoice.tag       = 'AutoROIchoice';
AutoROIchoice.values    = {ManualROI ManualEnterROI AutoROI}; 
AutoROIchoice.val       = {ManualROI}; 
AutoROIchoice.help      = {'Choose whether to generate ROI manually or'
        'automatically'}'; 
        
% Executable Branch
create_roi1      = cfg_exbranch;       % This is the branch that has information about how to run this module
create_roi1.name = 'Create ROI';             % The display name
create_roi1.tag  = 'create_roi1'; %Very important: tag is used when calling for execution
create_roi1.val  = {IOImat redo1 RemovePreviousROI IOImatCopyChoice select_names AutoROIchoice};    % The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs
create_roi1.prog = @ioi_create_roi_run;  % A function handle that will be called with the harvested job to run the computation
create_roi1.vout = @ioi_cfg_vout_create_roi; % A function handle that will be called with the harvested job to determine virtual outputs
create_roi1.help = {'Create regions of interest.'};

return

%make IOI.mat available as a dependency
function vout = ioi_cfg_vout_create_roi(job)
vout = cfg_dep;                     % The dependency object
vout.sname      = 'IOI.mat';       % Displayed dependency name
vout.src_output = substruct('.','IOImat'); %{1}); %,'IOImat');
%substruct('()',{1}); % The output subscript reference. This could be any reference into the output variable created during computation
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
