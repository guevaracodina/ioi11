function create_roi1 = ioi_create_roi_cfg
% Graphical interface configuration function for ioi_create_roi_run
%_______________________________________________________________________________
% Copyright (C) 2010 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________

% Select IOI.mat
IOImat                  = ioi_dfg_IOImat(1);
% Force processing
redo1                   = ioi_dfg_redo(0);
% Choose ROI selection method (all/selected)
ROImat                  = ioi_dfg_ROImat(1); %oops needs to be done for this unit
%(specifically: allow user to specify a new location for ROIs)
% IOI copy/overwrite method
IOImatCopyChoice        = ioi_dfg_IOImatCopyChoice('ROI');

% Keep/delete previous ROIs
RemovePreviousROI       = cfg_menu;
RemovePreviousROI.tag   = 'RemovePreviousROI';
RemovePreviousROI.name  = 'Treatment of previous ROIs';
RemovePreviousROI.labels = {'Keep','Remove'};
RemovePreviousROI.values = {0,1};
RemovePreviousROI.val   = {0};
RemovePreviousROI.help  = {'If previous ROIs had been selected, choose whether to keep'
    'or to remove this information.'}';

% Give names for ROIs
select_names            = cfg_menu;
select_names.tag        = 'select_names';
select_names.name       = 'Select names for ROIs';
select_names.labels     = {'Yes','No'};
select_names.values     = {1,0};
select_names.val        = {1};
select_names.help       = {'Option for user to manually enter names of ROIs.'
    'If No is selected, ROI names will be a number (enumeration).'}';

% Array of automatic ROIs
ArrayROI                = cfg_entry;
ArrayROI.name           = 'Size of array of ROIs';
ArrayROI.tag            = 'ArrayROI';       
ArrayROI.strtype        = 'r';
ArrayROI.val{1}         = [2 3];
ArrayROI.num            = [1 2];     
ArrayROI.help           = {'Enter array (number of rows by number of columns) '
    'for automatically generated array of ROIs.'
    'Make sure the product of rows by columns (total number of ROIs)'
    'is less than 100.'
    'Names of ROIs will include array positions, increasing from down to up'
    'vertically and from left to right horizontally.'}'; 

AutoROI                 = cfg_branch;
AutoROI.tag             = 'AutoROI';
AutoROI.name            = 'Automatic ROI selection'; 
AutoROI.val             = {ArrayROI};
AutoROI.help            = {'Automatic ROI selection.'}';
        
ManualROI               = cfg_branch;
ManualROI.tag           = 'ManualROI';
ManualROI.name          = 'Manual ROI selection: graphical tool'; 
ManualROI.val           = {};
ManualROI.help          = {'Manual ROI selection: specify polygonal ROI/seed.'}';

ManualROIspline         = cfg_branch;
ManualROIspline.tag     = 'ManualROIspline';
ManualROIspline.name    = 'Manual ROI selection: spline graphical tool'; 
ManualROIspline.val     = {};
ManualROIspline.help    = {'Manual ROI selection: specify spline ROI/seed.'}';

ManualEnterROI          = cfg_branch;
ManualEnterROI.tag      = 'ManualEnterROI';
ManualEnterROI.name     = 'Manual ROI selection: coordinate entry'; 
ManualEnterROI.val      = {};
ManualEnterROI.help     = {'Manual ROI selection: coordinate entry.'}';

ManualROIradius         = cfg_entry;
ManualROIradius.name    = 'Radius size of ROIs/seeds';
ManualROIradius.tag     = 'ManualROIradius';       
ManualROIradius.strtype = 'r';
ManualROIradius.val{1}  = 7;                    % Default value
ManualROIradius.num     = [1 1];     
ManualROIradius.help    = {'Enter radius of the ROIs/seeds in pixels.'}'; 

pointNclickROI          = cfg_branch;
pointNclickROI.tag      = 'pointNclickROI';
pointNclickROI.name     = 'Manual ROI selection: point & click'; 
pointNclickROI.val      = {ManualROIradius};
pointNclickROI.help     = {'Manual ROI selection: point & click the center of circular ROI/seed. Usually'
                                '  1,2: Frontal'
                                '  3,4: Motor'
                                '  5,6: Cingulate'
                                '  7,8: Somatosensory'
                                ' 9,10: Retrosplenial'
                                '11,12: Visual'}';

AutoROIchoice           = cfg_choice;
AutoROIchoice.name      = 'Choose ROI generation method';
AutoROIchoice.tag       = 'AutoROIchoice';
AutoROIchoice.values    = {ManualROI ManualROIspline ManualEnterROI pointNclickROI AutoROI}; 
AutoROIchoice.val       = {ManualROI}; 
AutoROIchoice.help      = {'Choose whether to generate ROI manually or'
        'automatically'}'; 

% Display Brain Mask    
displayBrainmask        = cfg_menu;
displayBrainmask.tag    = 'displayBrainmask';
displayBrainmask.name   = 'Display Brain Mask';
displayBrainmask.labels = {'Yes','No'};
displayBrainmask.values = {1,0};
displayBrainmask.val    = {0}; % Default value = 0
displayBrainmask.help   = {'Display network mask containing only brain pixels'};

use_gray_contrast      = cfg_menu;
use_gray_contrast.tag  = 'use_gray_contrast';
use_gray_contrast.name = 'Use gray contrast';
use_gray_contrast.labels = {'False','True'};
use_gray_contrast.values = {0,1};
use_gray_contrast.val  = {1};
use_gray_contrast.help = {'Use gray contrast.'}';

% Choose whether interface should ask whether to use a previously saved list of ROI
SelectPreviousROI       = cfg_menu;
SelectPreviousROI.tag   = 'SelectPreviousROI';
SelectPreviousROI.name  = 'Select whether interface should ask for a possible previous list of ROIs';
SelectPreviousROI.labels = {'No','Yes'};
SelectPreviousROI.values = {0,1};
SelectPreviousROI.val   = {0};
SelectPreviousROI.help  = {'If this option is selected, then before manual selection of ROIs,'
    'The interface will ask the user if they want to use an ROI list from elsewhere.'
    'This ROI list can then be selected by locating the IOI.mat structure that contains the information.'}';

% Executable Branch
create_roi1             = cfg_exbranch;       % This is the branch that has information about how to run this module
create_roi1.name        = 'Create ROI/seed';             % The display name
create_roi1.tag         = 'create_roi1'; %Very important: tag is used when calling for execution
create_roi1.val         = {IOImat redo1 RemovePreviousROI IOImatCopyChoice ...
    select_names AutoROIchoice displayBrainmask use_gray_contrast SelectPreviousROI};    % The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs
create_roi1.prog        = @ioi_create_roi_run;  % A function handle that will be called with the harvested job to run the computation
create_roi1.vout        = @ioi_cfg_vout_create_roi; % A function handle that will be called with the harvested job to determine virtual outputs
create_roi1.help        = {'Create regions of interest/seeds.'};

return

%make IOI.mat available as a dependency
function vout = ioi_cfg_vout_create_roi(job)
vout = cfg_dep;                     % The dependency object
vout.sname      = 'IOI.mat';       % Displayed dependency name
vout.src_output = substruct('.','IOImat'); %{1}); %,'IOImat');
%substruct('()',{1}); % The output subscript reference. This could be any reference into the output variable created during computation
vout.tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});

% EOF
