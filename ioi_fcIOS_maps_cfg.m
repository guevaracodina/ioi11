function fcIOS_maps         = ioi_fcIOS_maps_cfg
% Graphical interface configuration function for ioi_fcIOS_maps_run
%_______________________________________________________________________________
% Copyright (C) 2012 LIOM Laboratoire d'Imagerie Optique et Mol�culaire
%                    �cole Polytechnique de Montr�al
%_______________________________________________________________________________

% Select IOI.mat
IOImat                      = ioi_dfg_IOImat(1);
% IOI copy/overwrite method
IOImatCopyChoice            = ioi_dfg_IOImatCopyChoice('fcIOS_maps');
% Choose ROI selection method (all/selected)
ROI_choice                  = ioi_dfg_ROI_choice;
% Colors to include (OD,HbO,HbR,HbT,Flow,CMRO2)
IC                          = ioi_dfg_include_colors(0,1,1,1,1);
% Choose session selection method (all/selected)
session_choice              = ioi_dfg_session_choice;
% Default choice to draw a circle to indicate seed position and size
draw_circle_def_on          = true;

% Anonymous functions used later. They can be evaluated in transM
rotx = @(theta) [1 0 0 0; 0 cos(theta) -sin(theta) 0; 0 sin(theta) cos(theta) 0; 0 0 0 1];
roty = @(theta) [cos(theta) 0 sin(theta) 0; 0 1 0 0; -sin(theta) 0 cos(theta) 0; 0 0 0 1];
rotz = @(theta) [cos(theta) -sin(theta) 0 0; sin(theta) cos(theta) 0 0; 0 0 1 0; 0 0 0 1];
translate = @(a,b) [1 0 a 0; 0 1 b 0; 0 0 1 0; 0 0 0 1];

% Select directory to save global results
parent_results_dir          = cfg_files;
parent_results_dir.tag      = 'parent_results_dir';
parent_results_dir.name     = 'Top directory to save group results';
parent_results_dir.filter   = 'dir';
parent_results_dir.num      = [1 1];
parent_results_dir.help     = {'Select the directory where figures will be saved.'}';

% Figure size
figSize                     = cfg_entry;
figSize.tag                 = 'figSize';
figSize.name                = 'Figure size';
figSize.strtype             = 'r';
figSize.num                 = [1 2];
figSize.val{1}              = [1 1];
figSize.help                = {'Enter figure size in inches [Width Height].'};

% Figure resolution
figRes                      = cfg_entry;
figRes.tag                  = 'figRes';
figRes.name                 = 'Figure resolution';
figRes.strtype              = 'r';
figRes.num                  = [1 1];
figRes.val{1}               = 1200;
figRes.help                 = {'Enter figure resolution in dpi [150-1200]'};

% Colormap range
figRange                    = cfg_entry;
figRange.tag                = 'figRange';
figRange.name               = 'Colormap range';
figRange.strtype            = 'r';
figRange.num                = [1 2];
figRange.val{1}             = [-1 1];
figRange.help               = {'Enter colormap range. For correlation maps default is [-1 1]'};

% Alpha transparency values
figAlpha                    = cfg_entry;
figAlpha.tag                = 'figAlpha';
figAlpha.name               = 'Transparency';
figAlpha.strtype            = 'r';
figAlpha.num                = [1 1];
figAlpha.val{1}             = 0.8;
figAlpha.help               = {'Enter colormap transparency, between 0 and 1'};

% Intensity for the anatomical background
figIntensity                = cfg_entry;
figIntensity.tag            = 'figIntensity';
figIntensity.name           = 'Intensity';
figIntensity.strtype        = 'r';
figIntensity.num            = [1 1];
figIntensity.val{1}         = 1;
figIntensity.help           = {'Enter background (anatomical) intensity, between 0 and 1'};

% Colormap to use
figCmap                     = cfg_entry;
figCmap.tag                 = 'figCmap';
figCmap.name                = 'Colormap';
figCmap.strtype             = 'e';
figCmap.num                 = [Inf 3];
figCmap.val{1}              = jet(256);
figCmap.help                = {'Enter colormap to use. e.g. type jet(256), Input is evaluated'};

% ------------------------------------------------------------------------------
% Seed positions and sizes will be shown with black circles.
% ------------------------------------------------------------------------------
circleLW                    = cfg_entry;
circleLW.tag                = 'circleLW';
circleLW.name               = 'Circle LineWidth';
circleLW.val                = {0.8};
circleLW.strtype            = 'r';
circleLW.num                = [1 1];
circleLW.help               = {'Circle Line Width'};

circleLS                    = cfg_entry;
circleLS.tag                = 'circleLS';
circleLS.name               = 'Circle LineStyle';
circleLS.val                = {'-'};
circleLS.strtype            = 's';
circleLS.num                = [1 2];
circleLS.help               = {'Circle Line Style'};

drawCircle_On               = cfg_branch;
drawCircle_On.tag           = 'drawCircle_On';
drawCircle_On.name          = 'Draw circle';
drawCircle_On.val           = {circleLW circleLS};
drawCircle_On.help          = {'Draw a circle in seed position'};

drawCircle_Off              = cfg_branch;
drawCircle_Off.tag          = 'drawCircle_Off';
drawCircle_Off.name         = 'Draw nothing';
drawCircle_Off.val          = {};
drawCircle_Off.help         = {'Draw nothing'};

drawCircle                  = cfg_choice;
drawCircle.tag              = 'drawCircle';
drawCircle.name             = 'Draw circle in seed';
drawCircle.values           = {drawCircle_On drawCircle_Off};
if draw_circle_def_on
    drawCircle.val          = {drawCircle_On};
else
    drawCircle.val          = {drawCircle_Off};
end
drawCircle.help             = {'Seed positions and sizes will be shown with black circles.'};
% ------------------------------------------------------------------------------

% transM Reorientation Matrix
transM                      = cfg_entry;
transM.tag                  = 'transM';
transM.name                 = 'Reorientation Matrix';
transM.strtype              = 'e';
transM.val{1}               = rotz(pi);
transM.num                  = [4 4];
transM.help                 = {
                            'Enter a valid 4x4 matrix for reorientation.'
                            ''
                            'Example:'
                            'This will rotate the images to have rostral orientation up.'
                            ''
                            'rotz(pi)'
                            ''
                            'rotx(theta), roty(theta) and translate(a,b) can also be evaluated'
                            }';

% Executable Branch
fcIOS_maps                  = cfg_exbranch; % This is the branch that has information about how to run this module
fcIOS_maps.name             = 'Plot multiple correlation maps'; % The display name
fcIOS_maps.tag              = 'fcIOS_maps'; %Very important: tag is used when calling for execution
fcIOS_maps.val              = {IOImat IOImatCopyChoice ROI_choice ...
    IC session_choice parent_results_dir figSize figRes figRange figAlpha ...
    figIntensity figCmap drawCircle transM};    % The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs
fcIOS_maps.prog             = @ioi_fcIOS_maps_run; % A function handle that will be called with the harvested job to run the computation
fcIOS_maps.vout             = @ioi_cfg_vout_fcIOS_maps; % A function handle that will be called with the harvested job to determine virtual outputs
fcIOS_maps.help             = {'Plots multiple correlation maps at the same scale. Ideal to make a mosaique figure.'};


return

% Make IOI.mat available as a dependency
function vout = ioi_cfg_vout_fcIOS_maps(job)
vout                        = cfg_dep;                  % The dependency object
vout.sname                  = 'IOI.mat';                % Displayed dependency name
vout.src_output             = substruct('.','IOImat');  %{1}); %,'IOImat');
vout.tgt_spec               = cfg_findspec({{'filter','mat','strtype','e'}});

% EOF
