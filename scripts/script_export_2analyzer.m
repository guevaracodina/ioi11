%% script_export_2analyzer
%% Adding nirs10 toolbox to the path
% Run spm8 shortcut first!
% Run SPM 8
% cd c:\spm8
% addpath(genpath('c:\spm8'))
% spm fmri
clear; close all; clc
addpath(genpath('C:\spm8\toolbox\nirs10')); % After running spm8 shortcut

%% Parent folders with IOImat
pathNameParentList{2} = 'C:\Edgar\Data\IOIResults\j03\ROItest';


%% Folders with physiological data (in the same order as the IOImat list)
pathNamePhysioList{2} = 'C:\Edgar\Data\IOIData20141211\Physio Monitoring\121423';

%% Batch processing of data to be exported to analyzer 2. Do not forget to create the header files
OK = ioi_export_2analyzer (pathNameParent, pathNamePhysio)

% EOF