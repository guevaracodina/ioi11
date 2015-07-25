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
pathNameParentList{1} = 'C:\Edgar\Data\IOIResults\j02\ROItest';
pathNameParentList{2} = 'C:\Edgar\Data\IOIResults\j03\ROItest';
pathNameParentList{3} = 'C:\Edgar\Data\IOIResults\j04\ROItest';
pathNameParentList{4} = 'C:\Edgar\Data\IOIResults\j05s01\ROItest';
pathNameParentList{5} = 'C:\Edgar\Data\IOIResults\j06\ROItest';
pathNameParentList{6} = 'C:\Edgar\Data\IOIResults\j07\ROItest';

pathNameParentList{7} = 'C:\Edgar\Data\IOIResults\k01\ROItest';
pathNameParentList{8} = 'C:\Edgar\Data\IOIResults\k02\ROItest';
pathNameParentList{9} = 'C:\Edgar\Data\IOIResults\k03\ROItest';
pathNameParentList{10} = 'C:\Edgar\Data\IOIResults\k04\ROItest';
pathNameParentList{11} = 'C:\Edgar\Data\IOIResults\k05\ROItest';
pathNameParentList{12} = 'C:\Edgar\Data\IOIResults\k06s02\ROItest';

pathNameParentList{13} = 'C:\Edgar\Data\IOIResults\L01\ROItest';
pathNameParentList{14} = 'C:\Edgar\Data\IOIResults\L02\ROItest';
pathNameParentList{15} = 'C:\Edgar\Data\IOIResults\L03\ROItest';
pathNameParentList{16} = 'C:\Edgar\Data\IOIResults\L04\ROItest';
pathNameParentList{17} = 'C:\Edgar\Data\IOIResults\L05\ROItest';
pathNameParentList{18} = 'C:\Edgar\Data\IOIResults\L06\ROItest';
pathNameParentList{19} = 'C:\Edgar\Data\IOIResults\L07\ROItest';
pathNameParentList{20} = 'C:\Edgar\Data\IOIResults\L08\ROItest';

%% Folders with physiological data (in the same order as the IOImat list)
pathNamePhysioList{1} = 'C:\Edgar\Data\IOIData20141016\Physio_Monitoring\142813_J02';
pathNamePhysioList{2} = 'C:\Edgar\Data\IOIData20141211\Physio Monitoring\121423';
pathNamePhysioList{3} = 'C:\Edgar\Data\IOIData20141016\Physio_Monitoring\131700_J04';
pathNamePhysioList{4} = 'C:\Edgar\Data\IOIData20141016\Physio_Monitoring\163300_J05_acq1';
pathNamePhysioList{5} = 'C:\Edgar\Data\IOIData20141016\Physio_Monitoring\154029_J06';
pathNamePhysioList{6} = 'C:\Edgar\Data\IOIData20141016\Physio_Monitoring\160704_J07';

pathNamePhysioList{7} = 'C:\Edgar\Data\IOIData20141023\Physio Monitoring\103931_K01';
pathNamePhysioList{8} = 'C:\Edgar\Data\IOIData20141023\Physio Monitoring\111141_K04';
pathNamePhysioList{9} = 'C:\Edgar\Data\IOIData20141023\Physio Monitoring\113104_K03';
pathNamePhysioList{10} = 'C:\Edgar\Data\IOIData20141023\Physio Monitoring\111141_K04';
pathNamePhysioList{11} = 'C:\Edgar\Data\IOIData20141023\Physio Monitoring\114756_K05';
pathNamePhysioList{12} = 'C:\Edgar\Data\IOIData20141023\Physio Monitoring\100756_K06';

pathNamePhysioList{13} = 'C:\Edgar\Data\IOIData20141211\Physio Monitoring\115503';
pathNamePhysioList{14} = 'C:\Edgar\Data\IOIData20141211\Physio Monitoring\124912';
pathNamePhysioList{15} = 'C:\Edgar\Data\IOIData20141211\Physio Monitoring\104650';
pathNamePhysioList{16} = 'C:\Edgar\Data\IOIData20141211\Physio Monitoring\111837';
pathNamePhysioList{17} = 'C:\Edgar\Data\IOIData20141211\Physio Monitoring\113338';
pathNamePhysioList{18} = 'C:\Edgar\Data\IOIData20141211\Physio Monitoring\130833';
pathNamePhysioList{19} = 'C:\Edgar\Data\IOIData20141211\Physio Monitoring\123555';
pathNamePhysioList{20} = 'C:\Edgar\Data\IOIData20141211\Physio Monitoring\121423';

%% Batch processing of data to be exported to analyzer 2. Do not forget to create the header files
nSessions = 3;
% nSessions =  numel(pathNameParentList);
for iSessions = 3:nSessions,
    OK(iSessions) = ioi_export_2analyzer (pathNameParentList{iSessions}, pathNamePhysioList{iSessions});
end

%% Send e-mail
% List of open inputs
nrun = 1; % enter the number of runs here
jobfile = {'C:\Edgar\Dropbox\PostDoc\Newborn\script_batch_analyzer_export_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(0, nrun);
for crun = 1:nrun
end
spm('defaults', 'FMRI');
spm_jobman('serial', jobs, '', inputs{:});
% EOF