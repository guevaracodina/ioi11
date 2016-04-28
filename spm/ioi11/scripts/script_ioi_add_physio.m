%% script_ioi_add_physio
clear; close all; clc
topDir{1} = 'C:\Edgar\Data\IOIResults\j02\';
topDir{2} = 'C:\Edgar\Data\IOIResults\j03\';
topDir{3} = 'C:\Edgar\Data\IOIResults\j04\';
% topDir{4} = 'C:\Edgar\Data\IOIResults\j05s01\';
topDir{5} = 'C:\Edgar\Data\IOIResults\j06\';
% topDir{6} = 'C:\Edgar\Data\IOIResults\j07\';
topDir{7} = 'C:\Edgar\Data\IOIResults\k01\';
topDir{8} = 'C:\Edgar\Data\IOIResults\k02\';
topDir{9} = 'C:\Edgar\Data\IOIResults\k03\';
topDir{10} = 'C:\Edgar\Data\IOIResults\k04\';
topDir{11} = 'C:\Edgar\Data\IOIResults\k05\';
topDir{12} = 'C:\Edgar\Data\IOIResults\k06s01\';
topDir{13} = 'C:\Edgar\Data\IOIResults\L01\';
topDir{14} = 'C:\Edgar\Data\IOIResults\L02\';
topDir{15} = 'C:\Edgar\Data\IOIResults\L03\';
topDir{16} = 'C:\Edgar\Data\IOIResults\L04\';
topDir{17} = 'C:\Edgar\Data\IOIResults\L05\';
topDir{18} = 'C:\Edgar\Data\IOIResults\L06\';
topDir{19} = 'C:\Edgar\Data\IOIResults\L07\';
topDir{20} = 'C:\Edgar\Data\IOIResults\L08\';

currentDir = 'ROItest';

%%
for iDirs = 5:numel(topDir)
    % Choose heart rate data
    [dirPhysio, sts] = cfg_getfile([1 1],'dir','Select folder',{fullfile(topDir{iDirs},currentDir)}, topDir, '.*');
    dirListIOI = dir(fullfile(topDir{iDirs},[currentDir filesep 'IOI.mat']));
    [heartRateFile, sts] = cfg_getfile(1,'mat','Select heart rate file',[], dirPhysio{1}, 'heartRate*');
    
    % Backup IOI
    IOImat = fullfile(topDir{iDirs},[currentDir filesep dirListIOI.name]);
    copyfile(IOImat,...
        fullfile(topDir{iDirs},[currentDir filesep 'IOI.bak']));
    
    % Point to the heart rate file and save IOI mat
    load(IOImat);
    IOI.fcIOS.SPM(1).physioHeartFile = heartRateFile;
    save(IOImat, 'IOI')
end