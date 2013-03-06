%% Network analysis (first=level)
% Old script superseded by ioi_network_analyses_run
cd('D:\spm8')
addpath('D:\spm8')
spm fmri
spm quit
addpath(genpath('D:\Edgar\conn'))
pathName = 'E:\Edgar\Data\conn\DATA\conn_NYU\results\firstlevel\ANALYSIS_01\';
currentPath = pwd;
cd(pathName);
fileName{1} = fullfile(pathName, 'resultsROI_Condition001.mat');
fileName{2} = fullfile(pathName, 'resultsROI_Condition002.mat');
fileName{3} = fullfile(pathName, 'resultsROI_Condition003.mat');
% Index vector of ROIs
roiIndex = 1:99;
% 1: Efficiency/cost measures; 2: PathDistance/Clustering measures
measType = [1 2];
% 0: uses raw (Fisher-transformed) correlation coefficients; 
% 1: normalizes across subjects (transform to z-scores)
% 2: Cost
normType = [0 1];
% threshold value to compute adjacency matrix
thr = [0.5 1 1.5];

%% Perform connectivity analysis
clc
for iFiles = 1:numel(fileName)
    for iMeas = 1:numel(measType)
        for iNorm = 1:numel(normType)
            for iThr = 1:numel(thr)
                results = conn_network(fileName{iFiles}, roiIndex, measType(iMeas), normType(iNorm), thr(iThr));
                varName = sprintf('results_F%02d_M%02d_N%02d_T%02d', iFiles, measType(iMeas), normType(iNorm), iThr);
                assignin('base', varName, results);
                fclose('all');
            end
        end
    end
end

% Go back to working folder
cd(currentPath);

%% Save results
matFileName = 'results_ROI_all_conditions';
save(fullfile(pathName, matFileName), 'results_*', 'fileName', 'roiIndex', 'measType', 'normType', 'thr');

% EOF
