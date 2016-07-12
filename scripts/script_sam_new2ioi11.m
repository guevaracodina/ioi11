%% Pre-processing script before fc processing
% Run first script_ioi_realignment.m, then this script
moveData = true;        % Flag to move data to C:\ 
%% Path (after running spm8)
addpath(genpath('D:\spm8\toolbox\ioi'))

%% dir field
IOI.dir.dir_group_all = 'D:\Edgar\';
IOI.dir.dir_group_raw = 'D:\Edgar\OIS_Data\';
IOI.dir.dir_group_res = 'D:\Edgar\OIS_Results\';

%% Subject list for processing
clear subjectList
subjectList{1}  = '16_02_25,LP01a';
subjectList{2}  = '16_02_25,LP01b';
subjectList{3}  = '16_02_25,NC02';
subjectList{4}  = '16_02_25,NC03a';
subjectList{5}  = '16_02_25,NC03b';
subjectList{6}  = '16_02_26,NC05a';
subjectList{7}  = '16_02_26,NC05b';
subjectList{8}  = '16_02_26,NC06a';
subjectList{9}  = '16_02_26,NC06b';
subjectList{10} = '16_07_07,NC07';
subjectList{11} = '16_07_07,LP02a';
subjectList{12} = '16_07_07,LP02b';
subjectList{13} = '16_07_08,LP03';
subjectList{14} = '16_07_08,LP04';
subjectList{15} = '16_07_08,LP05a';
subjectList{16} = '16_07_08,LP05b';
subjectList{17} = '16_07_08,NC08';
subjectList{18} = '16_07_08,NC09';

%% Call function
subjects2Run = 11:18;      % List of subject numbers to run
for iSubjects = subjects2Run
    IOI.subj_name = subjectList{iSubjects};
    ioi_sam_new2ioi11(IOI)
    if moveData
        movefile(fullfile(IOI.dir.dir_group_raw,subjectList{iSubjects}),...
            fullfile('C:\Edgar\OIS_Data\',subjectList{iSubjects}))
        disp(strcat('Current folder moved to:',fullfile('C:\Edgar\OIS_Data\',subjectList{iSubjects})))
    end
end

%% Send e-mail
% List of open inputs
nrun = 1; % enter the number of runs here
jobfile = {'D:\Edgar\send_email_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(0, nrun);
for crun = 1:nrun
end
spm('defaults', 'FMRI');
spm_jobman('serial', jobs, '', inputs{:});
