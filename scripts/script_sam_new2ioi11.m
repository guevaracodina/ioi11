%% Pre-processing script before fc processing
% Run first script_ioi_realignment, then this script
%% Path (after running spm8)
addpath(genpath('D:\spm8\toolbox\ioi'))

%% dir field
IOI.dir.dir_group_all = 'D:\Edgar\';
IOI.dir.dir_group_raw = 'D:\Edgar\OIS_Data\';
IOI.dir.dir_group_res = 'D:\Edgar\OIS_Results\';
clear subjectList
subjectList{1} = '16_02_25,LP01a';
subjectList{2} = '16_02_25,LP01b';
subjectList{3} = '16_02_25,NC02';
subjectList{4} = '16_02_25,NC03a';
subjectList{5} = '16_02_25,NC03b';
subjectList{6} = '16_02_26,NC05a';
subjectList{7} = '16_02_26,NC05b';
subjectList{8} = '16_02_26,NC06a';
subjectList{9} = '16_02_26,NC06b';
subjectList{10} = '16_07_07,NC07';

%% Call function
subjects2Run = 10;      % List of subject numbers to run
for iSubjects = subjects2Run
    IOI.subj_name = subjectList{iSubjects};
    ioi_sam_new2ioi11(IOI)
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
