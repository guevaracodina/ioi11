%% Load IOI mat
addpath(genpath('D:\spm8\toolbox'))
load('D:\Edgar\OIS_Results\12_10_18,NC09\IOI.mat')

%% List names
IOInames = fieldnamesr(IOI);
IOI.warning = {};
IOI.subj_name = '12_10_18,NC09';
IOI.subj_OK = 1;

%% dir field

IOI.dir.dir_group_all = 'D:\Edgar\';
IOI.dir.dir_group_raw = 'D:\Edgar\OIS_Data\';
IOI.dir.dir_group_res = 'D:\Edgar\OIS_Results\';
IOI.dir.dir_subj_raw = 'D:\Edgar\OIS_Data\12_10_18,NC09\';
IOI.dir.dir_subj_res = 'D:\Edgar\OIS_Results\12_10_18,NC09';

%% Color field
IOI.color.eng = 'RGYLOD';
IOI.color.red = 'R';
IOI.color.green = 'G';
IOI.color.yellow = 'Y';
IOI.color.laser = 'L';
IOI.color.HbO = 'O';
IOI.color.HbR = 'D';

%% dev field
IOI.dev.TR = 0.2;

%% conc field
IOI.conc.baseline_hbt = 100;
IOI.conc.baseline_hbo = 60;
IOI.conc.baseline_hbr = 40;

%% ROIname field
IOI.ROIname = {'M_L'; 'M_R'};

%% ROI field
IOI.ROI.ROIfname = 'D:\Edgar\OIS_Results\12_10_18,NC09\ROI.mat';

%% sess_res field
IOI.sess_res{1, 1}.names = 'Stim_1';
IOI.sess_res{1, 1}.onsets = {[]};
IOI.sess_res{1, 1}.durations{1, 1} = [];
IOI.sess_res{1, 1}.parameters = {[]};
IOI.sess_res{1, 1}.n_frames = 4313;
IOI.sess_res{1, 1}.fname{1}
IOI.sess_res{1, 1}.si
IOI.sess_res{1, 1}.ei
IOI.sess_res{1, 1}.fname_median
IOI.sess_res{1, 1}.hasRGY
IOI.sess_res{1, 1}.availCol
IOI.sess_res{1, 1}.missingFrames



    {1x34 cell}    {1x34 cell}    {1x34 cell}    {1x34 cell}    {34x1 cell}    {34x1 cell}


ans = 

  Columns 1 through 11

    [1]    [129]    [257]    [385]    [513]    [641]    [769]    [897]    [1025]    [1153]    [1281]

  Columns 12 through 21

    [1409]    [1537]    [1665]    [1793]    [1921]    [2049]    [2177]    [2305]    [2433]    [2561]

  Columns 22 through 31

    [2689]    [2817]    [2945]    [3073]    [3201]    [3329]    [3457]    [3585]    [3713]    [3841]

  Columns 32 through 34

    [3969]    [4097]    [4225]


ans = 

  Columns 1 through 11

    [128]    [256]    [384]    [512]    [640]    [768]    [896]    [1024]    [1152]    [1280]    [1408]

  Columns 12 through 21

    [1536]    [1664]    [1792]    [1920]    [2048]    [2176]    [2304]    [2432]    [2560]    [2688]

  Columns 22 through 31

    [2816]    [2944]    [3072]    [3200]    [3328]    [3456]    [3584]    [3712]    [3840]    [3968]

  Columns 32 through 34

    [4096]    [4224]    [4313]


ans = 

    [1x68 char]    [1x68 char]    [1x68 char]


ans =

RGY


ans =

RGYL


ans =

   854


%% res fied

%% fcIOS field

















