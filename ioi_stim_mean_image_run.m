function out = ioi_stim_mean_image_run(job)
%select a subset of sessions
[all_sessions selected_sessions] = ioi_get_sessions(job);
%filters
HPF = ioi_get_HPF(job);
LPF = ioi_get_LPF(job);

%Spatial filter
if isfield(job.spatial_LPF,'spatial_LPF_On')
    radius = job.spatial_LPF.spatial_LPF_On.spatial_LPF_radius;
    spatial_LPF = 1;
else
    spatial_LPF = 0;
end
%shrinking of images
[shrinkage_choice SH] = ioi_get_shrinkage_choice(job);
[Z.rmi Z.ust Z.remove_stims_SD] = ioi_get_stims_choices(job);
Z.remove_segment_drift = job.remove_segment_drift;
%Other options
IC = job.IC; %colors to include
%save_figures
Z.save_figures = job.save_figures;
Z.generate_figures = job.generate_figures;
Z.normalize_choice = job.normalize_choice;
Z.interactive_mode = job.interactive_mode;
for SubjIdx=1:length(job.IOImat)
    try
        clear IOI ROI onsets_list M
        %Load IOI.mat information
        IOImat = job.IOImat{SubjIdx};
        load(IOImat);
        %check whether to skip any calculations, and create new IOI directory
        if ~isfield(IOI.res,'meanimageOK') || job.force_redo
            [dir_ioimat dummy] = fileparts(job.IOImat{SubjIdx});
            if isfield(job.IOImatCopyChoice,'IOImatCopy')
                newDir = job.IOImatCopyChoice.IOImatCopy.NewIOIdir;
                newDir = fullfile(dir_ioimat,newDir);
                if ~exist(newDir,'dir'),mkdir(newDir); end
                IOImat = fullfile(newDir,'IOI.mat');
            else
                newDir = dir_ioimat;
            end
            %create directory for figures
            if Z.save_figures
                dir_fig = fullfile(newDir,'fig');
                if ~exist(dir_fig,'dir'),mkdir(dir_fig);end
            end
            if ~isfield(IOI,'dev')
                IOI.dev.TR = 0.2;
                disp('Careful, assuming TR = 0.2 s')
            end
            %rescale windows to data points (round off to integers)
            Z.window_after = round(job.window_after/IOI.dev.TR);
            Z.window_before = round(job.window_before/IOI.dev.TR);
            Z.window_offset = round(job.window_offset/IOI.dev.TR);
            Z.window_start_delay = round(job.window_start_delay/IOI.dev.TR);
            %Include HbT
            if IC.include_HbT
                if ~isfield(IOI.color,'HbT')
                    IOI.color.HbT = 'T';
                    IOI.color.eng = [IOI.color.eng IOI.color.HbT];
                end
            end
            Nc = length(IOI.color.eng);
            %save shrunk images
            if shrinkage_choice
                IOI = ioi_save_shrunk_images(IOI,job,SH,dir_ioimat);
            end
            %restric onsets
            [IOI onsets_list] = ioi_restrict_onsets(IOI,job,Z.rmi,Z.ust);
            %Loop over sessions
            for s1=1:length(IOI.sess_res)
                if all_sessions || sum(s1==selected_sessions)
                    %loop over colors
                    for c1=1:Nc
                        %check whether this color should be done or not
                        doColor = ioi_doColor(IOI,c1,IC);
                        if doColor
                            %losad image data
                            y = ioi_get_images(IOI,1:IOI.sess_res{s1}.n_frames,c1,s1,dir_ioimat,shrinkage_choice);
                            [nx ny nt] = size(y);
                            if ~isempty(y)
                                y = reshape(permute(y,[3 1 2]),nt,nx*ny);
                                if HPF.hpf_butter_On
                                    y = ButterHPF(1/IOI.dev.TR,HPF.hpf_butter_freq,HPF.hpf_butter_order,y);
                                end
                                if LPF.lpf_gauss_On
                                    K = get_K(1:nt,LPF.fwhm1,IOI.dev.TR);
                                    y = ioi_filter_HPF_LPF_WMDL(K,y);
                                end
                                
                                %
                                %                                              %normalize flow
                                %                                              if isfield(IOI.color,'flow')
                                %                                                  if IOI.color.eng(c1)==IOI.color.flow
                                %                                                      tmp_d = tmp_d/mean(tmp_d); %or median
                                %                                                  end
                                %                                              end
                                
                                %                                             %remove jumps options:
                                %                                             OP.Sb = 4; %number of standard deviations
                                %                                             OP.Nr = 1/IOI.dev.TR; %number of points removed before and after
                                %                                             OP.Mp = 10/IOI.dev.TR; %size of gaps to be filled
                                %                                             OP.sf = 1/IOI.dev.TR; %sampling frequency
                                %                                             OP.ubf = 1; %use Butterworth filter
                                %                                             OP.bf = 0.01; %Butterworth HPF cutoff
                                %                                             OP.bo = 2; %Butterworth order
                                %
                                %                                             %if ~isempty(rmi{i0})
                                %                                             if IOI.color.eng(c1) == IOI.color.flow
                                %                                                 OP.ubf = 0;
                                %                                             else
                                %                                                 OP.ubf = 1;
                                %                                             end
                                %                                             tmp_d = ioi_remove_jumps(tmp_d,OP);
                                %end
                                
                                %y = permute(reshape(y,nt,nx,ny),[2 3 1]);
                                
                            end
                            %loop over onset types
                            for m1=1:length(onsets_list{s1})
                                if isempty(job.which_onset_type) || any(m1==job.which_onset_type)
                                    %fill structure to pass
                                    Z.ons = onsets_list{s1}{m1};
                                    Z.s1 = s1; Z.m1 = m1; Z.c1 = c1;
                                    [IOI,D,Z] = ioi_average_image_core(IOI,y,Z);
                                    Ma = squeeze(reshape(D.Ma,[nx ny]));
                                    Da = squeeze(reshape(D.Da,[nx ny]));
                                    Ta = squeeze(reshape(D.Ma./(D.Da/sqrt(D.ka)),[nx ny]));
                                    IOI.Avg.ka{s1}{c1,m1} = D.ka;
                                    IOI.Avg.kb{s1}{c1,m1} = D.kb;
                                    IOI.Avg.ka2{s1}{c1,m1} = D.ka2;
                                    IOI.Avg.kb2{s1}{c1,m1} = D.kb2;
                                    fname = fullfile(dir_fig,['Avg_S' int2str(s1) IOI.color.eng(c1) int2str(m1)]);
                                    IOI.Avg.fname{s1}{c1,m1} = fname;                                    
                                    tit = ['Average, Session' int2str(s1) ', Color ' IOI.color.eng(c1) ', Stimulus ' int2str(m1)];
                                    ioi_save_images(Ma,fname,[1 1 1],Z,tit);
                                    fname = fullfile(dir_fig,['Std_S' int2str(s1) IOI.color.eng(c1) int2str(m1)]);
                                    IOI.Avg.fname_std{s1}{c1,m1} = fname;                                    
                                    tit = ['Std, Session' int2str(s1) ', Color ' IOI.color.eng(c1) ', Stimulus ' int2str(m1)];
                                    ioi_save_images(Da,fname,[1 1 1],Z,tit);
                                    fname = fullfile(dir_fig,['T_S' int2str(s1) IOI.color.eng(c1) int2str(m1)]);
                                    IOI.Avg.fname_t{s1}{c1,m1} = fname;                                    
                                    tit = ['Tstat, Session' int2str(s1) ', Color ' IOI.color.eng(c1) ', Stimulus ' int2str(m1)];
                                    ioi_save_images(Ta,fname,[1 1 1],Z,tit);
%                                     try
%                                         %Superpose onto anatomical image
%                                         vol = spm_vol(IOI.res.file_anat);
%                                         A.im_anat = spm_read_vols(vol);
%                                         ioi_save_images_anat(Ta,fname,[1 1 1],Z,tit,A);
%                                     end
                                    %Filtered images
                                    if spatial_LPF
                                        Ks.k1 = nx;
                                        Ks.k2 = ny;
                                        Ks.radius = radius;
                                        Ks = ioi_spatial_LPF('set',Ks);
                                        %Gaussian spatial low pass filter
                                        Ma = ioi_spatial_LPF('lpf',Ks,Ma);
                                        Da = ioi_spatial_LPF('lpf',Ks,Da);
                                        Ta = Ma./(Da/sqrt(D.ka));
                                        fname = fullfile(dir_fig,['FiltT_S' int2str(s1) IOI.color.eng(c1) int2str(m1)]);
                                        IOI.Avg.fname_filt_t{s1}{c1,m1} = fname;
                                        tit = ['FiltTstat, Session' int2str(s1) ', Color ' IOI.color.eng(c1) ', Stimulus ' int2str(m1)];
                                        ioi_save_images(Ta,fname,[1 1 1],Z,tit);
                                    end
                                    
                                end
                            end
                        end
                    end
                end
            end
            
            IOI.res.meanimageOK = 1;            
            save(IOImat,'IOI');
        end
        
        disp(['Subject ' int2str(SubjIdx) ' complete']);
        out.IOImat{SubjIdx} = IOImat;
    catch exception
        disp(exception.identifier)
        disp(exception.stack(1))
        out.IOImat{SubjIdx} = job.IOImat{SubjIdx};
    end
end