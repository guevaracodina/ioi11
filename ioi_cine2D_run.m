function out = ioi_cine2D_run(job)
%select a subset of sessions
if isfield(job.session_choice,'select_sessions')
    all_sessions = 0;
    selected_sessions = job.session_choice.select_sessions.selected_sessions;
else
    all_sessions = 1;
end
%LPF
if isfield(job.lpf_choice,'lpf_gauss_On')
    LPF.lpf_gauss_On = 1;
    LPF.fwhm1 = job.lpf_choice.lpf_gauss_On.fwhm1;
else
    LPF.lpf_gauss_On = 0;
end
which_onset_type = job.which_onset_type;
include_OD = job.include_OD;
include_flow = job.include_flow;
include_HbT = job.include_HbT;
normalize_choice = job.normalize_choice;
show_movie = job.show_movie;
dir_cine = 'Cine_movies';
high_limit = job.high_limit/100;
low_limit = job.low_limit/100;
%get size of windows before and after stimulation to average on, in data points
for SubjIdx=1:length(job.IOImat)
    try
        tic
        clear IOI onsets_list M
        
        if ~isfield(IOI,'dev')
            IOI.dev.TR = 0.2;
        end
        window_after = round(job.window_after/IOI.dev.TR);
        window_before = round(job.window_before/IOI.dev.TR);
        
        %Load IOI.mat information
        IOImat = job.IOImat{SubjIdx};               
        [dir_ioimat dummy] = fileparts(job.IOImat{SubjIdx});
        if isfield(job.IOImatCopyChoice,'IOImatCopy')
            newDir = job.IOImatCopyChoice.IOImatCopy.NewIOIdir;
            newDir = fullfile(dir_ioimat,newDir);
            if ~exist(newDir,'dir'),mkdir(newDir); end
            IOImat = fullfile(newDir,'IOI.mat');
        else
            newDir = dir_ioimat;
        end
        try
            load(IOImat);
        catch
            load(job.IOImat{SubjIdx});
        end
        
        if ~isfield(IOI.res,'cineOK') || job.force_redo           
            cineDir = fullfile(newDir,dir_cine);
            if ~exist(cineDir,'dir'), mkdir(cineDir); end
            %get stimulation information - Careful, here onset duration is ignored!
            if include_HbT
                if ~isfield(IOI.color,'HbT')
                    IOI.color.HbT = 'T';
                    IOI.color.eng = [IOI.color.eng IOI.color.HbT];
                end
            end
            %loop over sessions
            for s1=1:length(IOI.sess_res);
                if all_sessions || sum(s1==selected_sessions)
                    onsets_list{s1} = IOI.sess_res{s1}.onsets;
                    %loop over colors
                    for c1=1:length(IOI.color.eng);
                        do_color = 0;
                        if include_flow
                            if isfield(IOI.color,'flow')
                                if IOI.color.eng(c1)==IOI.color.flow
                                    do_color = 1;
                                end
                            end
                        end
                        if include_OD
                            if any(IOI.color.eng(c1)==[IOI.color.green IOI.color.red IOI.color.yellow])
                                do_color = 1;
                            end
                        end
                        if isfield(IOI.color,'HbO')
                            if include_HbT
                                HbColors = [IOI.color.HbO IOI.color.HbR IOI.color.HbT];
                            else
                                HbColors = [IOI.color.HbO IOI.color.HbR];
                            end
                            
                            if any(IOI.color.eng(c1)==HbColors)
                                do_color = 1;
                            end
                        end
                        if do_color
                            m2 = 0;
                            %loop over onset types
                            for m1=1:length(onsets_list{s1})
                                if any(m1==which_onset_type)
                                    m2 = m2+1;
                                    arrays_assigned = 0;
                                    kb = 0; %counter of segments before onsets
                                    ka = 0; %counter of segments after onsets
                                    kb2 = 0; %counter of skipped segments before onsets
                                    ka2 = 0; %counter of skipped segments after onsets
                                    %loop over onsets for that session
                                    U = round(onsets_list{s1}{m1}/IOI.dev.TR); %in data points
                                    for u1=1:length(U)
                                        %Frames to get:
                                        fr_before = U(u1)-window_before:U(u1)-1;
                                        fr_after = U(u1):U(u1)+window_after-1;
                                        %Load images before and after stimulus
                                        Yb = ioi_get_images(IOI,fr_before,c1,s1,dir_ioimat);
                                        Ya = ioi_get_images(IOI,fr_after,c1,s1,dir_ioimat);
                                        if ~isempty(Ya) && ~isempty(Yb)
                                            if ~arrays_assigned %assign arrays
                                                tmp_array_before = zeros(size(Yb));
                                                tmp_array_after = zeros(size(Ya));
                                                arrays_assigned = 1;
                                            end
                                            
                                            switch normalize_choice
                                                case 1
                                                    tmp_median = median(Yb,3);
                                                case 2
                                                    tmp_median = squeeze(Yb(:,:,end));
                                            end
                                            tmp_median_b = repmat(tmp_median,[1 1 size(Yb,3)]);
                                            tmp_median_a = repmat(tmp_median,[1 1 size(Ya,3)]);
                                            tmp_array_before = tmp_array_before + Yb - tmp_median_b;
                                            tmp_array_after = tmp_array_after + Ya - tmp_median_a;
                                            ka = ka+1;
                                            kb = kb+1;
                                        else
                                            if isempty(Ya)
                                                ka2 = ka2+1;
                                                if ka2<3
                                                    disp(['Could not include segment due to incomplete data after onset ' int2str(u1) ' at time ' num2str(U(u1)*IOI.dev.TR) ...
                                                        ' for session ' int2str(s1) ...
                                                        ' for color ' int2str(c1) ' for onset type ' int2str(m1) ...
                                                        '... skipping ' int2str(ka2) ' so far']);
                                                    
                                                end
                                            end
                                            if isempty(Yb)
                                                kb2 = kb2+1;
                                                if kb2 < 3
                                                    disp(['Could not include segment due to incomplete data before onset ' int2str(u1) ' at time ' num2str(U(u1)*IOI.dev.TR) ...
                                                        ' for session ' int2str(s1) ...
                                                        ' for color ' int2str(c1) ' for onset type ' int2str(m1) ...
                                                        '... skipping ' int2str(kb2) ' so far']);
                                                end
                                            end
                                        end
                                    end
                                    %divide by number of found segments (ka = kb always in this module)
                                end
                                %save movie
                                d = tmp_array_after/ka;
                                m0 = min(d(:));
                                M0 = max(d(:));
                                %best is to save d
                                save(fullfile(cineDir,['cine_S' gen_num_str(s1,2) '_' IOI.color.eng(c1) '_onset' int2str(m1) '.Mdat']),'d');
                                dM = M0-m0;
                                Nf = size(d,3);
                                F(Nf) = struct('cdata',[],'colormap',[]);
                                h0 = figure;
                                %if IOI.color.eng(c1) == IOI.color.HbR
                                %    clims = [m0+dM*(1-5*high_limit) M0+dM*(1-low_limit)];
                                %else
                                clims = [m0+dM*low_limit M0+dM*high_limit];
                                %end
                                try
                                try
                                    fname_movie = fullfile(cineDir,['cine_S' gen_num_str(s1,2) '_' IOI.color.eng(c1) '_onset' int2str(m1) '.avi']);
                                    vidObj = VideoWriter(fname_movie);
                                    open(vidObj);
                                    VideoOK = 1;
                                catch
                                    VideoOK = 0;
                                    fname_movie = fullfile(cineDir,['cine_S' gen_num_str(s1,2) '_' IOI.color.eng(c1) '_onset' int2str(m1) '.mat']);
                                end
                                %set(gca,'NextPlot','replacechildren');
                                for i0=1:size(d,3); %20
                                    try
                                        imagesc(squeeze(d(:,:,i0)),clims);
                                    catch %in case there are Inf in d
                                        imagesc(squeeze(d(:,:,i0)));
                                    end
                                    F(i0) = getframe;
                                    if VideoOK
                                        writeVideo(vidObj,F(i0));
                                    end
                                end
                                if VideoOK
                                    close(vidObj);
                                else
                                    save(fname_movie,'F');
                                end
                                %movie(h0,F,1,1/IOI.dev.TR); %,[0 0 0 0]);
                                try close(h0); end
                                IOI.cine{s1,m1}{c1}.fname_movie = fname_movie;
                                end
                            end
                            if (ka2>0 || kb2 > 0)
                                disp(['Skipped ' int2str(ka2+kb2) ' segments']);
                            end
                        end
                    end
                end
                
            end
            IOI.res.cineOK = 1;
        end
        toc
        disp(['Subject ' int2str(SubjIdx) ' complete']);
        out.IOImat{SubjIdx} = IOImat;
    catch exception
        disp(exception.identifier)
        disp(exception.stack(1))
        out.IOImat{SubjIdx} = job.IOImat{SubjIdx};
    end
end
end