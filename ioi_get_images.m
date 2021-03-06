function Ya = ioi_get_images(IOI,frames,c1,s1,dir_ioimat,shrinkage_choice)
%IOI: the IOI structure
%frames: the frames to be loaded
%c1: the color desired
%s1: the session number
%dir_ioimat: directory for the IOI.mat of interest
%shrinkage_choice: whether the images should be shrunk for faster speed
try
    Ya = [];
    doHbT = 0;
    try
        if IOI.color.eng(c1) == IOI.color.HbT
            doHbT = 1;
            [cHbR cHbO] = ioi_find_HbRHbO(IOI,s1);
        end
    end
    first_pass = 1;
    fmod_previous = 1;
    frames_loaded = 0;
    image_block = 1;
    if ~any(frames < 1)
        for i=1:length(frames)
            %Loop through the values of ei to find the right block of images
            j0 = image_block;
            try %new format
                while j0<length(IOI.sess_res{s1}.ei)
                    if frames(i) <= IOI.sess_res{s1}.ei{j0}
                        %this is the right block
                        break
                    else
                        j0 = j0+1;
                    end
                end
                if j0 > image_block
                    frames_loaded = 0;
                end
                image_block = j0;
                if image_block>1
                    fmod = frames(i)-IOI.sess_res{s1}.ei{image_block-1};
                else
                    fmod = frames(i);
                end
                fct = image_block;
            catch %old format
                if i==1
                    %find number of frames per block %careful, this can vary
                    %if length(IOI.sess_res{s1}.si) > 1
                    % fpb = 500; %IOI.sess_res{s1}.si{2}-IOI.sess_res{s1}.si{1};
                    fpb = 2000; % New Sam's format (2016)
                    %else
                    %    fpb = IOI.sess_res{s1}.ei{1}-IOI.sess_res{s1}.si{1}+1;
                    %end
                end
                fmod = mod(frames(i),fpb);
                %find location of desired frame in this block
                
                if fmod == 0
                    fmod = fpb;
                end
                if fmod <= fmod_previous
                    frames_loaded = 0;
                end
                fct = ceil(frames(i)/fpb);
            end
            
            if ~frames_loaded
                try
                    if ~isfield(IOI,'sess_shrunk') || ~shrinkage_choice
                        if ~doHbT
                            V = spm_vol(IOI.sess_res{s1}.fname{c1}{fct});
                            Y = spm_read_vols(V);
                            frames_loaded = 1;
                        else
                            V1 = spm_vol(IOI.sess_res{s1}.fname{cHbO}{fct});
                            V2 = spm_vol(IOI.sess_res{s1}.fname{cHbR}{fct});
                            Y = spm_read_vols(V1)+spm_read_vols(V2);
                            frames_loaded = 1;
                        end
                    else
                        V = spm_vol(IOI.sess_shrunk{s1}.fname{c1}{fct});
                        Y = spm_read_vols(V);
                        frames_loaded = 1;
                    end
                catch
                    if ~isfield(IOI,'sess_shrunk') || ~shrinkage_choice
                        if ~doHbT
                            [dir0 fil0 ext0] = fileparts(IOI.sess_res{s1}.fname{c1}{fct});
                            fsep = strfind(dir0,filesep);
                            res = strfind(dir0,'Res');
                            fgr = fsep(fsep > res);
                            tdir = dir0(1:fgr(2));
                            fsep = strfind(dir_ioimat,filesep);
                            res = strfind(dir_ioimat,'Res');
                            fgr = fsep(fsep > res);
                            tdir = dir_ioimat(1:fgr(2));
                            
                            tfname = fullfile(tdir,['S' gen_num_str(s1,2)],[fil0 ext0]);
                            try
                                V = spm_vol(tfname);
                                Y = spm_read_vols(V);
                                frames_loaded = 1;
                            catch
                                %file does not exist
                                return
                            end
                        else
                            [dir0 fil0 ext0] = fileparts(IOI.sess_res{s1}.fname{cHbO}{fct});
                            [dir0 fil2 ext0] = fileparts(IOI.sess_res{s1}.fname{cHbR}{fct});
                            fsep = strfind(dir0,filesep);
                            res = strfind(dir0,'Res');
                            fgr = fsep(fsep > res);
                            tdir = dir0(1:fgr(2));
                            fsep = strfind(dir_ioimat,filesep);
                            res = strfind(dir_ioimat,'Res');
                            fgr = fsep(fsep > res);
                            tdir = dir_ioimat(1:fgr(2));
                            
                            tfname1 = fullfile(tdir,['S' gen_num_str(s1,2)],[fil0 ext0]);
                            tfname2 = fullfile(tdir,['S' gen_num_str(s1,2)],[fil2 ext0]);
                            try
                                V1 = spm_vol(tfname1);
                                V2 = spm_vol(tfname2);
                                Y = spm_read_vols(V1)+spm_read_vols(V2);
                                frames_loaded = 1;
                            catch
                                %file does not exist
                                return
                            end
                        end
                    else
                        [dir0 fil0 ext0] = fileparts(IOI.sess_shrunk{s1}.fname{c1}{fct});
                        fsep = strfind(dir0,filesep);
                        res = strfind(dir0,'Res');
                        fgr = fsep(fsep > res);
                        tdir = dir0(1:fgr(2));
                        
                        fsep = strfind(dir_ioimat,filesep);
                        res = strfind(dir_ioimat,'Res');
                        fgr = fsep(fsep > res);
                        tdir = dir_ioimat(1:fgr(2));
                        tfname = fullfile(tdir,['S' gen_num_str(s1,2)],[fil0 ext0]);
                        try
                            V = spm_vol(tfname);
                            Y = spm_read_vols(V);
                            frames_loaded = 1;
                        catch
                            %file does not exist
                            return
                        end
                    end
                end
                if first_pass
                    %initialize Ya
                    Ya = zeros(size(Y,1),size(Y,2),length(frames));
                    first_pass = 0;
                end
            end
            try
                try
                    Ya(:,:,i) = squeeze(Y(:,:,1,fmod));
                catch
                    Ya(:,:,i) = squeeze(Y(:,:,fmod));
                end
            catch
                %There is one frame missing?
                Ya = [];
                return
            end
        end
    end
catch exception
    disp(exception.identifier)
    disp(exception.stack(1))
end