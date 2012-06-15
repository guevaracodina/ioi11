function out = ioi_seeds_run(job)
% Lets the user choose seeds on the anatomical image either manually or in
% automatic fashion.
%_______________________________________________________________________
% Copyright (C) 2010 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%______________________________________________________________________

if isfield(job.AutoSeedChoice,'AutoSeed')
    seedSize = job.AutoSeedChoice.AutoSeed.ArraySeed;
    AutoSeed = 1;
else
    AutoSeed = 0;
    if isfield(job.AutoSeedChoice,'ManualSeed')
        graphicalSeed = true;
    else
        graphicalSeed = false;
    end
end
for SubjIdx=1:length(job.IOImat)
    try
        clear IOI
        %Load IOI.mat information
        IOImat = job.IOImat{SubjIdx};
        load(IOImat);   % load missing? EGC
        [dir_ioimat ~] = fileparts(job.IOImat{SubjIdx});
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
        catch exception % returns the contents of the MException object
            disp(exception.identifier)
            disp(exception.stack(1))
            load(job.IOImat{SubjIdx});
        end
        if ~isfield(IOI.fcIOS,'seedsOK') || job.force_redo
            if job.RemovePreviousSeed
                try
                    IOI.fcIOS = rmfield(IOI.fcIOS,'seedsOK');
                end
                try
                    for i1=1:length(IOI.fcIOS.seeds)
                        %clean up: delete seeds mask files
                        delete(IOI.fcIOS.seeds{i1}.fname);
                    end
                end
                try
                    IOI.fcIOS = rmfield(IOI.fcIOS,'seeds');
                end
                index = 0;
            else
                if isfield(IOI.fcIOS,'seeds')
                    index = length(IOI.fcIOS.seeds);
                else
                    index = 0;
                end
            end
            %display anatomical image
            vol = spm_vol(IOI.res.file_anat);
            vx = [1 1 1];
            [dir1 fil1] = fileparts(vol.fname);
            im_anat = spm_read_vols(vol);
            if ~AutoSeed
                %Display images of changes from 10th to 90th percentile for all sessions
                try
                    for i0=1:length(IOI.sess_res)
                        hs{i0} = figure;
                        V = spm_vol(IOI.sess_res{i0}.fname_change_90_10{1}); %color green
                        tmp_image = spm_read_vols(V);
                        imagesc(tmp_image);
                        title(['Session ' int2str(i0) ': ratio of 90th to 10th percentile']);
                    end
                end
                oneMoreSeed = true;
                % h2 = figure; 
                h2 = spm_figure('GetWin', 'Interactive');
                spm_figure('Clear', 'Interactive');
                spm_input(['Subject ' int2str(SubjIdx) ' (' IOI.subj_name ')' ],'-1','d');
                full_mask = ones(size(im_anat));
                linecount = 0;
                while oneMoreSeed
                    figure(h2);
                    oneMoreSeed = spm_input('Add another seed?',2*linecount+2,'y/n');
                    if oneMoreSeed == 'y', oneMoreSeed = true; else oneMoreSeed = false; end
                    if oneMoreSeed
%                         h1 = figure('Position',[20 50 3*size(im_anat,1) 3*size(im_anat,2)]);
                        h1 = spm_figure('GetWin', 'Graphics');
                        spm_figure('Clear', 'Graphics');
                        cmap = contrast(im_anat .* full_mask);
                        imagesc(im_anat .* full_mask);
                        colormap(cmap); 
                        if graphicalSeed
                            title('Make ROI polygon, then double click in it to create seed.');
                            mask = roipoly;
                        else
                            linecount = linecount + 1;
                            rc = spm_input('Enter [row,column] of center',2*linecount+1,'e',[],2);
                            linecount = linecount + 1;
                            radius = spm_input('Enter radius in pixels',2*linecount+1,'e',[],1);
                            radius = round(radius);
                            if radius < 0, radius = 0; end
                            mask = zeros(size(im_anat));
                            for x1=-radius:radius
                                for y1=-radius:radius
                                    if x1^2+y1^2 <= radius^2
                                        try %will skip pixels outside the image
                                            mask(rc(1)+y1,rc(2)+x1) = 1;
                                        end
                                    end
                                end
                            end
                        end
                        mask = single(mask);
                        full_mask = full_mask-mask;
                        index = index + 1; linecount = linecount + 1;
                        %imagesc(mask);
                        if job.select_names
                            figure(h2);
                            name = spm_input(['Enter name of seed' int2str(index)],2*linecount+1,'s');
                        else
                            name = int2str(index);
                        end
                        IOI.fcIOS.seeds{index}.name = name;
                        if index < 10, str0 = '0'; else str0 = ''; end
                        str = [str0 int2str(index)];
                        fname_mask = fullfile(dir1,[fil1 '_seed_' str '.nii']);
                        IOI.fcIOS.seeds{index}.fname = fname_mask;
                        ioi_save_nifti(mask, fname_mask, vx);
                    end
                end
%                 try close(h1); end
%                 try close(h2); end
                for i0=1:length(IOI.sess_res)
                    try close(hs{i0}); end
                end
            else
                %automatic seeds
                %h1 = figure('Position',[20 50 3*size(im_anat,1) 3*size(im_anat,2)]);
                sz = size(im_anat);
                sz = sz(1:2); %remove 3rd component, which is one
                sz = floor(sz./seedSize);
                for i1 = 1:seedSize(1) %N
                    for i2 = 1:seedSize(2) %M
                        %ROI created in following order
                        % 1      2       3 ...   seedSize(2)=M
                        % M+1   M+2     M+3 ...  2*M
                        %  .    .        .       .
                        %  .    .        .       N*M
                        index = index+1;
                        %vertices in the order: UL, UR, LR, LL
                        mask = roipoly(im_anat,...
                            [1+sz(2)*(i2-1) 1+sz(2)*i2 1+sz(2)*i2 1+sz(2)*(i2-1)],...
                            [1+sz(1)*(i1-1) 1+sz(1)*(i1-1) 1+sz(1)*i1 1+sz(1)*i1]);
                        mask = single(mask);
                        if index < 10, str0 = '0'; else str0 = ''; end
                        str = [str0 int2str(index)];
                        fname_mask = fullfile(dir1,[fil1 '_seed_' str '_' int2str(i1) 'x' int2str(i2) '.nii']);
                        IOI.fcIOS.seeds{index}.fname = fname_mask;
                        ioi_save_nifti(mask, fname_mask, vx);
                    end
                end
            end
            IOI.seedName = {};
            for i0=1:length(IOI.fcIOS.seeds)
                if isfield(IOI.fcIOS.seeds{i0},'name')
                    IOI.seedName = [IOI.seedName; IOI.fcIOS.seeds{i0}.name];
                else
                    IOI.seedName = [IOI.seedName; ['seed' gen_num_str(i0,3)]];
                end
            end
            IOI.fcIOS.seedsOK = true;
            if isfield(job.IOImatCopyChoice,'IOImatCopy')
                newDir = job.IOImatCopyChoice.IOImatCopy.NewIOIdir;
                newDir = fullfile(dir_ioimat,newDir);
                if ~exist(newDir,'dir'),mkdir(newDir); end
                IOImat = fullfile(newDir,'IOI.mat');
            end
            save(IOImat,'IOI');
        end
         
        disp(['Subject ' int2str(SubjIdx) ' (' IOI.subj_name ')' ' complete']);
        out.IOImat{SubjIdx} = IOImat;
    catch exception
        disp(exception.identifier)
        disp(exception.stack(1))
        out.IOImat{SubjIdx} = job.IOImat{SubjIdx};
    end
end
