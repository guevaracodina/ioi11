function out = ioi_extract_roi_time_series_run(job)
% Extract the time series for the ROIs (seeds), loops along subjects, sessions,
% colors and files. Needs ioi_extract_core
% Added brain mask extraction if needed
%_______________________________________________________________________________
% Copyright (C) 2012 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________

for SubjIdx=1:length(job.IOImat)
    try
        tic
        clear IOI ROI
        %Load IOI.mat information
        IOImat = job.IOImat{SubjIdx};
        
        % IOI copy/overwrite
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
        
        if ~isfield(IOI.res,'ROIOK')
            disp(['No ROI available for subject ' int2str(SubjIdx) ' ... skipping series extraction']);
        else
            if ~isfield(IOI.res,'seriesOK') || job.force_redo
                % Get mask for each ROI
                [IOI mask] = ioi_get_ROImask(IOI,job);
                % Extract ROI
                [ROI IOI] = ioi_extract_core(IOI,job,mask);
                IOI.res.seriesOK = 1;
                if isfield(job.IOImatCopyChoice,'IOImatCopy')
                    newDir = job.IOImatCopyChoice.IOImatCopy.NewIOIdir;
                    newDir = fullfile(dir_ioimat,newDir);
                    if ~exist(newDir,'dir'),mkdir(newDir); end
                    IOImat = fullfile(newDir,'IOI.mat');
                end
                [dir1 dummy] = fileparts(IOImat);
                ROIfname = fullfile(dir1,'ROI.mat');
                save(ROIfname,'ROI');
                IOI.ROI.ROIfname = ROIfname;
                save(IOImat,'IOI');
            end
        end
        if job.extractBrainMask
            if ~isfield(IOI.fcIOS.mask, 'maskOK')
                disp(['No brain mask available for subject ' int2str(SubjIdx) ' (' IOI.subj_name ')'  ' ... skipping series extraction']);
            elseif job.extractBrainMask
                % ------------------ Extract brain mask here -----------------------
                if ~isfield(IOI.fcIOS.mask,'seriesOK') || job.force_redo
                    disp(['Extracting time series of global brain signal for subject ' int2str(SubjIdx) ' (' IOI.subj_name ')']);
                    % It means we will extract only 1 ROI
                    job.extractingBrainMask = true;
                    % Get brain mask
                    [IOI mask] = ioi_get_brain_mask(IOI,job);
                    % Extract brain mask here
                    [brainMaskSeries IOI] = ioi_extract_core(IOI,job,mask);
                    % Reset flag
                    job.extractingBrainMask = false;
                    % Brain mask time series extraction succesful!
                    IOI.fcIOS.mask.seriesOK = true;
                    if isfield(job.IOImatCopyChoice,'IOImatCopy')
                        newDir = job.IOImatCopyChoice.IOImatCopy.NewIOIdir;
                        newDir = fullfile(dir_ioimat,newDir);
                        if ~exist(newDir,'dir'),mkdir(newDir); end
                        IOImat = fullfile(newDir,'IOI.mat');
                    end
                    [dir1 dummy] = fileparts(IOImat);
                    fnameSeries = fullfile(dir1,'brainMaskSeries.mat');
                    save(fnameSeries,'brainMaskSeries');
                    % identify in IOI the file name of the time series
                    IOI.fcIOS.mask.fnameSeries = fnameSeries;
                    save(IOImat,'IOI');
                end
            end
        end
        disp(['Elapsed time: ' datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS')]);
        disp(['Subject ' int2str(SubjIdx) ' complete']);
        out.IOImat{SubjIdx} = IOImat;
    catch exception
        disp(exception.identifier)
        disp(exception.stack(1))
        out.IOImat{SubjIdx} = job.IOImat{SubjIdx};
    end
end
