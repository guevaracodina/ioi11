function out = ioi_extract_roi_time_series_run(job)
% Extract the time series for the ROIs (seeds), loops along subjects, sessions,
% colors and files. Needs ioi_extract_core
% Added brain mask extraction if needed
%_______________________________________________________________________________
% Copyright (C) 2012 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________

if isfield(job,'save_figures')
    %save_figures
    save_figures = job.save_figures;
    generate_figures = job.generate_figures;
else
    save_figures = 0;
    generate_figures = 0;
end
for SubjIdx=1:length(job.IOImat)
    try
        tic
        clear ROI
        %Load IOI.mat information
        [IOI IOImat dir_ioimat]= ioi_get_IOI(job,SubjIdx);
        
        if ~isfield(IOI.res,'ROIOK')
            disp(['No ROI available for subject ' int2str(SubjIdx) ' ... skipping series extraction']);
        else
            if ~isfield(IOI.res,'seriesOK') || job.force_redo
                % Get mask for each ROI
                [IOI mask] = ioi_get_ROImask(IOI,job);
                Amask = []; %Initialize activation mask                
                if isfield(job.activMask_choice,'activMask')
                    try
                        mask_image = job.activMask_choice.activMask.mask_image{1};
                        threshold = job.activMask_choice.activMask.threshold;
                        two_sided = job.activMask_choice.activMask.two_sided;
                        
                        h=hgload(mask_image);
                        ch=get(h,'Children');
                        l=get(ch,'Children');
                        z=get(l{3},'cdata');
                        if two_sided
                            Amask = z > abs(threshold) | z < -abs(threshold);
                        else
                            if z > 0
                                Amask = z > threshold;
                            else
                                Amask =  z < threshold;
                            end
                        end
                        close(h);
                        clear z l ch threshold two_sided mask_image
                    catch
                        disp('Could not mask by specified mask -- no masking by activation will be done')
                        Amask = [];
                    end
                end
                % We are not extracting brain mask here
                job.extractingBrainMask = false;
                % Extract ROI
                [ROI IOI] = ioi_extract_core(IOI,job,mask,Amask);
                IOI.res.seriesOK = 1;
                
                ROIfname = fullfile(dir_ioimat,'ROI.mat');
                save(ROIfname,'ROI');
                IOI.ROI.ROIfname = ROIfname;
                save(IOImat,'IOI');
            end
        end
        if job.extractBrainMask
            if ~isfield(IOI.fcIOS.mask, 'maskOK')
                disp(['No brain mask available for subject ' int2str(SubjIdx) ...
                    ' (' IOI.subj_name ')'  ' ... skipping series extraction']);
            elseif job.extractBrainMask
                % ------------------ Extract brain mask here -----------------------
                if ~isfield(IOI.fcIOS.mask,'seriesOK') || job.force_redo
                    disp(['Extracting time series of global brain signal for subject ' ...
                        int2str(SubjIdx) ' (' IOI.subj_name ')']);
                    % It means we will extract only 1 ROI
                    job.extractingBrainMask = true;
                    % Get brain mask
                    [IOI mask] = ioi_get_brain_mask(IOI,job);
                    % Extract brain mask here
                    [brainMaskSeries IOI] = ioi_extract_core(IOI,job,mask,Amask);
                    % Reset flag
                    job.extractingBrainMask = false;
                    % Brain mask time series extraction succesful!
                    IOI.fcIOS.mask.seriesOK = true;
                    
                    fnameSeries = fullfile(dir_ioimat,'brainMaskSeries.mat');
                    save(fnameSeries,'brainMaskSeries');
                    % identify in IOI the file name of the time series
                    IOI.fcIOS.mask.fnameSeries = fnameSeries;
                    save(IOImat,'IOI');
                end
            end
        end
        if generate_figures || save_figures
            ioi_plot_roi_time_series(IOI,ROI,generate_figures,save_figures);
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
