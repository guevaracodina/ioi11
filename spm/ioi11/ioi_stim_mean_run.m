function out = ioi_stim_mean_run(job)
[all_sessions selected_sessions] = ioi_get_sessions(job);
if isfield(job,'remove_stims')
    rmi = job.remove_stims;
else
    rmi = '';
end
if isfield(job,'use_stims')
    ust = job.use_stims;
else
    ust = '';
end
if isfield(job,'write_excel_text')
    if job.write_excel_text
        write_excel = 1;
    else
        write_excel = 0;
    end
else
    write_excel = 0;
end
IC = job.IC; %colors to include
PGM = [];

for SubjIdx=1:length(job.IOImat)
    try
        tic
        clear ROI onsets_list M
        %Load IOI.mat information
        [IOI IOImat dir_ioimat]= ioi_get_IOI(job,SubjIdx);
        
        if ~isfield(IOI,'dev')
            IOI.dev.TR = 0.2;
        end
        
        if ~isfield(IOI.res,'seriesOK')
            disp(['No extracted time series available for subject ' int2str(SubjIdx) ' ... skipping series extraction']);
        else
            if ~isfield(IOI.res,'meanOK') || job.force_redo
                
                %load ROI
                if ~isempty(job.ROImat)
                    load(job.ROImat{SubjIdx});
                else
                    try
                        load(IOI.ROI.ROIfname);
                    catch
                        ts1 = strfind(dir_ioimat,'Res');
                        ts2 = strfind(dir_ioimat,filesep);
                        ts3 = ts2(ts2>ts1);
                        tdir = dir_ioimat(1:ts3(2));
                        load(fullfile(tdir,'ROI.mat'));
                    end
                end
                if job.save_figures
                    dir_fig = fullfile(dir_ioimat,'fig');
                    if ~exist(dir_fig,'dir'),mkdir(dir_fig);end
                end
                
                %get stimulation information - Careful, here onset duration is ignored!
                Ns = length(IOI.sess_res);
                if IC.include_HbT
                    if ~isfield(IOI.color,'HbT')
                        IOI.color.HbT = 'T';
                        IOI.color.eng = [IOI.color.eng IOI.color.HbT];
                    end
                end
                %restric onsets
                [IOI onsets_list pars_list] = ioi_restrict_onsets(IOI,job,rmi,ust);
                %                 onsets_list=onsets_list{2};
                %Check whether there is the same number of onset types in
                %each session; this is a
                %check if averaging over all selected sessions is sensible
                %in which case, global_M = 1 (M is the number of onset types)
                if job.generate_global
                    global_M = 1; found_first_M = 0;
                    for s1=1:Ns
                        if all_sessions || sum(s1==selected_sessions)
                            M{s1} = length(onsets_list{s1});
                            if ~found_first_M
                                PGM = M{s1};
                                found_first_M = 1;
                            else
                                if ~(M{s1}==PGM) %PGM = possible_global_M;
                                    %no global M possible
                                    global_M = 0;
                                end
                            end
                        end
                    end
                else
                    global_M = 0;
                end
                
                %find max M to loop over
                maxM = length(onsets_list{1});
                if Ns > 1
                    for s1=2:Ns
                        if all_sessions || sum(s1==selected_sessions)
                            if length(onsets_list{s1}) > maxM
                                maxM = length(onsets_list{s1});
                            end
                        end
                    end
                end
                                               
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Compute mean over stimulations
                [IOI U0] = ioi_stim_mean_call(job,IOI,ROI,maxM,global_M,...
                    PGM,onsets_list,pars_list);
                IOI.res.meanOK = 1;
                save(IOImat,'IOI');
                
                %Define M structure for Expectation-Maximization algorithm
                %Fit difference of two gamma functions
                if job.extract_HRF
                    Ma = IOI.res.Ma; %mean of "after segments", by region, stimulus type, and by color and region
                    GMa = IOI.res.GMa;                
                    IOI = ioi_extract_HRF(job,IOI,ROI,maxM,global_M,Ma,GMa);
                    save(IOImat,'IOI');
                end
                
                if job.generate_figures || job.save_figures 
                    ioi_stim_figures(job,IOI,onsets_list,ROI,dir_fig);
                end
                if write_excel
                    ioi_stim_write_excel(job,IOI,onsets_list,ROI,dir_fig);
                end
            end
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