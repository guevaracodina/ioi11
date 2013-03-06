function [ROI IOI] = ioi_extract_core(IOI,job,mask,Amask,varargin)
% ROI time course extraction over colors and files. Calls ioi_extract_main
% SYNTAX:
% [ROI IOI] = ioi_extract_core(IOI, job, mask, Amask, dataType)
% INPUTS:
% IOI       Matrix with the IOI structure
% job       Matlab batch job structure
% mask      ROI/seed binary mask
% Amask     Activation mask
% [dataType]
%           'rawData'       - Default: extract ROIs/seeds raw data
%           'filtData'      - extract ROIs/seeds temporally BPF data
%           'regressData'   - extract ROIs/seeds GLM-regressed data
% OUTPUTS:
% ROI       Cell with ROIs/seeds time traces
% IOI       Matrix with the IOI structure
%_______________________________________________________________________________
% Copyright (C) 2012 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________

% only want 1 optional input at most
numVarArgs = length(varargin);
if numVarArgs > 1
    error('ioi_extract_core:TooManyInputs', ...
        'requires at most 1 optional inputs: dataType');
end
% set defaults for optional inputs ()
optArgs = {'rawData'};
% skip any new inputs if they are empty
newVals = cellfun(@(x) ~isempty(x), varargin);
% now put these defaults into the optArgs cell array, and overwrite the ones
% specified in varargin.
optArgs(newVals) = varargin(newVals);
% Place optional args in memorable variable names
[dataType] = optArgs{:};

[all_sessions selected_sessions] = ioi_get_sessions(job);
[all_ROIs selected_ROIs] = ioi_get_ROIs(job);
IC = job.IC;
%loop over sessions
for s1=1:length(IOI.sess_res)
    if all_sessions || sum(s1==selected_sessions)
        %loop over available colors
        for c1=1:length(IOI.sess_res{s1}.fname)
            doColor = ioi_doColor(IOI,c1,IC);
            if doColor
                colorOK = 1;
                if ~(IOI.color.eng(c1)==IOI.color.laser)
                    % 
                    switch dataType
                        case 'rawData'
                            % raw data filenames
                            fname_list = IOI.sess_res{s1}.fname{c1};
                        case 'filtData'
                            % filtered data filenames
                            fname_list = IOI.fcIOS.filtNdown.fnameWholeImage(s1,c1);
                        case 'regressData'
                            % GLM regressed data filenames
                            fname_list = IOI.fcIOS.SPM.fname(s1,c1);
                        otherwise
                            % raw data filenames
                            fname_list = IOI.sess_res{s1}.fname{c1};
                    end
                    
                    %skip laser - only extract for flow
                    % fname_list = IOI.sess_res{s1}.fname{c1}; %//EGC
                    %initialize
                    
                    if job.extractBrainMask && job.extractingBrainMask
                        nROI = 1; % Only 1 brain mask
                    else
                        nROI = 1:length(IOI.res.ROI); % All the ROIs
                    end

                    for r1 = nROI
                        if all_ROIs || sum(r1==selected_ROIs)
                            ROI{r1}{s1,c1} = [];
                        end
                    end
                    %loop over files
                    for f1=1:length(fname_list)
                        try
                            fname = fname_list{f1};
                            vols = spm_vol(fname);
                            d = spm_read_vols(vols);
                            [d1 d2 d3 d4] = size(d);
                            if d1 <= 1 || d2 <= 1
                                colorOK = 0;
                            end
                        catch
                            colorOK = 0;
                        end
                        %time dimension in 3rd dimension for colors
                        %R, G, Y, but in 4th dimension for O, D, F
                        %Loop over ROIs
                        [IOI ROI] = ioi_extract_main(IOI,ROI,job,d,d3,d4,c1,s1,colorOK,mask,Amask);
                    end
                    if colorOK
                        disp(['ROIs for session ' int2str(s1) ' and color ' IOI.color.eng(c1) ' completed']);
                    end
                end
            end
        end
    end
end
