function IOI = ioi_get_elinfo(IOI)
% Gets .mat file names with elinfo data (ECG, pressure, etc.)
% SYNTAX
% IOI = ioi_get_elinfo(IOI)
% INPUT 
% IOI       Matrix with no elinfo saved names
% OUTPUT 
% IOI       Matrix with updated elinfo filenames
% 
%_______________________________________________________________________________
% Copyright (C) 2012 LIOM Laboratoire d'Imagerie Optique et Moléculaire
%                    École Polytechnique de Montréal
%_______________________________________________________________________________

% Look for elinfo files in the results top folder
dirStruct = dir(fullfile(IOI.dir.dir_subj_res,'elinfo_S*.mat'));
% Exclude those files that dir cannot query
dirStruct = dirStruct(~cellfun(@isempty,{dirStruct(:).date}));

% If no elinfo files were found search directory recursively
if isempty(dirStruct)
    % Save current folder
    currentFolder = pwd;
    % Go to results folder
    cd(IOI.dir.dir_subj_res);
    % Look inside the sub-folders
    [status, dirStruct] = system( 'dir elinfo_S*.mat /s /b' );
    % Break text into cell of rows
    dirStruct = textscan( dirStruct, '%s', 'delimiter', '\n' );
    % Parse output into the same format as IOI.res.elinfo
    elinfoList = dirStruct{1}';
    % Go back to current folder
    cd(currentFolder)
else
    % Get cell with the elinfo filenames
    elinfoList = cellfun(@fullfile, repmat(...
        {IOI.dir.dir_subj_res}, [1 length(dirStruct)]),...
        {dirStruct(:).name},'UniformOutput',false);
end

% Update elinfo in the IOI matrix
IOI.res.elinfo = elinfoList;

% EOF
