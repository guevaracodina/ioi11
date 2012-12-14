function [IOI IOImat dir_ioimat] = ioi_get_IOI(job,SubjIdx)
IOI = [];
IOImat = job.IOImat{SubjIdx};
[dir_ioimat dummy] = fileparts(job.IOImat{SubjIdx});
if isfield(job.IOImatCopyChoice,'IOImatCopy')
    newDir = job.IOImatCopyChoice.IOImatCopy.NewIOIdir;
    newDir = fullfile(dir_ioimat,newDir);
    if ~exist(newDir,'dir'),mkdir(newDir); end
    IOImat = fullfile(newDir,'IOI.mat');
    dir_ioimat = newDir;
end
try
    load(IOImat);
    display([IOImat ' now loaded']);
catch
    % XXX: This is not true, if we create a new dir, then there should be no
    % error... 
    load(job.IOImat{SubjIdx});
    display([IOImat ' not found -- ' job.IOImat{SubjIdx} ' now loaded']);
end