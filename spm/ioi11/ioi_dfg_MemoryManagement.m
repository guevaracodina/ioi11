function MemoryManagementMenu = ioi_dfg_MemoryManagement(varargin)
% Graphical interface configuration function for memory managament options.
% This code is part of a batch job configuration system for MATLAB. See help
% matlabbatch for a general overview.
% SYNTAX:
% MemoryManagementMenu      = ioi_dfg_MemoryManagement(loadAll)
% INPUTS:
% loadAll                   [OPTIONAL] If true, load all at once is chosen.
% OUTPUTS:
% MemoryManagementMenu      Batch job configuration menu
%_______________________________________________________________________________
% Copyright (C) 2012 LIOM Laboratoire d'Imagerie Optique et Moleculaire
%                    Ecole Polytechnique de Montreal
%_______________________________________________________________________________

% ------------------------------------------------------------------------------
% Optional inputs handling
% ------------------------------------------------------------------------------
% only want 1 optional input at most
numvarargs                  = length(varargin);
if numvarargs > 1
    error('ioi_dfg_MemoryManagement:TooManyInputs', ...
        'Requires at most 1 optional inputs');
end
% set defaults for optional inputs
optargs                     = {1};
% now put these defaults into the optargs cell array, and overwrite the ones
% specified in varargin.
optargs(1:numvarargs)       = varargin;
% Place optional args in memorable variable names
[loadAll]                   = optargs{:};
% ------------------------------------------------------------------------------

MemoryManagementMenu        = cfg_menu;
MemoryManagementMenu.tag    = 'MemoryManagementMenu';
MemoryManagementMenu.name   = 'Memory Management';
MemoryManagementMenu.labels = {'Load all at once','Load one at a time'};
MemoryManagementMenu.values = {1,0};
MemoryManagementMenu.val    = {loadAll};
MemoryManagementMenu.help   = {'Load all images at once (faster but requires'
    'more memory, or load one image at a time, to compute concentrations.'}';

% EOF
