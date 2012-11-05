function varargout = ioi_get_defaults(defstr, varargin)
% Get/set the defaults values associated with an identifier
% FORMAT defval = ioi_get_defaults(defstr)
% Return the defaults value associated with identifier "defstr". 
% Currently, this is a '.' subscript reference into the global  
% "defaults" variable defined in spm_defaults.m.
%
% FORMAT ioi_get_defaults(defstr, defval)
% Sets the nirs8 value associated with identifier "defstr". The new
% ioi8 value applies immediately to:
% * new modules in batch jobs
% * modules in batch jobs that have not been saved yet
% This value will not be saved for future sessions of SPM. To make
% persistent changes, edit cg_vbm8_defaults.m.
%_______________________________________________________________________
% Copyright (C) 2010 Laboratoire d'Imagerie Optique et Moleculaire
%
% based on Christian Gaser version of
% cg_vbm8_get_defaults $Id: cg_vbm8_get_defaults.m 2696 2009-02-05
% 20:29:48Z guillaume $
%
% Clément Bonnéry

global ioi11;
if isempty(ioi11)
    %load various default files
    ioi_defaults;
    ioi_msioi_defaults;
    ioi_concentration_defaults;
    ioi_flow_defaults;
    ioi_cmro2_defaults;
end

% construct subscript reference struct from dot delimited tag string
tags = textscan(defstr,'%s', 'delimiter','.');
subs = struct('type','.','subs',tags{1}');

if nargin == 1
    varargout{1} = subsref(ioi11, subs);
else
    ioi11 = subsasgn(ioi11, subs, varargin{1});
end
