% __________________________________________________________________________  
% SPM8-Compatible Intrinsic Optical Imaging Toolbox  
% Version0(IOI11)2011-09-01  
% __________________________________________________________________________  
% Copyright (C) 2011 LIOM - Ecole Polytechnique Montreal  
%   
% ==========================================================================  
% Description  
% ==========================================================================  
% This toolbox is available to the scientific   
% community under the terms of the GNU General Public License.  
%   
% ==========================================================================  
% How to use  
% ==========================================================================  
%  
% The initial reading module is only designed for intrinsic optical imaging   
% acquisitions with the system designed at LIOM.   
%  
% To use IOI11 on data acquired with the LIOM IOI system(s), some steps  
% must be taken to ensure integrity of the data prior to processing:  
% an approximately equal number of images must be present in the colors   
% available (Green, Red, Yellow); the images must have the same size  
%  
% Prior to proper processing, images of different sizes can be brought to   
% one size by using the IOI_shrink module - this is the only module that  
% writes directly over raw data   
%  
% Raw data should be placed according to the following directory structure:  
% first a group level directory, inside which a group directory for raw data,  
% inside which a directory for each subject, inside which a pair of directories   
% for each session, of which one is for all the images in binary format  
%  
% Be careful not to put extraneous files in the raw data directories; as the  
% code attempts automatic detection of files, it will get confused if more  
% files are present  
%  
% Color names for raw image files are in French, later translated into English  
% in the code  
%  
% ==========================================================================  
% Technical information  
% ==========================================================================  
%  
% Modules are called by a .m file ending with _run and have a graphical  
% interface configured by a file ending with _cfg.   
