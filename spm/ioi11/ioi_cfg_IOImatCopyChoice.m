function IOImatCopyChoice = ioi_cfg_IOImatCopyChoice(dir_name)
IOImatOverwrite         = cfg_branch;
IOImatOverwrite.tag     = 'IOImatOverwrite';
IOImatOverwrite.name    = 'Overwrite IOI.mat structure'; 
IOImatOverwrite.help    = {'Will not copy IOI structure.'
            'This will write over the previous NIRS.mat'}';

NewIOIdir         = cfg_entry;
NewIOIdir.name    = 'New directory for IOI.mat';
NewIOIdir.tag     = 'NewIOIdir';       
NewIOIdir.strtype = 's';
NewIOIdir.val{1}    = dir_name; %'Flow';
NewIOIdir.num     = [1 Inf];     
NewIOIdir.help    = {'Directory for IOI.mat.'}'; 

IOImatCopy         = cfg_branch;
IOImatCopy.tag     = 'IOImatCopy';
IOImatCopy.name    = 'Create new directory and copy IOI structure there'; 
IOImatCopy.val     = {NewIOIdir};
IOImatCopy.help    = {'Create new directory and copy IOI structure there.'}';
        
%Common to most modules: for creating a new directory and copying IOI.mat
IOImatCopyChoice           = cfg_choice;
IOImatCopyChoice.name      = 'Choose IOI copy method';
IOImatCopyChoice.tag       = 'IOImatCopyChoice';
IOImatCopyChoice.values    = {IOImatOverwrite IOImatCopy}; 
IOImatCopyChoice.val       = {IOImatOverwrite}; 
IOImatCopyChoice.help      = {'Choose whether to overwrite the IOI.mat structure'
            'or to create a new directory'
            'and copy the IOI.mat structure there'}'; 

        IOImatOverwrite         = cfg_branch;
IOImatOverwrite.tag     = 'IOImatOverwrite';
IOImatOverwrite.name    = 'Overwrite IOI.mat structure'; 
IOImatOverwrite.help    = {'Will not copy IOI structure.'
            'This will write over the previous NIRS.mat'}';

NewIOIdir         = cfg_entry;
NewIOIdir.name    = 'New directory for IOI.mat';
NewIOIdir.tag     = 'NewIOIdir';       
NewIOIdir.strtype = 's';
NewIOIdir.val{1}    = 'GLM';
NewIOIdir.num     = [1 Inf];     
NewIOIdir.help    = {'Directory for IOI.mat.'}'; 

IOImatCopy         = cfg_branch;
IOImatCopy.tag     = 'IOImatCopy';
IOImatCopy.name    = 'Create new directory and copy IOI structure there'; 
IOImatCopy.val     = {NewIOIdir};
IOImatCopy.help    = {'Create new directory and copy IOI structure there.'}';
        
%Common to most modules: for creating a new directory and copying IOI.mat
IOImatCopyChoice           = cfg_choice;
IOImatCopyChoice.name      = 'Choose IOI copy method';
IOImatCopyChoice.tag       = 'IOImatCopyChoice';
IOImatCopyChoice.values    = {IOImatOverwrite IOImatCopy}; 
IOImatCopyChoice.val       = {IOImatOverwrite}; 
IOImatCopyChoice.help      = {'Choose whether to overwrite the IOI.mat structure'
            'or to create a new directory'
            'and copy the IOI.mat structure there'}'; 

        