function doColor = ioi_doColor(IOI,c1,IC)
% Determines if color c1 is to be processed for this batch if selected in
% ioi_dfg_include_colors options.
%_______________________________________________________________________________
% Copyright (C) 2012 LIOM Laboratoire d'Imagerie Optique et Moleculaire
%                    Ecole Polytechnique de Montreal
%_______________________________________________________________________________
doColor = 1;
%potentially exclude various colors, to save time
if ~IC.include_OD
    if IOI.color.eng(c1) == IOI.color.red || ...
            IOI.color.eng(c1) == IOI.color.green || ...
            IOI.color.eng(c1) == IOI.color.yellow
        doColor = 0;
    end
end
try
    if IOI.color.eng(c1) == IOI.color.flow
        if ~IC.include_flow
            doColor = 0;
        end
    end
end
try
    if IOI.color.eng(c1) == IOI.color.laser
        doColor = 0; %for debugging, will turn it off eventually
    end
end
try
    if IOI.color.eng(c1) == IOI.color.contrast
        doColor = 0;
    end
end
try
    if ~IC.include_HbT
        if IOI.color.eng(c1) == IOI.color.HbT
            doColor = 0;
        end
    end
end
try
    if ~IC.include_HbR
        if IOI.color.eng(c1) == IOI.color.HbR
            doColor = 0;
        end
    end
end
try
    if ~IC.include_HbO
        if IOI.color.eng(c1) == IOI.color.HbO
            doColor = 0;
        end
    end
end
try
    if ~IC.include_CMRO2
        if IOI.color.eng(c1) == IOI.color.CMRO2
            doColor = 0;
        end
    end
end

% EOF
