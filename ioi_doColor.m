function doColor = ioi_doColor(IOI,c1,IC)
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
        doColor = 0;
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