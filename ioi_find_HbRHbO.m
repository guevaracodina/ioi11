function [cHbR cHbO] = ioi_find_HbRHbO(IOI,s1)
%find colors for HbO and HbR
for c0=1:length(IOI.sess_res{s1}.fname)
    if IOI.color.eng(c0) == IOI.color.HbO
        cHbO = c0;
    end
    if IOI.color.eng(c0) == IOI.color.HbR
        cHbR = c0;
    end
end