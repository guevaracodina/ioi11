function O = ioi_get_PS(IOI,O)
%Loop to find HbO, HbR and Flow positions in ROI
for c1=1:length(IOI.color.eng)
    if IOI.color.eng(c1)==IOI.color.HbO, cHbO = c1; end
    if IOI.color.eng(c1)==IOI.color.HbR, cHbR = c1; end
    if IOI.color.eng(c1)==IOI.color.flow, cFlow = c1; end
end
O.cHbO = cHbO;
O.cHbR = cHbR;
O.cFlow = cFlow;