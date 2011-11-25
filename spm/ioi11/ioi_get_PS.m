function PS = ioi_get_PS(IOI,includeHbR,includeHbT,includeFlow,PhysioModel_Choice)
%Loop to find HbO, HbR and Flow positions in ROI
for c1=1:length(IOI.color.eng)
    if IOI.color.eng(c1)==IOI.color.HbO, cHbO = c1; end
    if IOI.color.eng(c1)==IOI.color.HbR, cHbR = c1; end
    if IOI.color.eng(c1)==IOI.color.flow, cFlow = c1; end
end
PS.xY.cHbO = cHbO;
PS.xY.cHbR = cHbR;
PS.xY.cFlow = cFlow;
PS.xY.includeHbR = includeHbR;
PS.xY.includeHbT = includeHbT;
PS.xY.includeFlow = includeFlow;
PS.PhysioModel_Choice = PhysioModel_Choice;