function S = ioi_save_simu(M)
%for simulations: store results in place of subjects
if M.S.simuOn
    S.EpS = M.P;
    S.pE = M.pE;
    S.Ep = M.Ep;
    S.Cp = M.Cp;
    S.F = M.F;
end