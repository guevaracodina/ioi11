function U = ioi_get_U(IOI,name,ons,dur,s1)
SPM = [];
SPM.xBF.dt = IOI.dev.TR;
SPM.xBF.T = 1;
SPM.xBF.T0 = 1;
SPM.xBF.UNITS = 'secs';
% Get inputs, neuronal causes or stimulus functions U
%------------------------------------------------------------------
SPM.nscan = IOI.sess_res{s1}.n_frames;
P.name = 'none';
P.h    = 0;
if isempty(name)
    SPM.Sess.U.name = {'Spk'};
else
    SPM.Sess.U.name = {name};
end
SPM.Sess.U.ons = ons;
SPM.Sess.U.dur = dur;
SPM.Sess.U.P = P;
U = spm_get_ons(SPM,1);