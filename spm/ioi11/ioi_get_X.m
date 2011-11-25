function [X U] = ioi_get_X(IOI,name,ons,dur,s1,bases,volt)
SPM.xBF.dt = IOI.dev.TR;
SPM.xBF.T = 1;
SPM.xBF.T0 = 0;
SPM.xBF.UNITS = 'secs';
nambase = fieldnames(bases);
if ischar(nambase)
    nam=nambase;
else
    nam=nambase{1};
end
if strcmp(fieldnames(bases),'hrf') || strcmp(fieldnames(bases),'rat') ...
        || strcmp(fieldnames(bases),'mouse')
    %canonical HRF or rat or mouse
    if all(bases.(nam).derivs == [0 0])
        SPM.xBF.name = 'hrf';
    elseif all(bases.(nam).derivs == [1 0])
        SPM.xBF.name = 'hrf (with time derivative)';
    elseif all(bases.(nam).derivs == [1 1])
        SPM.xBF.name = 'hrf (with time and dispersion derivatives)';
    else
        disp('Unrecognized hrf derivative choices.')
    end
    SPM.xBF.nam = nam;
else
    switch nam,
        case 'fourier',
            SPM.xBF.name = 'Fourier set';
        case 'fourier_han',
            SPM.xBF.name = 'Fourier set (Hanning)';
        case 'gamma',
            SPM.xBF.name = 'Gamma functions';
        case 'fir',
            SPM.xBF.name = 'Finite Impulse Response';
        otherwise
            error('Unrecognized hrf derivative choices.')
    end
    SPM.xBF.length = bases.(nam).length;
    SPM.xBF.order  = bases.(nam).order;
end

SPM.xBF = spm_get_bf_rat_mouse(SPM.xBF);
if size(SPM.xBF.bf,1) == 1 || size(SPM.xBF.bf,2) == 1
    SPM.xBF.bf = SPM.xBF.bf/sum(SPM.xBF.bf); %normalize
end

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
% Convolve stimulus functions with basis functions
%------------------------------------------------------------------
[X,Xn,Fc] = spm_Volterra(U,SPM.xBF.bf,volt);

% Resample regressors at acquisition times (32 bin offset)
%-------------------------------------------------
X = X((0:(SPM.nscan - 1))*SPM.xBF.T + SPM.xBF.T0 + 32,:);

% and orthogonalise (within trial type)
%--------------------------------------
for i = 1:length(Fc)
    X(:,Fc(i).i) = spm_orth(X(:,Fc(i).i));
end
end