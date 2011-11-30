function H = ioi_HDM_hrf(TR,x,d,color,c1,include_flow,fit_3_gamma)
H = [];
warning('off')
%nonlinear fit - two steps
%fit the hemodynamic increase (decrease for HbR)
% if strcmp(color.eng(c1),color.HbR)
%     d = -d;
% end
if ~strcmp(color.eng(c1),color.flow) || include_flow  
    H.U.dt = TR;
    H.Y.dt = TR;
    H.Y.y = d'/sum(d);
    u = zeros(1,length(x));
    u(1) = 1;
    H.U.u = u';
    H.f     = 'ioi_fx_hrf';
    if fit_3_gamma   
        p = [5.8 10 0.3 0.5 8 15 0.5 8 0];
        H.g     = 'ioi_gx_hrf_3gamma';
    else  
        p = [5.8 10 0.3 0.5 8 0];
        H.g     = 'ioi_gx_hrf';
    end
    H.x     = 0;
    H.n     = 1;
    %Number of inputs of direct model
    H.m=1;
    %choose priors
    H.pE = p;
    H.pC = eye(length(p));
    Y = H.Y;
    U = H.U;
    % nonlinear system identification
    %--------------------------------------------------------------------------
    [H.Ep,H.Cp,H.Eh,H.F] = spm_nlsi_GN(H,U,Y);
    %calculate the fitting curve
    H.yp = sum(d)*feval(H.g,x,[],H.Ep,[]);
end
warning('on')
end