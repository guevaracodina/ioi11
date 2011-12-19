function M = ioi_set_simu(M,it1)
%Buxton-Friston balloon and other models
M.P = M.S.pA(it1,:)';
%efficacies
M.P(end-size(M.U.u,2)+M.S.simuS) = 1;
%which dt is used in there?
ys = spm_int_J(M.P,M,M.U); %6.7 times slower than spm_int_D, but produces a very different result
Y0 = M.Y0;
%Upsample
if M.S.simuUpsample > 1
    Y0.y = interp(Y0.y,M.S.simuUpsample);
    for iX0=1:size(Y0.X0,2)
        tmp0(:,iX0) = interp(Y0.X0(:,iX0),M.S.simuUpsample);
    end
    Y0.X0 = tmp0;
end
ns = size(Y0.y,1);
Y0.dt = Y0.dt/M.S.simuUpsample;
ys = ys(round((0:(ns - 1))*Y0.dt/U.dt)+1);

if M.S.simuNoise
    Y.y = Y0.y + ys*M.S.simuA/100;
else
    Y.y = ys*M.S.simuA/100;
end
%Low pass filtering of the data after downsampling -- otherwise there will be aliasing
if M.S.simuUpsample < 16 %here Y0 might get LPF twice...
    Y.y = ButterLPF(1/Y.dt,0.95*1/(2*Y.dt),3,Y.y);
end
M.Y = Y;
end