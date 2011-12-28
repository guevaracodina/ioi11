function varargout = ioi_spatial_LPF(mode,K,data2D)
%Gaussian spatial LPF - based on NIRS_SPM function spm_filter_HPF_LPF_WMDL
switch mode
    case 'set'        
        %FWHM = 2*K.radius;
        sigma   = K.radius/2;
        h       = round(4*sigma);
        h       = exp(-(-h:h).^2/(2*sigma^2));
        n       = length(h);
        d       = (1:n) - (n + 1)/2;
        if      n == 1, h = 1; end
        
        k = K.k1;
        L = spdiags(ones(k,1)*h, d, k,k);
        K.K1 = spdiags(1./sum(L')',0,k,k)*L;
        k = K.k2;
        L = spdiags(ones(k,1)*h, d, k,k);
        K.K2 = spdiags(1./sum(L')',0,k,k)*L;
        K.Ks1 = K.K1*K.K1';
        K.Ks2 = K.K2*K.K2';
        varargout{1} = K;
    case 'lpf'
        varargout{1} =  K.Ks1 * data2D * K.Ks2';
    otherwise
end
end