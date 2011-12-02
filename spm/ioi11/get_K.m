function K = get_K(svec,fwhm,TR)
K.HParam.type = 'none';
K.LParam.FWHM = fwhm;
K.LParam.type = 'Gaussian';
K = struct( 'HParam', K.HParam,...
    'row', svec ,...
    'RT', TR,...
    'LParam', K.LParam);
K = spm_filter_HPF_LPF_WMDL(K);
K.row = 1:length(svec);
end
