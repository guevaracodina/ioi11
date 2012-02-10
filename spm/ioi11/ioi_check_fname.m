function fname = ioi_check_fname(fname0,dir_ioimat,s1)
try
    V = spm_vol(fname0{1});
    fname = fname0;
catch
    for f1=1:length(fname0)
        [dir0 fil0 ext0] = fileparts(fname0{f1});
        fsep = strfind(dir_ioimat,filesep);
        res = strfind(dir_ioimat,'Res');
        fgr = fsep(fsep > res);
        tdir = dir_ioimat(1:fgr(2));
        fname{f1} = fullfile(tdir,['S' gen_num_str(s1,2)],[fil0 ext0]);
    end
end
                            
                            