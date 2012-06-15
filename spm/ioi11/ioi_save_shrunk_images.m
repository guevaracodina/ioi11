function IOI = ioi_save_shrunk_images(IOI,job,SH,dir_ioimat)
shrink_x = SH.shrink_x;
shrink_y = SH.shrink_y;
IC = job.IC;
%select a subset of sessions
[all_sessions selected_sessions] = ioi_get_sessions(job);

for s1=1:length(IOI.sess_res)
    if all_sessions || sum(s1==selected_sessions)
        try
            IOI.sess_shrunk{s1};
            %shrunk images
            if SH.force_shrink_recompute
                %force recompute
                IOI.sess_shrunk{1000}; %will certainly break
            end
        catch
            for c1=1:length(IOI.color.eng)
                doColor = ioi_doColor(IOI,c1,IC);
                fname0 = {};
                if doColor
                    if IOI.color.eng(c1) == 'T'
                        doHbT = 1;
                        [cHbR cHbO] = ioi_find_HbRHbO(IOI,s1);
                        fname = IOI.sess_res{s1}.fname{cHbR};
                        fname2 = IOI.sess_res{s1}.fname{cHbO};
                        fname2 = ioi_check_fname(fname2,dir_ioimat,s1);
                    else
                        doHbT = 0;
                        fname = IOI.sess_res{s1}.fname{c1};
                    end
                    fname = ioi_check_fname(fname,dir_ioimat,s1);
                    for f1=1:length(fname)
                        
                        V = spm_vol(fname{f1});
                        Y = spm_read_vols(V);
                        if doHbT
                            V2 = spm_vol(fname2{f1});
                            Y2 = spm_read_vols(V2);
                            Y = Y+Y2;
                        end
                        
                        %shrink by averaging
                        Y0 = zeros(size(Y(1:shrink_x:(end-shrink_x+1),1:shrink_y:(end-shrink_y+1),:)));
                        for i1=1:shrink_x
                            for i2=1:shrink_y
                                Y0 = Y0 + Y(i1:shrink_x:(end-shrink_x+i1),i2:shrink_y:(end-shrink_y+i2),:);
                            end
                        end
                        %save images
                        [dir0 fil0 ext0] = fileparts(fname{f1});
                        if doHbT
                            fil0 = regexprep(fil0, ['_' IOI.color.HbR '_'],  ['_' IOI.color.HbT '_']);
                        end
                        tn = fullfile(dir0,[fil0 '_shrunk_' int2str(shrink_x) 'x' int2str(shrink_y) ext0]);
                        fname0 = [fname0; tn];
                        Y0 = reshape(Y0,[size(Y0,1) size(Y0,2) 1 size(Y0,3)]);
                        ioi_save_nifti(Y0,tn,[1 1 1]);
                    end
                end
                IOI.sess_shrunk{s1}.fname{c1} = fname0;
            end
        end
    end
end