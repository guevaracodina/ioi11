function ioi_show_stats(handles)
Y = handles.Movie.Y;
[nx ny nF] = size(Y);
scalex = handles.Info.scalex;
scaley = handles.Info.scaley;
stat_mode_selected = handles.Info.stat_mode_selected;
sf = handles.Movie.AcqSamplingFreq;
switch stat_mode_selected
    case 'Time-to-Peak (+)'
        [values ind] = max(Y,[],3);
        pY = ind/sf; %in seconds
        tit = stat_mode_selected;
        clims = [0 nF/sf];
    case 'Time-to-Peak (-)'
        [values ind] = min(Y,[],3);
        pY = ind/sf;
        tit = stat_mode_selected;
        clims = [0 nF/sf];
    case 'Min'
        pY = min(Y,[],3);
        tit = 'Minimum along time direction';
        [clims lmin lmax] = ioi_get_clims(handles);    
    case 'Max'
        %just plot the max, time-wise
        pY = max(Y,[],3);
        tit = 'Maximum along time direction';
        [clims lmin lmax] = ioi_get_clims(handles);    
    case 'FWHM (+)'
        %start by finding max and time-to-peak
        [values dummy_ind] = max(Y,[],3);
        %find times to half-peaks
        half_peaks = values/2;
        ind_hp = Y>=repmat(half_peaks,[1 1 nF]);
        %find indices of first half-peaks found 
        [dummy_values_hp ind_hp_first] = max(ind_hp,[],3);
        %find indices of last half-peaks found
        ind_hp = ind_hp(:,:,end:-1:1);
        [dummy_values_hp ind_hp_last] = max(ind_hp,[],3);
        ind_hp_last = nF-ind_hp_last;
        %calculate FWHM
        pY = (ind_hp_last-ind_hp_first)/sf;               
        tit = stat_mode_selected;
    case 'FWHM (-)'
        [values dummy_ind] = min(Y,[],3);
        half_peaks = values/2;
        ind_hp = Y<=repmat(half_peaks,[1 1 nF]);
        %find indices of first half-peaks found 
        [dummy_values_hp ind_hp_first] = max(ind_hp,[],3);
        %find indices of last half-peaks found
        ind_hp = ind_hp(:,:,end:-1:1);
        [dummy_values_hp ind_hp_last] = max(ind_hp,[],3);
        ind_hp_last = nF-ind_hp_last;
        %calculate FWHM
        pY = (ind_hp_last-ind_hp_first)/sf;    
        tit = stat_mode_selected;
end
        
set(handles.axes5,'Position',[4+100+1000 12+400 scalex*nx scaley*ny]); 
axes(handles.axes5) %select movie axes
try 
    imagesc(pY,clims);
catch
    imagesc(pY);
end
axis(handles.axes5, 'off')
colorbar('location','WestOutside')
title(tit)  

end