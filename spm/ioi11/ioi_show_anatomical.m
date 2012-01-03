function ioi_show_anatomical(handles)
subject_selected = handles.Info.subject_selected;
[filesRec,dummy] = spm_select('FPListRec',subject_selected,'_anat.nii');
S = ioi_get_common_settings(handles);
V = spm_vol(filesRec{1});
Y = spm_read_vols(V);
[nx ny] = size(Y);
set(handles.axes6,'Position',[S.XOffset+2*S.XShift S.YOffset S.scalex*nx S.scaley*ny]); %only to ensure no old axes are left behind
axes(handles.axes6) %select movie axes
imagesc(Y);
%axis(handles.axes6, 'off')
set(handles.axes6,'FontSize',handles.Movie.CommonFontSize);
%colormap(handles.axes6,gray) - does not work
title('Anatomical image (green)')  
end