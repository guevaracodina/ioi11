function ioi_display_frame(handles)
%get the desired frame:
frame = str2double(get(handles.edit_frame,'String'));
Y = handles.Movie.Y;
[nx ny nF] = size(Y);
S = ioi_get_common_settings(handles);
set(handles.axes1,'Position',[S.XOffset S.YOffset+S.YShift S.scalex*nx S.scaley*ny]); %only to ensure no old axes are left behind
axes(handles.axes1) %select movie axes
%produce one frame to obtain a colorbar for the movie
frame = ioi_check_frame(frame,nF);
[clims lmin lmax] = ioi_get_clims(handles);  
imagesc(squeeze(Y(:,:,frame)),clims);  
%axis(handles.axes1, 'off')
set(handles.axes1,'FontSize',handles.Movie.CommonFontSize);
colorbar('location','EastOutside')
sf = handles.Movie.AcqSamplingFreq;
title(['Frame at time ' num2str(frame/sf) ' s']) 
drawnow
end 