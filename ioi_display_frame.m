function ioi_display_frame(handles)
%get the desired frame:
frame = str2double(get(handles.edit_frame,'String'));
Y = handles.Movie.Y;
[nx ny nF] = size(Y);
scalex = handles.Info.scalex;
scaley = handles.Info.scaley;
set(handles.axes1,'Position',[4+100 12+400 scalex*nx scaley*ny]); %only to ensure no old axes are left behind
axes(handles.axes1) %select movie axes
%produce one frame to obtain a colorbar for the movie
frame = ioi_check_frame(frame,nF);
[clims lmin lmax] = ioi_get_clims(handles);  
imagesc(squeeze(Y(:,:,frame)),clims);  
axis(handles.axes1, 'off')
colorbar('location','WestOutside')
end