function ioi_show_xy_profiles(handles)
Y = handles.Movie.Y;
[nx ny nF] = size(Y);
scalex = handles.Info.scalex;
scaley = handles.Info.scaley;
%x-cut
ypos = str2double(get(handles.edit_ypos,'String'));
xpos = str2double(get(handles.edit_xpos,'String'));
[xpos ypos] = ioi_check_pos(xpos,ypos,nx,ny);
set(handles.axes_cutx,'Position',[4+100 12 scalex*nx scaley*ny]); %only to ensure no old axes are left behind
axes(handles.axes_cutx) %select movie axes
imagesc(squeeze(Y(:,ypos,:)));
axis(handles.axes_cutx, 'off')
colorbar('location','WestOutside')
title(['x-cut at y=' int2str(ypos)])  
%y-cut
set(handles.axes_cuty,'Position',[4+100+500 12 scalex*nx scalex*ny]); %only to ensure no old axes are left behind
axes(handles.axes_cuty) %select movie axes
imagesc(squeeze(Y(xpos,:,:)));
axis(handles.axes_cuty, 'off')
colorbar('location','WestOutside')
title(['y-cut at x=' int2str(xpos)])
end