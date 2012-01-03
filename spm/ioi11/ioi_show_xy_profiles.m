function ioi_show_xy_profiles(handles)
%Unfortunately, due to Matlab (row, column) convention,
%the variables x and y are inverted in the code compared to their usual meaning
%Namely: x is rows and therefore corresponds to vertical on the images
%y is columns, and is horizontal on the images
%Also, the images increment from left to right, but from top to bottom
Y = handles.Movie.Y;
[nx ny nF] = size(Y);
S = ioi_get_common_settings(handles);
%vertical cut
ypos = str2double(get(handles.edit_ypos,'String'));
xpos = str2double(get(handles.edit_xpos,'String'));
[xpos ypos] = ioi_check_pos(xpos,ypos,nx,ny);
set(handles.axes_cutx,'Position',[S.XOffset S.YOffset S.scalex*nx S.scaley*ny]); %only to ensure no old axes are left behind
axes(handles.axes_cutx) %select movie axes
imagesc(squeeze(Y(:,ypos,:)));
set(handles.axes_cutx,'FontSize',handles.Movie.CommonFontSize);
%axis(handles.axes_cutx, 'off')
colorbar('location','EastOutside')
title(['Vertical cut, at x=' int2str(ypos)])  
drawnow
%horizontal cut
%make it a little smaller
S.scaley = S.scaley*0.95;
set(handles.axes_cuty,'Position',[S.XOffset+S.XShift S.YOffset S.scalex*nx S.scaley*ny]); %only to ensure no old axes are left behind
axes(handles.axes_cuty) %select movie axes
imagesc(squeeze(Y(xpos,:,:)));
set(handles.axes_cuty,'FontSize',handles.Movie.CommonFontSize);
%axis(handles.axes_cuty, 'off')
colorbar('location','EastOutside')
title(['Horizontal cut, at y=' int2str(xpos)])
drawnow
end