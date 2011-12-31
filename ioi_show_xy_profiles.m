function ioi_show_xy_profiles(handles)
F = handles.Movie.F;
Y = handles.Movie.Y;
[nx ny nF] = size(Y);
%x-cut
ypos = round(ny/2);
set(handles.axes_cutx,'Position',[4+100 12 nx ny]); %only to ensure no old axes are left behind
axes(handles.axes_cutx) %select movie axes
imagesc(squeeze(Y(:,ypos,:)));
axis(handles.axes_cutx, 'off')
colorbar('location','WestOutside')
title(['x-cut at y=' int2str(ypos)])  
%y-cut
xpos = round(nx/2);
set(handles.axes_cuty,'Position',[4+100+600 12 nx ny]); %only to ensure no old axes are left behind
axes(handles.axes_cuty) %select movie axes
imagesc(squeeze(Y(xpos,:,:)));
axis(handles.axes_cuty, 'off')
colorbar('location','WestOutside')
title(['y-cut at x=' int2str(xpos)])
end