function ioi_show_time_plot(handles)
Y = handles.Movie.Y;
[nx ny nF] = size(Y);
S = ioi_get_common_settings(handles);
%time plot
xpos = str2double(get(handles.edit_xpos,'String')); 
ypos = str2double(get(handles.edit_ypos,'String'));
[xpos ypos] = ioi_check_pos(xpos,ypos,nx,ny);
set(handles.axes_plot,'Position',[S.XOffset+S.XShift S.YOffset+S.YShift S.scalex*nx S.scaley*ny]); %only to ensure no old axes are left behind
axes(handles.axes_plot) %select movie axes
lp = linspace(0,nF/handles.Movie.AcqSamplingFreq,nF);
drawnow expose update
plot(lp,squeeze(Y(xpos,ypos,:)))
set(handles.axes_plot,'FontSize',handles.Movie.CommonFontSize);
radius = round(str2double(get(handles.edit_spatialLPF,'String')));
title(['Time series at (' int2str(xpos) ',' int2str(ypos) ')' ' Radius: ' int2str(radius)]) 
sf = handles.Movie.AcqSamplingFreq;
set(handles.axes_plot,'xtick',0:nF/sf/5:nF/sf)
%drawnow expose update
%xlabel('Time (s)')
end