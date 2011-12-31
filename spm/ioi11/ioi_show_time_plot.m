function ioi_show_time_plot(handles)
Y = handles.Movie.Y;
[nx ny nF] = size(Y);
scalex = handles.Info.scalex;
scaley = handles.Info.scaley;
%time plot
xpos = str2double(get(handles.edit_xpos,'String')); 
ypos = str2double(get(handles.edit_ypos,'String'));
[xpos ypos] = ioi_check_pos(xpos,ypos,nx,ny);
set(handles.axes_plot,'Position',[4+100+500 12+400 scalex*nx scaley*ny]); %only to ensure no old axes are left behind
axes(handles.axes_plot) %select movie axes
lp = linspace(0,nF/handles.Movie.AcqSamplingFreq,nF);
plot(lp,squeeze(Y(xpos,ypos,:)))
%axis(handles.axes5, 'off')
title(['Time series at (' int2str(xpos) ',' int2str(ypos) ')']) 
xlabel('Time (s)')
end