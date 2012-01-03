function ioi_play_movie(handles)
F = handles.Movie.F;
% Play back the movie once at set frame rate.
Y = handles.Movie.Y;
[nx ny nF] = size(Y);
S = ioi_get_common_settings(handles);
set(handles.axes1,'Position',[S.XOffset S.YOffset+S.YShift S.scalex*nx S.scaley*ny]); %only to ensure no old axes are left behind
axes(handles.axes1) %select movie axes
%produce one frame to obtain a colorbar for the movie
[clims lmin lmax] = ioi_get_clims(handles);    
imagesc(squeeze(Y(:,:,1)),clims);
%axis(handles.axes1, 'off')
set(handles.axes1,'FontSize',handles.Movie.CommonFontSize);
%set(handles.axes1,'YAxisLocation','Right');
colorbar('location','EastOutside')
movie(handles.figure1, F, 1, handles.Movie.FrameRate,[S.XOffset S.YOffset+S.YShift 0 0]); 
ioi_show_xy_profiles(handles);
try ioi_show_anatomical(handles); end
ioi_show_time_plot(handles);
ioi_show_stats(handles);
end