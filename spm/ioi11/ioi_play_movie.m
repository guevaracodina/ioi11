function ioi_play_movie(handles)
F = handles.Movie.F;
% Play back the movie once at set frame rate.
Y = handles.Movie.Y;
[nx ny nF] = size(Y);
scalex = handles.Info.scalex;
scaley = handles.Info.scaley;
set(handles.axes1,'Position',[4+100 12+400 scalex*nx scaley*ny]); %only to ensure no old axes are left behind
axes(handles.axes1) %select movie axes
%produce one frame to obtain a colorbar for the movie
[clims lmin lmax] = ioi_get_clims(handles);    
imagesc(squeeze(Y(:,:,1)),clims);
axis(handles.axes1, 'off')
colorbar('location','WestOutside')
movie(handles.figure1, F, 1, handles.Movie.FrameRate,[10+90 10+400 0 0]); 
ioi_show_xy_profiles(handles);
ioi_show_anatomical(handles);
ioi_show_time_plot(handles);
ioi_show_stats(handles);
end