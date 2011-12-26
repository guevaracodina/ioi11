function ioi_play_movie(handles)
obj = handles.Movie.obj;
F = handles.Movie.F;
%not using obj.FrameRate
% Play back the movie once at the video's frame rate.
set(handles.axes1,'Position',[10 10 obj.Width obj.Height]);
movie(handles.figure1, F, 1, handles.Movie.FrameRate,[10 10 0 0]); %obj.Width obj.Height]);
end