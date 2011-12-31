function ioi_open_movie(hObject,handles)
%check that there is a valid movie there
try
    movie_selected = handles.Info.movie_selected;
    load(movie_selected,'-mat');
    %store data
    handles.Movie.Y = d;
    Y = d; clear d
    guidata(hObject, handles);
    [clims lmin lmax] = ioi_get_clims(handles);   
    Nf = size(Y,3);
    F(Nf) = struct('cdata',[],'colormap',[]);
    h0 = figure;
    for i0=1:Nf
        imagesc(squeeze(Y(:,:,i0)),clims);
        F(i0) = getframe;
    end
    try close(h0); end
    %store raw data (Y) and movie data (F)
    handles.Movie.F = F;
    set(handles.text_min,'String',sprintf('%4.3s',lmin));
    set(handles.text_max,'String',sprintf('%4.3s',lmax));
    guidata(hObject, handles);
    handles.Movie.FrameRate = 90; %play the movie quickly
    ioi_play_movie(handles);
    %restore desired movie frequency
    handles.Movie.FrameRate = str2double(get(handles.movie_frequency,'String'));
    set(handles.message_box,'String','');
    set(handles.message_box,'BackgroundColor','White');
catch
    set(handles.message_box,'String','Movie cannot be opened');
    set(handles.message_box,'BackgroundColor','Red');
end
guidata(hObject, handles);