function ioi_open_movie(hObject,handles)
%check that there is a valid movie there
try
    movie_selected = handles.Info.movie_selected;
    load(movie_selected,'-mat');
    %store data
    handles.Movie.Y = d;
    Y = d; clear d
    if isfield(handles.Movie,'Y0')
        handles.Movie = rmfield(handles.Movie,'Y0'); %clean up
    end
    guidata(hObject, handles);
    [clims lmin lmax] = ioi_get_clims(handles);   
    [nx ny nF] = size(Y);
    F(nF) = struct('cdata',[],'colormap',[]);
    h0 = figure;
    for i0=1:nF
        imagesc(squeeze(Y(:,:,i0)),clims);
        F(i0) = getframe;
    end
    try close(h0); end
    %store raw data (Y) and movie data (F)
    scalex = 400/nx;
    scaley = 300/ny;
    handles.Info.scalex = scalex;
    handles.Info.scaley = scaley;
    handles.Movie.F = F;
    %initialize the min and max boxes 
%     set(handles.text_min,'String',sprintf('%4.3s',lmin));
%     set(handles.text_max,'String',sprintf('%4.3s',lmax));
    set(handles.text_min,'String',num2str(lmin));
    set(handles.text_max,'String',num2str(lmax));
    guidata(hObject, handles);
    ioi_update_min(hObject, handles);
    ioi_update_max(hObject, handles);
    %sliders
    set(handles.slider_frame,'max',nF);
    set(handles.slider_xpos,'max',nx);
    set(handles.slider_ypos,'max',ny);
    set(handles.slider_frame,'min',1);
    set(handles.slider_xpos,'min',1);
    set(handles.slider_ypos,'min',1);
    frame = str2double(get(handles.edit_frame,'String'));
    xpos = str2double(get(handles.edit_xpos,'String'));
    ypos = str2double(get(handles.edit_ypos,'String'));
    radius = str2double(get(handles.edit_spatialLPF,'String'));
    set(handles.slider_frame,'value',frame);
    set(handles.slider_xpos,'value',xpos);
    set(handles.slider_ypos,'value',ypos);
    set(handles.slider_spatialLPF,'value',radius);
    set(handles.slider_frame,'SliderStep',[1/nF 10/nF]);
    set(handles.slider_xpos,'SliderStep',[1/nx 10/nx]);
    set(handles.slider_ypos,'SliderStep',[1/ny 10/ny]);
    handles.Movie.FrameRate = 90; %play the movie quickly
    handles.Movie.AcqSamplingFreq = 5; %acquisition sampling frequency
    guidata(hObject, handles);
    %Filter - spatial LPF
    ioi_spatialLPF_call(hObject,handles);
    ioi_play_movie(handles);
    %restore desired movie frequency
    handles.Movie.FrameRate = str2double(get(handles.movie_frequency,'String'));
    %set sliders
    set(handles.message_box,'String','');
    set(handles.message_box,'BackgroundColor','White');
catch
    set(handles.message_box,'String','Movie cannot be opened');
    set(handles.message_box,'BackgroundColor','Red');
end
guidata(hObject, handles);