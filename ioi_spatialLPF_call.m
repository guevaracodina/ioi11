function ioi_spatialLPF_call(hObject,handles)
%first pass, Y0 does not exist
if isfield(handles.Movie,'Y0')
    Y0 = handles.Movie.Y0;
else
    Y0 = handles.Movie.Y;
    handles.Movie.Y0 = Y0;
end
Y = Y0;
[nx ny nF] = size(Y0);
radius = round(str2double(get(handles.edit_spatialLPF,'String')));
if radius > 0
    K.k1 = nx;
    K.k2 = ny;
    K.radius = radius;
    K = ioi_spatial_LPF('set',K);
    %Gaussian spatial low pass filter
    for i1=1:nF
        Y(:,:,i1) = ioi_spatial_LPF('lpf',K,squeeze(Y0(:,:,i1)));
    end
end
handles.Movie.Y = Y;
guidata(hObject, handles);
