function ioi_update_min(hObject, handles)
lmin = str2double(get(handles.text_min,'String'));
lmax = str2double(get(handles.text_max,'String'));
pmin = get(handles.edit_pmin,'value');
min = lmin + (lmax-lmin)*pmin/100;
set(handles.edit_min,'String',min);
guidata(hObject, handles);