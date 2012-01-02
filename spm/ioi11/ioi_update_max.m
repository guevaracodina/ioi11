function ioi_update_max(hObject, handles)
lmin = str2double(get(handles.text_min,'String'));
lmax = str2double(get(handles.text_max,'String'));
pmax = get(handles.edit_pmax,'value');
max = lmin + (lmax-lmin)*pmax/100;
set(handles.edit_max,'String',max);
guidata(hObject, handles);