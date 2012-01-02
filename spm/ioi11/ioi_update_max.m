function ioi_update_max(hObject, handles)
lmin = str2double(get(handles.text_min,'String'));
lmax = str2double(get(handles.text_max,'String'));
pmax = str2double(get(handles.edit_pmax,'String'));
max = lmin + (lmax-lmin)*pmax/100;
set(handles.edit_max,'String',max);
guidata(hObject, handles);