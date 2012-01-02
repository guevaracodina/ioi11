function ioi_update_pmax(hObject, handles)
lmin = str2double(get(handles.text_min,'String'));
lmax = str2double(get(handles.text_max,'String'));
max = get(handles.edit_max,'value');
pmax = round(100*(max - lmin)/(lmax-lmin));
set(handles.edit_pmax,'String',pmax);
guidata(hObject, handles);