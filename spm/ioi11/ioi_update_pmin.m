function ioi_update_pmin(hObject, handles)
lmin = str2double(get(handles.text_min,'String'));
lmax = str2double(get(handles.text_max,'String'));
min = get(handles.edit_min,'value');
pmin = round(100*(min - lmin)/(lmax-lmin));
set(handles.edit_pmin,'String',pmin);
guidata(hObject, handles);