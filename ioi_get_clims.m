function [clims m0 M0] = ioi_get_clims(handles)
Y = handles.Movie.Y;
m0 = min(Y(isfinite(Y(:))));
M0 = max(Y(isfinite(Y(:))));
dM = M0-m0;
pmin = get(handles.edit_pmin,'String');
pmax = get(handles.edit_pmax,'String');
clims = [m0+dM*pmin/100 m0+dM*pmax/100];