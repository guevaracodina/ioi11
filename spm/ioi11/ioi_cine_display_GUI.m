function varargout = ioi_cine_display_GUI(varargin)
% IOI_CINE_DISPLAY_GUI MATLAB code for ioi_cine_display_GUI.fig
% Last Modified by GUIDE v2.5 26-Dec-2011 13:11:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ioi_cine_display_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @ioi_cine_display_GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before ioi_cine_display_GUI is made visible.
function ioi_cine_display_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ioi_cine_display_GUI (see VARARGIN)

% Choose default command line output for ioi_cine_display_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = ioi_cine_display_GUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
axes(handles.axes1);
cla;

popup_sel_index = get(handles.popupmenu1, 'Value');
switch popup_sel_index
    case 1
        %Play movie
        ioi_play_movie(handles);
end


% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
[file pathname] = uigetfile({'*.avi;*.mat','Select movies (.avi or .mat format)'});
if ~isequal(file, 0)
    file = fullfile(pathname,file);
    [dir0 fil0 ext0] = fileparts(file);
    switch ext0
        case '.avi'
            obj = VideoReader(file);
            nF = obj.NumberOfFrames;
            % Preallocate movie structure.
            F(1:nF) = struct('cdata',...
                zeros(obj.Height,obj.Width,3,'uint8'),'colormap',[]);
            % Read one frame at a time.
            for k=1:nF
                F(k).cdata = read(obj,k);
            end
        case '.mat'
            open(file);
    end
    %try data = guidata(handles.output); end
    %Use data to store everything
    handles.Movie.obj = obj;
    handles.Movie.F = F;
    handles.Movie.FrameRate = str2double(get(handles.movie_frequency,'String'));
    guidata(hObject, handles);
    %Play movie
    %ioi_play_movie(handles);
end


% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end
delete(handles.figure1)


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
end
set(hObject, 'String', {'Movie'});


% --- Executes on button press in select_subject.
function select_subject_Callback(hObject, eventdata, handles)
[t sts] = spm_select(Inf,'dir','select the top level res directories for the subjects');
if sts && ~isempty(t)
    handles.Info.subjects = t;
    handles.Info.subject_selected = 1; %?
    %write list to GUI
    set(handles.subject_list,'String',t);
    guidata(hObject, handles);
end

% --- Executes on selection change in subject_list.
function subject_list_Callback(hObject, eventdata, handles)
% Hints: contents = cellstr(get(hObject,'String')) returns subject_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from subject_list


% --- Executes during object creation, after setting all properties.
function subject_list_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function movie_frequency_Callback(hObject, eventdata, handles)
handles.Movie.FrameRate =  str2double(get(hObject,'String')); 
%Update
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function movie_frequency_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function play_movie_Callback(hObject, eventdata, handles)
ioi_play_movie(handles);
