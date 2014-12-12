function varargout = MEPHunter_Export(varargin)
% MEPHUNTER_EXPORT M-file for MEPHunter_Export.fig
%      MEPHUNTER_EXPORT, by itself, creates a new MEPHUNTER_EXPORT or raises the existing
%      singleton*.
%
%      H = MEPHUNTER_EXPORT returns the handle to a new MEPHUNTER_EXPORT or the handle to
%      the existing singleton*.
%
%      MEPHUNTER_EXPORT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MEPHUNTER_EXPORT.M with the given input arguments.
%
%      MEPHUNTER_EXPORT('Property','Value',...) creates a new MEPHUNTER_EXPORT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MEPHunter_Export_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MEPHunter_Export_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MEPHunter_Export

% Last Modified by GUIDE v2.5 16-Sep-2011 15:02:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MEPHunter_Export_OpeningFcn, ...
                   'gui_OutputFcn',  @MEPHunter_Export_OutputFcn, ...
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


% --- Executes just before MEPHunter_Export is made visible.
function MEPHunter_Export_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MEPHunter_Export (see VARARGIN)
if length (varargin) ~= 0
    handles.emg = varargin{1};
end

handles.nargin = length (varargin);

% Choose default command line output for MEPHunter_Export
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MEPHunter_Export wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MEPHunter_Export_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% --- Executes on button press in checkbox_ExportMean.
function checkbox_ExportMean_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_ExportMean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_ExportMean


% --- Executes on button press in checkbox_ExportMEPs.
function checkbox_ExportMEPs_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_ExportMEPs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_ExportMEPs


% --- Executes on button press in checkbox_ExportAmplitude.
function checkbox_ExportAmplitude_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_ExportAmplitude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_ExportAmplitude


% --- Executes on button press in checkbox_ExportData.
function checkbox_ExportData_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_ExportData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_ExportData



function edit_ExportMean_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ExportMean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ExportMean as text
%        str2double(get(hObject,'String')) returns contents of edit_ExportMean as a double


% --- Executes during object creation, after setting all properties.
function edit_ExportMean_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ExportMean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ExportAmplitude_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ExportAmplitude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ExportAmplitude as text
%        str2double(get(hObject,'String')) returns contents of edit_ExportAmplitude as a double


% --- Executes during object creation, after setting all properties.
function edit_ExportAmplitude_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ExportAmplitude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ExportMEPs_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ExportMEPs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ExportMEPs as text
%        str2double(get(hObject,'String')) returns contents of edit_ExportMEPs as a double


% --- Executes during object creation, after setting all properties.
function edit_ExportMEPs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ExportMEPs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ExportData_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ExportData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ExportData as text
%        str2double(get(hObject,'String')) returns contents of edit_ExportData as a double


% --- Executes during object creation, after setting all properties.
function edit_ExportData_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ExportData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in button_Export.
function button_Export_Callback(hObject, eventdata, handles)
% hObject    handle to button_Export (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.nargin == 0
    MEPHunter_ErrorExport
else
    
    export_dir = uigetdir;
    
    if get(handles.checkbox_ExportMean,'Value') == 1
        name = get(handles.edit_ExportMean,'String');
        export_name = strcat(export_dir,'/', name);
        savefile1 = [handles.emg.sum_pp, length(handles.emg.pp)];
        savefile2 = [handles.emg.mean_pp, handles.emg.std];
        savefile = [savefile1;savefile2];
        save(export_name,'savefile','-ASCII')
    end
    
    if get(handles.checkbox_ExportAmplitude,'Value') == 1
        name = get(handles.edit_ExportAmplitude,'String');
        export_name = strcat(export_dir,'/', name);
        savefile = handles.emg.pp';
        save(export_name,'savefile','-ASCII')
    end
    
    if get(handles.checkbox_ExportMEPs,'Value') == 1
        
        name = get(handles.edit_ExportMEPs,'String');
        export_name = strcat(export_dir,'/', name);
        
        for i = 1:length(handles.emg.MEPs)
            for j = 1:length(handles.emg.MEPs{i})
                savefile(j,i)=handles.emg.MEPs{i}(j);
            end
        end
        
        save(export_name,'savefile','-ASCII')
    end
    
    if get(handles.checkbox_ExportData,'Value') == 1
        name = get(handles.edit_ExportData,'String');
        export_name = strcat(export_dir,'/', name);
        savefile = [handles.emg.data, handles.emg.filter];
        save(export_name,'savefile','-ASCII')
    end
    
    close
end


% --- Executes on button press in button_ExportQuit.
function button_ExportQuit_Callback(hObject, eventdata, handles)
% hObject    handle to button_ExportQuit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close
