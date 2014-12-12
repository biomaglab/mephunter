function varargout = MEPHunter_Main(varargin)
% MEPHUNTER_MAIN M-file for MEPHunter_Main.fig
%      MEPHUNTER_MAIN, by itself, creates a new MEPHUNTER_MAIN or raises the existing
%      singleton*.
%
%      H = MEPHUNTER_MAIN returns the handle to a new MEPHUNTER_MAIN or the handle to
%      the existing singleton*.
%
%      MEPHUNTER_MAIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MEPHUNTER_MAIN.M with the given input arguments.
%
%      MEPHUNTER_MAIN('Property','Value',...) creates a new MEPHUNTER_MAIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MEPHunter_Main_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MEPHunter_Main_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MEPHunter_Main

% Last Modified by GUIDE v2.5 18-Oct-2011 15:05:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MEPHunter_Main_OpeningFcn, ...
                   'gui_OutputFcn',  @MEPHunter_Main_OutputFcn, ...
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


% --- Executes just before MEPHunter_Main is made visible.
function MEPHunter_Main_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MEPHunter_Main (see VARARGIN)

handles.emg.data = 0;
handles.emg.filter = handles.emg.data;
handles.emg.MEPs = {};
handles.emg.pp = 0;
handles.emg.mean_pp = 0;
handles.emg.std = 0;


% Choose default command line output for MEPHunter_Main
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MEPHunter_Main wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MEPHunter_Main_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_MainPath_Callback(hObject, eventdata, handles)
% hObject    handle to edit_MainPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_MainPath as text
%        str2double(get(hObject,'String')) returns contents of edit_MainPath as a double


% --- Executes during object creation, after setting all properties.
function edit_MainPath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_MainPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in button_MainPath.
function button_MainPath_Callback(hObject, eventdata, handles)
% hObject    handle to button_MainPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fn, pn, index] = uigetfile({'*.txt';'*.mat';'*.dat'}, 'Selecione um arquivo do emg');

if index ~= 0
    emg_path = strcat(pn,fn);
    set(handles.edit_MainPath,'String', emg_path); 
end

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in button_MainOpen.
function button_MainOpen_Callback(hObject, eventdata, handles)
% hObject    handle to button_MainOpen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
EMGstring = get(handles.popup_MainOpen,'String');
EMGvalue = get(handles.popup_MainOpen,'Value');
EMGname = lower(EMGstring(EMGvalue));

switch lower(EMGname{1})
    

    case 'myosystem'    
    emg_path = get(handles.edit_MainPath,'String');
    tam = length(emg_path);
    ext(1:3) = emg_path(tam-2:tam);
    
    switch lower(ext)
        case 'mat'
            loadfile = load (emg_path);
            handles.emg = loadfile.emg;
            axes(handles.axes_Main)
            plot (handles.emg.data);
            
        case {'txt','dat'}
            fid = fopen(emg_path);
            aux_emg = textscan(fid, '%f', 'headerlines', 130);
            handles.emg.data = aux_emg{1};
            handles.emg.filter = handles.emg.data;
            handles.emg.MEPs = {};
            fclose(fid);
            clear fid
            
            axes(handles.axes_Main)
            plot (handles.emg.data);
            
        otherwise
            MEPHunter_ErrorOpen()
    end
    
    case 'biopac'
        figure
        global Sinal
        BioMecanica
        uiwait
        handles.emg.data = Sinal.Dado(:,1);
        handles.emg.filter = handles.emg.data;
        handles.emg.MEPs = {};
        axes(handles.axes_Main)
        plot (handles.emg.data);
        
    case 'biopac map'
        Output = BiopacMap_BIN;
        aux = Output.Sinal.Mapa{1};
        handles.emg.data = aux(:,1);
        handles.emg.filter = handles.emg.data;
        handles.emg.MEPs = {};
        handles.emg.OutMap = Output;
        axes(handles.axes_Main)
        plot (handles.emg.data);
        
    case 'hdsemg'
        output = HDsEMG_BIN;
        aux = output.signal.emg_map{5};
        handles.emg.data = aux;
        handles.emg.filter = handles.emg.data;
        handles.emg.MEPs = {};
        handles.emg.OutMap = output;
        axes(handles.axes_Main)
        plot (handles.emg.data);
        MapHDsEMG(handles.emg.OutMap)
        
    case 'hdsemg_angles'
%         output = HDsEMG_diff_angles_BIN;
% %         aux = output.signal.emg_map{1}{5};
% %         aux = output.data_meps.mepmean{1}{5};
% %         handles.emg.data = aux;
% %         handles.emg.filter = handles.emg.data;
%         handles.emg.MEPs = {};
%         handles.emg.OutMap = output;
%         axes(handles.axes_Main)
% %         plot (handles.emg.data);
%         MapHDsEMG(handles.emg.OutMap)
        HDsEMG_diff_angles_BIN;
        
    case 'hdsemg_diff'
        output = HDsEMG_diff_BIN;
        aux = output.signal.emg_map{5};
        handles.emg.data = aux;
        handles.emg.filter = handles.emg.data;
        handles.emg.MEPs = {};
        handles.emg.OutMap = output;
        axes(handles.axes_Main)
        plot (handles.emg.data);
        MapHDsEMG(handles.emg.OutMap)
        
end
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in button_MainLine.
function button_MainLine_Callback(hObject, eventdata, handles)
% hObject    handle to button_MainLine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes_Main) 
plot (handles.emg.filter,'b-');


% --- Executes on button press in button_MainDot.
function button_MainDot_Callback(hObject, eventdata, handles)
% hObject    handle to button_MainDot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes_Main) 
plot (handles.emg.filter,'b.');


% --- Executes on button press in button_MainViewer.
function button_MainViewer_Callback(hObject, eventdata, handles)
% hObject    handle to button_MainViewer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
MEPHunter_Viewer(handles.emg);



% --- Executes on button press in button_MainFilter.
function button_MainFilter_Callback(hObject, eventdata, handles)
% hObject    handle to button_MainFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.emg = MEPHunter_Filter(handles.emg);

axes(handles.axes_Main) 
plot (handles.emg.filter);

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in button_MainProcessing.
function button_MainProcessing_Callback(hObject, eventdata, handles)
% hObject    handle to button_MainProcessing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% handles.emg = MEPHunter_Processing(handles.emg);
MapHDsEMG_mono_angles;

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in button_MainExport.
function button_MainExport_Callback(hObject, eventdata, handles)
% hObject    handle to button_MainExport (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
MEPHunter_Export(handles.emg);


% --- Executes on selection change in popup_MainOpen.
function popup_MainOpen_Callback(hObject, eventdata, handles)
% hObject    handle to popup_MainOpen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_MainOpen contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_MainOpen


% --- Executes during object creation, after setting all properties.
function popup_MainOpen_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_MainOpen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
