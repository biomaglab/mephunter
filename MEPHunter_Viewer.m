function varargout = MEPHunter_Viewer(varargin)
% MEPHUNTER_VIEWER M-file for MEPHunter_Viewer.fig
%      MEPHUNTER_VIEWER, by itself, creates a new MEPHUNTER_VIEWER or raises the existing
%      singleton*.
%
%      H = MEPHUNTER_VIEWER returns the handle to a new MEPHUNTER_VIEWER or the handle to
%      the existing singleton*.
%
%      MEPHUNTER_VIEWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MEPHUNTER_VIEWER.M with the given input arguments.
%
%      MEPHUNTER_VIEWER('Property','Value',...) creates a new MEPHUNTER_VIEWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MEPHunter_Viewer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MEPHunter_Viewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MEPHunter_Viewer

% Last Modified by GUIDE v2.5 31-Aug-2011 14:59:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MEPHunter_Viewer_OpeningFcn, ...
                   'gui_OutputFcn',  @MEPHunter_Viewer_OutputFcn, ...
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


% --- Executes just before MEPHunter_Viewer is made visible.
function MEPHunter_Viewer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MEPHunter_Viewer (see VARARGIN)
if length (varargin) ~= 0
    
    handles.emg = varargin{1};
    
    if length(handles.emg.MEPs) ~= 0
        
        handles.cont = 1;
        handles.pos = 1:length(handles.emg.MEPs);
        
        axes(handles.axes_Viewer)
        hold on
        
        for i = 1:length(handles.emg.MEPs)
            
            plot(handles.emg.MEPs{i},'-','color',[200 255 170]/255);
            
        end
        
        plot(handles.emg.MEPs{handles.cont},'k-');
        
        pos = strcat('MEP','_',num2str(handles.pos(handles.cont)));
        set(handles.text_Viewer,'String',pos)
        
        set(handles.text_ViewerAmplitudeMeanValue,'String',num2str(handles.emg.mean_pp))
        
        set(handles.text_ViewerAmplitudeMEPValue,'String',num2str(handles.emg.pp(handles.cont)))
        
    end
end

% Choose default command line output for MEPHunter_Viewer
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MEPHunter_Viewer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MEPHunter_Viewer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_ViewerPath_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ViewerPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ViewerPath as text
%        str2double(get(hObject,'String')) returns contents of edit_ViewerPath as a double


% --- Executes during object creation, after setting all properties.
function edit_ViewerPath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ViewerPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in button_ViewerPath.
function button_ViewerPath_Callback(hObject, eventdata, handles)
% hObject    handle to button_ViewerPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fn, pn, index] = uigetfile('*.mat', 'Selecione um arquivo do emg');

if index ~= 0
emg_path = strcat(pn,fn);
set(handles.edit_ViewerPath,'String', emg_path);
end

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in button_ViewerOpen.
function button_ViewerOpen_Callback(hObject, eventdata, handles)
% hObject    handle to button_ViewerOpen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
emg_path = get(handles.edit_ViewerPath,'String');
loadfile = load (emg_path);
handles.emg = loadfile.emg;

handles.cont = 1;
handles.pos = 1:length(handles.emg.MEPs);
axes(handles.axes_Viewer)
hold off
plot(handles.emg.MEPs{1},'k-');
hold on

for i = 1:length(handles.emg.MEPs)
    
    plot(handles.emg.MEPs{i},'-','color',[200 255 170]/255);  
   
end

plot(handles.emg.MEPs{handles.cont},'k-');

pos = strcat('MEP','_',num2str(handles.pos(handles.cont)));
set(handles.text_Viewer,'String',pos)

set(handles.text_ViewerAmplitudeMeanValue,'String',num2str(handles.emg.mean_pp))

set(handles.text_ViewerAmplitudeMEPValue,'String',num2str(handles.emg.pp(handles.cont)))

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in button_ViewerNext.
function button_ViewerNext_Callback(hObject, eventdata, handles)
% hObject    handle to button_ViewerNext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.cont >= length(handles.emg.MEPs)
    handles.cont = 1;
    plot(handles.emg.MEPs{length(handles.emg.MEPs)},'-','color',[200 255 170]/255);
    plot(handles.emg.MEPs{handles.cont},'k-');
    
    pos = strcat('MEP','_',num2str(handles.pos(handles.cont)));
    set(handles.text_Viewer,'String',pos)
    set(handles.text_ViewerAmplitudeMEPValue,'String',num2str(handles.emg.pp(handles.cont)))

else
    handles.cont = handles.cont+1;
    plot(handles.emg.MEPs{handles.cont-1},'-','color',[200 255 170]/255);
    plot(handles.emg.MEPs{handles.cont},'k-');
    
    pos = strcat('MEP','_',num2str(handles.pos(handles.cont)));
    set(handles.text_Viewer,'String',pos)
    set(handles.text_ViewerAmplitudeMEPValue,'String',num2str(handles.emg.pp(handles.cont)))

end


% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in button_ViewerPrevious.
function button_ViewerPrevious_Callback(hObject, eventdata, handles)
% hObject    handle to button_ViewerPrevious (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.cont <= 1
    handles.cont = length(handles.emg.MEPs);
    plot(handles.emg.MEPs{1},'-','color',[200 255 170]/255);
    plot(handles.emg.MEPs{handles.cont},'k-');
    
    pos = strcat('MEP','_',num2str(handles.pos(handles.cont)));
    set(handles.text_Viewer,'String',pos)
    set(handles.text_ViewerAmplitudeMEPValue,'String',num2str(handles.emg.pp(handles.cont)))
    
else
    handles.cont = handles.cont-1;
    plot(handles.emg.MEPs{handles.cont+1},'-','color',[200 255 170]/255);
    plot(handles.emg.MEPs{handles.cont},'k-');
    
    pos = strcat('MEP','_',num2str(handles.pos(handles.cont)));
    set(handles.text_Viewer,'String',pos)
    set(handles.text_ViewerAmplitudeMEPValue,'String',num2str(handles.emg.pp(handles.cont)))

end

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in button_ViewerQuit.
function button_ViewerQuit_Callback(hObject, eventdata, handles)
% hObject    handle to button_ViewerQuit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close

