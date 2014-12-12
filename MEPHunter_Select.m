function varargout = MEPHunter_Select(varargin)
% MEPHUNTER_SELECT M-file for MEPHunter_Select.fig
%      MEPHUNTER_SELECT, by itself, creates a new MEPHUNTER_SELECT or raises the existing
%      singleton*.
%
%      H = MEPHUNTER_SELECT returns the handle to a new MEPHUNTER_SELECT or the handle to
%      the existing singleton*.
%
%      MEPHUNTER_SELECT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MEPHUNTER_SELECT.M with the given input arguments.
%
%      MEPHUNTER_SELECT('Property','Value',...) creates a new MEPHUNTER_SELECT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MEPHunter_Select_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MEPHunter_Select_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MEPHunter_Select

% Last Modified by GUIDE v2.5 16-Sep-2011 15:08:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MEPHunter_Select_OpeningFcn, ...
                   'gui_OutputFcn',  @MEPHunter_Select_OutputFcn, ...
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


% --- Executes just before MEPHunter_Select is made visible.
function MEPHunter_Select_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MEPHunter_Select (see VARARGIN)
if length (varargin) ~= 0
    
    handles.emg = varargin{1};
    handles.cont = 1;
    handles.pos = 1:length(handles.emg.MEPs);
    
    axes(handles.axes_Select)
    hold on
    
    for i = 1:length(handles.emg.MEPs)
        
        plot(handles.emg.MEPs{i},'-','color',[200 255 170]/255);
        
        pp(i) = max(handles.emg.MEPs{i}) - min(handles.emg.MEPs{i});
        
    end
    
    plot(handles.emg.MEPs{handles.cont},'k-');
    
    pos = strcat('MEP','_',num2str(handles.pos(handles.cont)));
    set(handles.text_Select,'String',pos)
    
    handles.emg.pp = pp;
    handles.emg.mean_pp = mean(handles.emg.pp);
    handles.emg.std = std(handles.emg.pp);
    set(handles.text_SelectAmplitudeMeanValue,'String',num2str(handles.emg.mean_pp))
    
    set(handles.text_SelectAmplitudeMEPValue,'String',num2str(handles.emg.pp(handles.cont)))
    
end

% Choose default command line output for MEPHunter_Select
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MEPHunter_Select wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MEPHunter_Select_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiwait(hObject);
handles = guidata(hObject);

delete(hObject);

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_SelectPath_Callback(hObject, eventdata, handles)
% hObject    handle to edit_SelectPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_SelectPath as text
%        str2double(get(hObject,'String')) returns contents of edit_SelectPath as a double


% --- Executes during object creation, after setting all properties.
function edit_SelectPath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_SelectPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in button_SelectPath.
function button_SelectPath_Callback(hObject, eventdata, handles)
% hObject    handle to button_SelectPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fn, pn, index] = uigetfile('*.mat', 'Selecione um arquivo do emg pre-processado');

if index ~= 0
emg_path = strcat(pn,fn);
set(handles.edit_SelectPath,'String', emg_path);
end

handles.emg_path = strcat(pn,fn);

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in button_SelectOpen.
function button_SelectOpen_Callback(hObject, eventdata, handles)
% hObject    handle to button_SelectOpen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
emg_path = get(handles.edit_SelectPath,'String');
loadfile = load (emg_path);
handles.emg = loadfile.emg;

handles.cont = 1;
handles.pos = 1:length(handles.emg.MEPs);
axes(handles.axes_Select)
hold off
plot(handles.emg.MEPs{1},'k-');
hold on

for i = 1:length(handles.emg.MEPs)
    
    plot(handles.emg.MEPs{i},'-','color',[200 255 170]/255);  
   
end

plot(handles.emg.MEPs{handles.cont},'k-');

pos = strcat('MEP','_',num2str(handles.pos(handles.cont)));
set(handles.text_Select,'String',pos)

set(handles.text_SelectAmplitudeMeanValue,'String',num2str(handles.emg.mean_pp))

set(handles.text_SelectAmplitudeMEPValue,'String',num2str(handles.emg.pp(handles.cont)))

% Update handles structure
guidata(hObject, handles);




% --- Executes on button press in button_SelectNext.
function button_SelectNext_Callback(hObject, eventdata, handles)
% hObject    handle to button_SelectNext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.cont >= length(handles.emg.MEPs)
    handles.cont = 1;
    plot(handles.emg.MEPs{length(handles.emg.MEPs)},'-','color',[200 255 170]/255);
    plot(handles.emg.MEPs{handles.cont},'k-');
    
    pos = strcat('MEP','_',num2str(handles.pos(handles.cont)));
    set(handles.text_Select,'String',pos)
    set(handles.text_SelectAmplitudeMEPValue,'String',num2str(handles.emg.pp(handles.cont)))

else
    handles.cont = handles.cont+1;
    plot(handles.emg.MEPs{handles.cont-1},'-','color',[200 255 170]/255);
    plot(handles.emg.MEPs{handles.cont},'k-');
    
    pos = strcat('MEP','_',num2str(handles.pos(handles.cont)));
    set(handles.text_Select,'String',pos)
    set(handles.text_SelectAmplitudeMEPValue,'String',num2str(handles.emg.pp(handles.cont)))

end


% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in button_SelectPrevious.
function button_SelectPrevious_Callback(hObject, eventdata, handles)
% hObject    handle to button_SelectPrevious (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.cont <= 1
    handles.cont = length(handles.emg.MEPs);
    plot(handles.emg.MEPs{1},'-','color',[200 255 170]/255);
    plot(handles.emg.MEPs{handles.cont},'k-');
    
    pos = strcat('MEP','_',num2str(handles.pos(handles.cont)));
    set(handles.text_Select,'String',pos)
    set(handles.text_SelectAmplitudeMEPValue,'String',num2str(handles.emg.pp(handles.cont)))

else
    handles.cont = handles.cont-1;
    plot(handles.emg.MEPs{handles.cont+1},'-','color',[200 255 170]/255);
    plot(handles.emg.MEPs{handles.cont},'k-');
    
    pos = strcat('MEP','_',num2str(handles.pos(handles.cont)));
    set(handles.text_Select,'String',pos)
    set(handles.text_SelectAmplitudeMEPValue,'String',num2str(handles.emg.pp(handles.cont)))

end


% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in button_SelectReject.
function button_SelectReject_Callback(hObject, eventdata, handles)
% hObject    handle to button_SelectReject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


handles.emg.MEPs(handles.cont) = [];
handles.pos(handles.cont) = [];
handles.emg.pp(handles.cont) = [];


handles.emg.mean_pp = mean(handles.emg.pp);
handles.emg.std = std(handles.emg.pp);


if length(handles.emg.MEPs) == 0;
    uiresume
else
    
    if handles.cont > length(handles.emg.MEPs)
        handles.cont = 1;
        
        hold off
        plot(handles.emg.MEPs{handles.cont},'k-');
        hold on
        
        for i = 1:length(handles.emg.MEPs)
            
            plot(handles.emg.MEPs{i},'-','color',[200 255 170]/255);
            
        end
        plot(handles.emg.MEPs{handles.cont},'k-');
        
        pos = strcat('MEP','_',num2str(handles.pos(handles.cont)));
        set(handles.text_Select,'String',pos)
        set(handles.text_SelectAmplitudeMeanValue,'String',num2str(handles.emg.mean_pp))
        set(handles.text_SelectAmplitudeMEPValue,'String',num2str(handles.emg.pp(handles.cont)))
        
    else
        %handles.cont = handles.cont+1;
        
        hold off
        plot(handles.emg.MEPs{handles.cont},'k-');
        hold on
        
        for i = 1:length(handles.emg.MEPs)
            
            plot(handles.emg.MEPs{i},'-','color',[200 255 170]/255);
            
        end
        plot(handles.emg.MEPs{handles.cont},'k-');
        
        pos = strcat('MEP','_',num2str(handles.pos(handles.cont)));
        set(handles.text_Select,'String',pos)
        set(handles.text_SelectAmplitudeMeanValue,'String',num2str(handles.emg.mean_pp))
        set(handles.text_SelectAmplitudeMEPValue,'String',num2str(handles.emg.pp(handles.cont)))
    end
    
end

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in button_SelectQuit.
function button_SelectQuit_Callback(hObject, eventdata, handles)
% hObject    handle to button_SelectQuit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.emg.sum_pp = sum(handles.emg.pp);

[file,path,index] = uiputfile('Filter.mat');

if index == 1
    savefile = strcat(path,file);
    emg = handles.emg;
    save(savefile,'emg')
end

handles.output = handles.emg;

uiresume

% Update handles structure
guidata(hObject, handles);
