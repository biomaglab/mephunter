function varargout = MEPHunter_Processing(varargin)
% MEPHunter_Processing M-file for MEPHunter_Processing.fig
%      MEPHunter_Processing, by itself, creates a new MEPHunter_Processing or raises the existing
%      singleton*.
%
%      H = MEPHunter_Processing returns the handle to a new MEPHunter_Processing or the handle to
%      the existing singleton*.
%
%      MEPHunter_Processing('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MEPHunter_Processing.M with the given input arguments.
%
%      MEPHunter_Processing('Property','Value',...) creates a new MEPHunter_Processing or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MEPHunter_Processing_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MEPHunter_Processing_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MEPHunter_Processing

% Last Modified by GUIDE v2.5 26-Apr-2013 16:33:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MEPHunter_Processing_OpeningFcn, ...
                   'gui_OutputFcn',  @MEPHunter_Processing_OutputFcn, ...
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


% --- Executes just before MEPHunter_Processing is made visible.
function MEPHunter_Processing_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MEPHunter_Processing (see VARARGIN)
set(handles.uipanel_ProcessingTrigger,'visible','off')
set(handles.uipanel_ProcessingCorrelation,'visible','off')
set(handles.axes_ProcessingPreview,'visible','off')

handles.cc.data = 0;
handles.cc.filter = 0;

if length (varargin) ~= 0
    handles.emg = varargin{1};
    
    axes(handles.axes_Processing)
    plot (handles.emg.filter);
end

% Choose default command line output for MEPHunter_Processing
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MEPHunter_Processing wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MEPHunter_Processing_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiwait(hObject); 
handles = guidata(hObject);

delete(hObject);

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_ProcessingPath_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ProcessingPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ProcessingPath as text
%        str2double(get(hObject,'String')) returns contents of edit_ProcessingPath as a double


% --- Executes during object creation, after setting all properties.
function edit_ProcessingPath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ProcessingPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in button_ProcessingPath.
function button_ProcessingPath_Callback(hObject, eventdata, handles)
% hObject    handle to button_ProcessingPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fn, pn, index] = uigetfile('*.mat', 'Selecione um arquivo do emg pre-processado');

if index ~= 0
    emg_path = strcat(pn,fn);
    set(handles.edit_ProcessingPath,'String', emg_path);
end

% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in button_ProcessingOpen.
function button_ProcessingOpen_Callback(hObject, eventdata, handles)
% hObject    handle to button_ProcessingOpen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

emg_path = get(handles.edit_ProcessingPath,'String');
loadfile = load (emg_path);
handles.emg = loadfile.emg;
axes(handles.axes_Processing)
plot (handles.emg.filter);

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in button_ProcessingPredictor.
function button_ProcessingPredictor_Callback(hObject, eventdata, handles)
% hObject    handle to button_ProcessingPredictor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.axes_Processing)
[x y] = ginput (2);
x = int32(x);

axes(handles.axes_ProcessingPreview)
handles.predictor = handles.emg.filter(x(1):x(2));
plot(handles.predictor);


% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in button_ProcessingRun.
function button_ProcessingRun_Callback(hObject, eventdata, handles)
% hObject    handle to button_ProcessingRun (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h = waitbar(0,'processando');

delta = length(handles.predictor);
handles.cc.data = zeros (1,length(handles.emg)-delta);


maxemg = find(handles.emg.filter);
diffmax = diff(maxemg);

posemg = find (diffmax>500);
posemg_final = maxemg(posemg);
posemg_inicial = maxemg(posemg+1);
handles.posemg_inicial = [maxemg(1); posemg_inicial];
handles.posemg_final = [posemg_final; maxemg(length(maxemg))];


for j = 1:length(handles.posemg_inicial)
      
    for i = handles.posemg_inicial(j)-50:handles.posemg_final(j)+50

        r = corrcoef(handles.predictor,handles.emg.filter(i:i+delta-1));
        handles.cc.data(i)=r(1,2);      
      
    end
        
    cont =  double(j)/double(length(handles.posemg_inicial));
    waitbar(cont);
   
    
end
close (h)

vetnan = isnan(handles.cc.data);
posnan = find(vetnan==1);
handles.cc.data(posnan)=0;

axes(handles.axes_Processing)
plot (handles.cc.data);
handles.cc.filter = handles.cc.data;


% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in button_ProcessingMax.
function button_ProcessingMax_Callback(hObject, eventdata, handles)
% hObject    handle to button_ProcessingMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


[X Y] = ginput (1);

threshold = find(handles.cc.filter < Y);
handles.cc.filter(threshold) = 0;
axes(handles.axes_Processing)
plot (handles.cc.filter);


maxlocal = find(handles.cc.filter > Y);

diffmax = diff(maxlocal);

pos = find (diffmax>5);
pos = [0,pos,length(maxlocal)];


for i = 1:length(pos)-1
    high{i} = handles.cc.filter(maxlocal(pos(i)+1:pos(i+1)));
    [ymax, xmax] = max(high{i});
    xmax_vet(i) = maxlocal(pos(i)+1)+xmax-1;
    ymax_vet(i) = ymax;
end

axes(handles.axes_Processing)
plot (handles.cc.filter,'b.');
hold on
plot (xmax_vet, ymax_vet,'ro')
hold off

figure
hold on
cores = rand(length(high),3);
for i = 1:length(high)
    MEPs{i} = handles.emg.data(xmax_vet(i):xmax_vet(i)+length(handles.predictor));
    plot(MEPs{i},'-','color',cores(i,:));
    
    pp(i) = max(MEPs{i}) - min(MEPs{i});
end

hold off

handles.emg.MEPs = MEPs;
handles.emg.pp = pp;
handles.emg.pp_mean = mean(pp);
handles.emg.std = std(pp);


% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in button_ProcessingLoad.
function button_ProcessingLoad_Callback(hObject, eventdata, handles)
% hObject    handle to button_ProcessingLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fn, pn, index] = uigetfile('*.mat', 'Selecione um arquivo do emg de Correlação Cruzada');

if index ~= 0
emg_path = strcat(pn,fn);
end

loadfile = load (emg_path);
handles.cc = loadfile.cc;


axes(handles.axes_Processing) 
plot (handles.cc.data);

% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in button_ProcessingAmplitude.
function button_ProcessingAmplitude_Callback(hObject, eventdata, handles)
% hObject    handle to button_ProcessingAmplitude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[x y] = ginput (2);

handles.threshold_amp = find(handles.emg.filter>y(2) & handles.emg.filter<y(1));
filter_aux = handles.emg.filter;
handles.emg.filter(handles.threshold_amp) = 0;


maxemg = find(handles.emg.filter);
diffmax = diff(maxemg);

posemg = find (diffmax>500);
posemg_final = maxemg(posemg);
posemg_inicial = maxemg(posemg+1);
handles.posemg_inicial = [maxemg(1); posemg_inicial];
handles.posemg_final = [posemg_final; maxemg(length(maxemg))];

for i = 1:length(handles.posemg_final);
    handles.emg.filter(handles.posemg_inicial(i)-0:handles.posemg_final(i)+0) = filter_aux(handles.posemg_inicial(i)-0:handles.posemg_final(i)+0);
end

figure
hold on
cores = rand(length(handles.posemg_final),3);
for i = 1:length(handles.posemg_final)
    
    MEPs{i} = handles.emg.filter(handles.posemg_inicial(i)-0:handles.posemg_final(i)+0);
    plot(MEPs{i},'-','color',cores(i,:));
    
end

handles.emg.MEPs = MEPs;

hold off

% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in button_ProcessingSelect.
function button_ProcessingSelect_Callback(hObject, eventdata, handles)
% hObject    handle to button_ProcessingSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.emg = MEPHunter_Select(handles.emg);

handles.output = handles.emg;

uiresume

% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in button_ProcessingCorrelation.
function button_ProcessingCorrelation_Callback(hObject, eventdata, handles)
% hObject    handle to button_ProcessingCorrelation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.uipanel_ProcessingTrigger,'visible','off')
set(handles.uipanel_ProcessingCorrelation,'visible','on')
set(handles.axes_ProcessingPreview,'visible','on')



% --- Executes on button press in button_ProcessingSave.
function button_ProcessingSave_Callback(hObject, eventdata, handles)
% hObject    handle to button_ProcessingSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file,path,index] = uiputfile('Correlacao.mat');

if index == 1
    savefile = strcat(path,file);
    cc = handles.cc;
    save(savefile,'cc')
end


% --- Executes on button press in button_ProcessingQuit.
function button_ProcessingQuit_Callback(hObject, eventdata, handles)
% hObject    handle to button_ProcessingQuit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if exist('handles.emg') == 0;
    handles.output = 0;
else
    handles.output = handles.emg;
end

uiresume

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in button_ProcessingTrigger.
function button_ProcessingTrigger_Callback(hObject, eventdata, handles)
% hObject    handle to button_ProcessingTrigger (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.uipanel_ProcessingCorrelation,'visible','off')
set(handles.axes_ProcessingPreview,'visible','off')

set(handles.uipanel_ProcessingTrigger,'visible','on')


% --- Executes on button press in button_ProcessingTriggerUnit.
function button_ProcessingTriggerUnit_Callback(hObject, eventdata, handles)
% hObject    handle to button_ProcessingTriggerUnit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in button_ProcessingTriggerMap.
function button_ProcessingTriggerMap_Callback(hObject, eventdata, handles)
% hObject    handle to button_ProcessingTriggerMap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Map5x5(handles.emg.OutMap)


% --- Executes on button press in button_ProcessingTriggerHDsEMG.
function button_ProcessingTriggerHDsEMG_Callback(hObject, eventdata, handles)
% hObject    handle to button_ProcessingTriggerHDsEMG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

MapHDsEMG(handles.emg.OutMap)
