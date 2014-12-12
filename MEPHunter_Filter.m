function varargout = MEPHunter_Filter(varargin)
% MEPHUNTER_FILTER M-file for MEPHunter_Filter.fig
%      MEPHUNTER_FILTER, by itself, creates a new MEPHUNTER_FILTER or raises the existing
%      singleton*.
%
%      H = MEPHUNTER_FILTER returns the handle to a new MEPHUNTER_FILTER or the handle to
%      the existing singleton*.
%
%      MEPHUNTER_FILTER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MEPHUNTER_FILTER.M with the given input arguments.
%
%      MEPHUNTER_FILTER('Property','Value',...) creates a new MEPHUNTER_FILTER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MEPHunter_Filter_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MEPHunter_Filter_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MEPHunter_Filter

% Last Modified by GUIDE v2.5 16-Sep-2011 17:22:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MEPHunter_Filter_OpeningFcn, ...
                   'gui_OutputFcn',  @MEPHunter_Filter_OutputFcn, ...
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


% --- Executes just before MEPHunter_Filter is made visible.
function MEPHunter_Filter_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MEPHunter_Filter (see VARARGIN)


handles.emg = varargin{1};
handles.markers = {};

axes(handles.axes_Filter) 
plot (handles.emg.filter);

% Choose default command line output for MEPHunter_Filter
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MEPHunter_Filter wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MEPHunter_Filter_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

uiwait(hObject);
handles = guidata(hObject);

delete(hObject);
% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in button_FilterAmplitude.
function button_FilterAmplitude_Callback(hObject, eventdata, handles)
% hObject    handle to button_FilterAmplitude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[x y] = ginput (2);

handles.limiar_amp = find(handles.emg.filter>y(2) & handles.emg.filter<y(1));
filter_aux = handles.emg.filter;
handles.emg.filter(handles.limiar_amp) = 0;


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

for i = 1:length(handles.markers)
    handles.emg.filter(handles.markers{i}(1):handles.markers{i}(2)) = handles.emg.data(handles.markers{i}(1):handles.markers{i}(2));
end

handles.output = handles.emg;

axes(handles.axes_Filter) 
plot (handles.emg.filter);

hold on
for i = 1:length(handles.markers)
plot (handles.markers{i}(1), handles.emg.filter(handles.markers{i}(1)), 'gv');
plot (handles.markers{i}(2), handles.emg.filter(handles.markers{i}(2)), 'rv');
end
hold off

% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in button_FilterReject.
function button_FilterReject_Callback(hObject, eventdata, handles)
% hObject    handle to button_FilterReject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[x y] = ginput (2);
x = int32(x);

handles.emg.filter(x(1):x(2)) = 0;

for i = 1:length(handles.markers)
    handles.emg.filter(handles.markers{i}(1):handles.markers{i}(2)) = handles.emg.data(handles.markers{i}(1):handles.markers{i}(2));
end
    

handles.output = handles.emg;

axes(handles.axes_Filter) 
plot (handles.emg.filter);

hold on
for i = 1:length(handles.markers)
plot (handles.markers{i}(1), handles.emg.filter(handles.markers{i}(1)), 'gv');
plot (handles.markers{i}(2), handles.emg.filter(handles.markers{i}(2)), 'rv');
end
hold off

% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in button_FilterSave.
function button_FilterSave_Callback(hObject, eventdata, handles)
% hObject    handle to button_FilterSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file,path] = uiputfile('Filter.mat');
savefile = strcat(path,file);
emg = handles.emg;
save(savefile,'emg')



% --- Executes on button press in button_FilterReset.
function button_FilterReset_Callback(hObject, eventdata, handles)
% hObject    handle to button_FilterReset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.emg.filter = handles.emg.data;
handles.markers = {};

axes(handles.axes_Filter) 
plot (handles.emg.filter);

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in button_FilterQuit.
function button_FilterQuit_Callback(hObject, eventdata, handles)
% hObject    handle to button_FilterQuit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output = handles.emg;

uiresume

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in button_FilterMarker.
function button_FilterMarker_Callback(hObject, eventdata, handles)
% hObject    handle to button_FilterMarker (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[x y] = ginput (2);
x = int32(x);

handles.markers{length(handles.markers)+1} = x;

hold on
axes(handles.axes_Filter)

for i = 1:length(handles.markers)
plot (handles.markers{i}(1), handles.emg.filter(handles.markers{i}(1)), 'gv');
plot (handles.markers{i}(2), handles.emg.filter(handles.markers{i}(2)), 'rv');
end


hold off

% Update handles structure
guidata(hObject, handles);
