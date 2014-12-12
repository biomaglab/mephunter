function varargout = MEPDetect_sim(varargin)
% MEPDETECT_SIM MATLAB code for MEPDetect_sim.fig
%      MEPDETECT_SIM, by itself, creates a new MEPDETECT_SIM or raises the existing
%      singleton*.
%
%      H = MEPDETECT_SIM returns the handle to a new MEPDETECT_SIM or the handle to
%      the existing singleton*.
%
%      MEPDETECT_SIM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MEPDETECT_SIM.M with the given input arguments.
%
%      MEPDETECT_SIM('Property','Value',...) creates a new MEPDETECT_SIM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MEPDetect_sim_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MEPDetect_sim_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MEPDetect_sim

% Last Modified by GUIDE v2.5 22-May-2013 17:47:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MEPDetect_sim_OpeningFcn, ...
                   'gui_OutputFcn',  @MEPDetect_sim_OutputFcn, ...
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


% --- Executes just before MEPDetect_sim is made visible.
function MEPDetect_sim_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MEPDetect_sim (see VARARGIN)

if ~isempty (varargin)
    handles.data = varargin{1};
    handles.signal = varargin{2};
    handles.MEPStart = varargin{3};
    handles.xs = varargin{4};
        
    axes(handles.axes1);
    hold on   
end

% get current values of min max positions and plot them
x_max_init =  handles.data.sim_pos_max/handles.data.configuration.fsample;
x_min_init =  handles.data.sim_pos_min/handles.data.configuration.fsample;

handles.hdata = plot(handles.xs+handles.MEPStart,handles.signal);
if max(handles.data.mep_sim_diff) ~= 0 && min(handles.data.mep_sim_diff) ~= 0
    handles.hmepmax_init = plot(x_max_init + handles.MEPStart, max(handles.data.mep_sim_diff),'+g');
    handles.hmepmin_init = plot(x_min_init + handles.MEPStart, min(handles.data.mep_sim_diff),'+g');
end
% -------------------------------------

% Choose default command line output for MEPDetect
handles.sim_minmax = [x_min_init + handles.MEPStart min(handles.data.mep_sim_diff);...
    x_max_init + handles.MEPStart max(handles.data.mep_sim_diff)];
handles.sim_latency = [];

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MEPDetect_sim wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MEPDetect_sim_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiwait(hObject);
handles = guidata(hObject);

delete(hObject);

% Get default command line output from handles structure
varargout{1} = handles.sim_minmax;
varargout{2} = handles.sim_latency;


% --- Executes on button press in button_initial_values.
function button_initial_values_Callback(hObject, eventdata, handles)
% hObject    handle to button_initial_values (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'hdata')
    if ishandle(handles.hdata)
        delete(handles.hdata)
    end
end
if isfield(handles,'hsim_minmax')
    if ishandle(handles.hsim_minmax)
        delete(handles.hsim_minmax)
    end
end
if isfield(handles,'hmepmax_init') && isfield(handles,'hmepmin_init')
    if ishandle(handles.hmepmax_init) && ishandle(handles.hmepmin_init)
        delete(handles.hmepmax_init)
        delete(handles.hmepmin_init)
    end
end

handles.signal = handles.data.mep_sim_diff;
handles.data.sim_pos_max = find(handles.signal == max(handles.signal));
handles.data.sim_pos_min = find(handles.signal == min(handles.signal));

x_max = handles.data.sim_pos_max/handles.data.configuration.fsample;
x_min = handles.data.sim_pos_min/handles.data.configuration.fsample;

handles.sim_minmax = [x_min + handles.MEPStart min(handles.signal);...
    x_max + handles.MEPStart max(handles.signal)];

handles.hdata = plot(handles.xs+handles.MEPStart,handles.signal);
handles.hmepmax = plot(x_max+handles.MEPStart, max(handles.signal),'+g');
handles.hmepmin = plot(x_min+handles.MEPStart, min(handles.signal),'+g');

handles.latency = [];

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in button_ok.
function button_ok_Callback(hObject, eventdata, handles)
% hObject    handle to button_ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles,'sim_minmax');
    handles.sim_minmax = [];
end

if ~isfield(handles,'sim_latency');
    handles.sim_latency = [];
end

uiresume

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in button_latency.
function button_latency_Callback(hObject, eventdata, handles)
% hObject    handle to button_latency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[xlat ylat] = getpts(handles.axes1);
handles.sim_latency = [xlat ylat];


if ~isempty(handles.sim_latency)
    axes(handles.axes1);
    handles.hlatencystart = plot(handles.sim_latency(1,1),handles.sim_latency(1,2),'gv');
end

if length(handles.sim_latency) > 1
    axes(handles.axes1);
    handles.hlatencystop = plot(handles.sim_latency(2,1),handles.sim_latency(2,2),'r^');
end

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in button_manual_select.
function button_manual_select_Callback(hObject, eventdata, handles)
% hObject    handle to button_manual_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'sim_minmax')
    if ishandle(handles.sim_minmax)
        delete(handles.sim_minmax)
    end
end
% trick to shift the x coordinate using MEPStart
% necessery to combine the plots on this MEP window and the Main one
[x y] = getpts(handles.axes1);
sim_minmax = [x y];

if ~isempty(sim_minmax)
    axes(handles.axes1);
    handles.hsim_minmax = plot(sim_minmax(:,1),sim_minmax(:,2),'ro');
end
handles.sim_minmax = [x y];

% Update min and max selected
handles.data.sim_pos_min = round((handles.sim_minmax(1,1)-handles.MEPStart)*handles.data.configuration.fsample);
handles.data.sim_pos_max = round((handles.sim_minmax(2,1)-handles.MEPStart)*handles.data.configuration.fsample);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in button_none.
function button_none_Callback(hObject, eventdata, handles)
% hObject    handle to button_none (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'hmepmax_init') && isfield(handles,'hmepmin_init')
    if ishandle(handles.hmepmax_init) && ishandle(handles.hmepmin_init)
        delete(handles.hmepmax_init)
        delete(handles.hmepmin_init)
    end
end

handles.sim_minmax = [0 0; 0 0];

handles.data.sim_pos_max = [];
handles.data.sim_pos_min = [];
handles.sim_latency = [];

% Update handles structure
guidata(hObject, handles);
