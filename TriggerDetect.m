function varargout = TriggerDetect(varargin)
% TRIGGERDETECT MATLAB code for TriggerDetect.fig
%      TRIGGERDETECT, by itself, creates a new TRIGGERDETECT or raises the existing
%      singleton*.
%
%      H = TRIGGERDETECT returns the handle to a new TRIGGERDETECT or the handle to
%      the existing singleton*.
%
%      TRIGGERDETECT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TRIGGERDETECT.M with the given input arguments.
%
%      TRIGGERDETECT('Property','Value',...) creates a new TRIGGERDETECT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TriggerDetect_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TriggerDetect_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TriggerDetect

% Last Modified by GUIDE v2.5 11-Mar-2013 17:33:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TriggerDetect_OpeningFcn, ...
                   'gui_OutputFcn',  @TriggerDetect_OutputFcn, ...
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


% --- Executes just before TriggerDetect is made visible.
function TriggerDetect_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to TriggerDetect (see VARARGIN)

set(handles.radiobutton_Auto,'Value',[1]);
set(handles.radiobutton_Manual,'Value',[0]);


set(handles.button_ManualSelect,'Enable','off')
set(handles.edit_Auto,'Enable','on')
set(handles.button_AutoOK,'Enable','on')



if ~isempty (varargin)
    handles.data = varargin{1};
    handles.trigger = varargin{2};
    handles.xs = varargin{3};
    handles.pos = varargin{4};
    
    axes(handles.axes1)
    hold on
    handles.hdata = plot(handles.xs,handles.data.signal.emg_map{handles.pos}(:,handles.data.channel));
    
     
    
    if ~isempty(handles.trigger)
        axes(handles.axes1)
        handles.htrigger = plot(handles.xs(handles.trigger(:,1)), handles.trigger(:,2),'ro');
    end
end

% if length (varargin) == 2
%     handles.trigger = varargin{2};
%     
%     if length(handles.trigger) ~= 0
%         axes(handles.axes1)
%         handles.htrigger = plot (handles.trigger(:,1), handles.trigger(:,2),'ro');
%     end
%     
% end


% Choose default command line output for TriggerDetect
handles.output_trigger = hObject;
handles.output_data = handles.data;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes TriggerDetect wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = TriggerDetect_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiwait(hObject);
handles = guidata(hObject);

delete(hObject);

% Get default command line output from handles structure
varargout{1} = handles.output_trigger;
varargout{2} = handles.output_data;


function edit_Auto_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Auto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Auto as text
%        str2double(get(hObject,'String')) returns contents of edit_Auto as a double


% --- Executes during object creation, after setting all properties.
function edit_Auto_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Auto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in button_AutoOK.
function button_AutoOK_Callback(hObject, eventdata, handles)
% hObject    handle to button_AutoOK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'htrigger')
    if ishandle(handles.htrigger) == 1
        delete(handles.htrigger)
    end
end

clear trigger handles.htrigger
handles.trigger = BiopacTrigger(handles.data.signal.emg_map{handles.pos}(:,handles.data.channel),...
    str2num(get(handles.edit_Auto,'String')));

if ~isempty(handles.trigger)
    axes(handles.axes1)
    handles.htrigger = plot(handles.xs(handles.trigger(:,1)),handles.trigger(:,2),'ro');
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes TriggerDetect wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Executes on button press in button_ManualSelect.
function button_ManualSelect_Callback(hObject, eventdata, handles)
% hObject    handle to button_ManualSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if isfield(handles,'htrigger')
    if ishandle(handles.htrigger) == 1
        delete(handles.htrigger)
    end
end

clear trigger handles.htrigger



[X y] = getpts(handles.axes1);

r = 1/handles.xs(1);

x = round(X*r);

handles.trigger = [x y];


if ~isempty(handles.trigger)
    axes(handles.axes1)
    handles.htrigger = plot(handles.xs(handles.trigger(:,1)),handles.trigger(:,2),'ro');
end

% Update handles structure
guidata(hObject, handles);




%--------------------------------------------------------------------------

% --- Executes on button press in radiobutton_Auto.
function radiobutton_Auto_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_Auto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_Auto

set(handles.radiobutton_Auto,'Value',[1])
set(handles.radiobutton_Manual,'Value',[0])
set(handles.button_ManualSelect,'Enable','off')
set(handles.edit_Auto,'Enable','on')
set(handles.button_AutoOK,'Enable','on')



% --- Executes on button press in radiobutton_Manual.
function radiobutton_Manual_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_Manual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_Manual

set(handles.radiobutton_Manual,'Value',[1])
set(handles.radiobutton_Auto,'Value',[0])
set(handles.button_ManualSelect,'Enable','on')
set(handles.edit_Auto,'Enable','off')
set(handles.button_AutoOK,'Enable','off')




% --- Executes on button press in button_OK.
function button_OK_Callback(hObject, eventdata, handles)
% hObject    handle to button_OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'trigger') == 0;
    handles.output_trigger = 0;
else
    handles.output_trigger = handles.trigger;
end

uiresume

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in button_outlier.
function button_outlier_Callback(hObject, eventdata, handles)
% hObject    handle to button_outlier (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

x = mod(handles.pos, 13);
if x == 0;
    x = 13;
end
y = ceil(handles.pos/13);

if isfield(handles.data.signal,'emg_data_mono') && handles.data.channel == 1
    handles.data.signal.emg_data_mono = reshape(interpolate_chs_nan(reshape(handles.data.signal.emg_data_mono',...
        13,5,size(handles.data.signal.emg_data_mono,1)),[x y]),65,size(handles.data.signal.emg_data_mono,1))';
    for i = 1:size(handles.data.signal.emg_data_mono, 2)
        handles.data.signal.emg_map{i}(:, handles.data.channel) = handles.data.signal.emg_data_mono(:, i);
    end
else
    handles.data.signal.emg_data_diff = reshape(interpolate_chs_nan(reshape(handles.data.signal.emg_data_diff',...
        13,5,size(handles.data.signal.emg_data_diff,1)),[x y]),65,size(handles.data.signal.emg_data_diff,1))';
    for i = 1:size(handles.data.signal.emg_data_diff, 2)
        handles.data.signal.emg_map{i}(:, handles.data.channel) = handles.data.signal.emg_data_diff(:, i);
    end
end

% finding samples denoting trigger onset
trigger_raw = handles.data.signal.raw_data(:,65) * 5/(2^11);
samples_triggeron = find(trigger_raw>4);
samples_triggeron = [samples_triggeron(diff([-inf;samples_triggeron])>1)];
handles.trigger = [samples_triggeron,...
    handles.data.signal.emg_map{handles.pos}(samples_triggeron, handles.data.channel)];
% handles.htrigger = plot(handles.xs(handles.trigger(:,1)),handles.trigger(:,2),'ro');

if isfield(handles,'trigger') == 0;
    handles.output_trigger = 0;
else
    handles.output_trigger = handles.trigger;
end

handles.output_data = handles.data;

uiresume

guidata(hObject, handles);

%Function to interpolate outliers
function data = interpolate_chs_nan(data,outliers)

for n=1:size(outliers,1)
    data(outliers(n,1),outliers(n,2),:)=nan;
end

nan_mat=nan(size(data,1)+2,size(data,2)+2,size(data,3));
nan_mat(2:size(data,1)+1,2:size(data,2)+1,:)=data;

for n=1:size(outliers,1)
    row=outliers(n,1)+1; col=outliers(n,2)+1;
    data2interp=squeeze([nan_mat(row-1,col-1,:); nan_mat(row-1,col,:); nan_mat(row-1,col+1,:);...
        nan_mat(row,col-1,:); nan_mat(row,col+1,:); nan_mat(row+1,col-1,:); nan_mat(row+1,col,:);...
        nan_mat(row+1,col+1,:)]);
    data(outliers(n,1),outliers(n,2),:)=nanmean(data2interp);
end
