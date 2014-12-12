function varargout = MAPAmplitudes(varargin)
% MAPAMPLITUDES MATLAB code for MAPAmplitudes.fig
%      MAPAMPLITUDES, by itself, creates a new MAPAMPLITUDES or raises the existing
%      singleton*.
%
%      H = MAPAMPLITUDES returns the handle to a new MAPAMPLITUDES or the handle to
%      the existing singleton*.
%
%      MAPAMPLITUDES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAPAMPLITUDES.M with the given input arguments.
%
%      MAPAMPLITUDES('Property','Value',...) creates a new MAPAMPLITUDES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MAPAmplitudes_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MAPAmplitudes_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MAPAmplitudes

% Last Modified by GUIDE v2.5 03-Oct-2013 09:07:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MAPAmplitudes_OpeningFcn, ...
                   'gui_OutputFcn',  @MAPAmplitudes_OutputFcn, ...
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


% --- Executes just before MAPAmplitudes is made visible.
function MAPAmplitudes_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MAPAmplitudes (see VARARGIN)

if ~isempty (varargin)
    handles.data = varargin{1};
    handles.signal = varargin{2};
end

% center the figure window on the screen
movegui(hObject, 'center');

if isfield(handles.data, 'iz_row')
    set(handles.edit_iz_row,'String', num2str(handles.data.iz_row));
else
    handles.data.iz_row = 7;
    set(handles.edit_iz_row,'String', num2str(handles.data.iz_row));
end

handles.amp_map = cell(1, handles.signal.n_conditions);
aux_map_min = zeros(1, handles.signal.n_conditions);
aux_map_max = zeros(1, handles.signal.n_conditions);

for i = 1:handles.signal.n_conditions
    handles.amp_map{i} = handles.data.amp_pp_map{i}(:,:);
    handles.fmed_map{i} = handles.data.amp_fmed_map{i}(:,:);
    aux_map_min(i) = min(min(handles.amp_map{i}));
    aux_map_max(i) = max(max(handles.amp_map{i}));
end

% Plotting meps' RMS image interpolated
handles.map_min_init = min(aux_map_min);
handles.map_max_init = max(aux_map_max);
set(handles.edit_image_color_min, 'String', num2str(handles.map_min_init,'%.1f'))
set(handles.edit_image_color_max, 'String', num2str(handles.map_max_init,'%.1f'))

for i = 1:handles.signal.n_conditions
    axes(eval(strcat('handles.axes',num2str(i))));
    handles.hmap_amp{i} = image(interp2(handles.amp_map{i},8), 'CDataMapping','scaled', 'parent', gca);
    
    %  maps with one scale for all
%     set(gca,'xtick',x_tick, 'ytick',y_tick, 'xticklabel',[], 'yticklabel', [],...
%         'CLimMode', 'manual','CLim', [handles.map_min_init handles.map_max_init])
    
    % maps with individual scale
%     set(gca,'xtick',x_tick, 'ytick',y_tick, 'xticklabel',[], 'yticklabel', [],...
%         'CLimMode', 'manual','CLim', [aux_map_min(i) aux_map_max(i)])
    
    if i == 1 || i == 5
        y_ticklabel = {'1';'2';'3';'4';'5';'6';'7';'8';'9';'10';'11';'12';'13'};
        if i == 5
            x_ticklabel = {'1';'2';'3';'4';'5'};
            set(gca,'xticklabel',x_ticklabel, 'yticklabel', y_ticklabel)
        else
            set(gca,'yticklabel', y_ticklabel)
        end
    end
    if i > 5
        x_ticklabel = {'1';'2';'3';'4';'5'};
        set(gca,'xtick',x_tick, 'xticklabel',x_ticklabel)
    end
    % hide all colorbars
    colorbar('hide', 'peer', gca);
    % show one colorbar for each map
%     colorbar('peer', gca);
    set(eval(strcat('handles.text_angle',num2str(i))), 'String',...
        [num2str(handles.signal.angle_stim(i)) 'º'])
end

handles.img_color_bar = image(interp2(handles.amp_map{1}, 8),...
    'Visible','off', 'CDataMapping','scaled', 'parent', handles.axes9);
% maps with one scale for all
set(handles.axes9,'Visible','off',...
    'CLimMode', 'manual','CLim', [handles.map_min_init handles.map_max_init])
% maps with individual scale
% set(handles.axes9,'Visible','off',...
%     'CLimMode', 'manual','CLim', [aux_map_min(2) aux_map_max(2)])
colorbar('peer', handles.axes9);

handles.output = handles.data.iz_row;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MAPAmplitudes wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MAPAmplitudes_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiwait(hObject);
handles = guidata(hObject);

delete(hObject);

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in button_amplitude_image.
function button_amplitude_image_Callback(hObject, eventdata, handles)
% hObject    handle to button_amplitude_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'hmap_amp')
    for i = 1:handles.signal.n_conditions
        if ishandle(handles.hmap_amp{i})
            delete(handles.hmap_amp{i})
        end
    end
end

if isfield(handles,'img_color_bar')
    if ishandle(handles.img_color_bar)
        delete(handles.img_color_bar)
    end
end

handles = amplitude_map(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in button_default_values.
function button_default_values_Callback(hObject, eventdata, handles)
% hObject    handle to button_default_values (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'hmap_amp')
    for i = 1:handles.signal.n_conditions
        if ishandle(handles.hmap_amp{i})
            delete(handles.hmap_amp{i})
        end
    end
end

if isfield(handles,'img_color_bar')
    if ishandle(handles.img_color_bar)
        delete(handles.img_color_bar)
    end
end

set(handles.edit_image_color_max,'String', num2str(handles.map_max_init,'%.1f'));
set(handles.edit_image_color_min,'String', num2str(handles.map_min_init,'%.1f'));

set(handles.radiobutton_map_rms,'Value', 1);
set(handles.radiobutton_map_pp,'Value', 0);
button_amplitude_image_Callback(hObject, eventdata, handles)

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in button_ok.
function button_ok_Callback(hObject, eventdata, handles)
% hObject    handle to button_ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiresume

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in button_save_figure.
function button_save_figure_Callback(hObject, eventdata, handles)
% hObject    handle to button_save_figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure_name_default = [num2str(handles.signal.patient_id) '_' handles.signal.emg_side];

[figure_name, figure_path] = uiputfile({'*.emf','Enhanced metafile (*.emf)';...
    '*.bmp','Windows bitmap (*.bmp)'},...
    'Save figure as...', figure_name_default);

saveas(handles.figure1, [figure_path figure_name])

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in radiobutton_map_pp.
function radiobutton_map_pp_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_map_pp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_map_pp


if isfield(handles,'hmap_amp')
    for i = 1:handles.signal.n_conditions
        if ishandle(handles.hmap_amp{i})
            delete(handles.hmap_amp{i})
        end
    end
end

if isfield(handles,'img_color_bar')
    if ishandle(handles.img_color_bar)
        delete(handles.img_color_bar)
    end
end

handles.map_pp_value = get(handles.radiobutton_map_pp,'Value');
aux_map_min = zeros(1, handles.signal.n_conditions);
aux_map_max = zeros(1, handles.signal.n_conditions);

if handles.map_pp_value
    set(handles.radiobutton_map_rms,'Value', 0);
    for i = 1:handles.signal.n_conditions
        handles.amp_map{i} = handles.data.amp_pp_map{i}(:,:);
        aux_map_min(i) = min(min(handles.amp_map{i}));
        aux_map_max(i) = max(max(handles.amp_map{i}));
    end
else
    set(handles.radiobutton_map_rms,'Value', 1);
    for i = 1:handles.signal.n_conditions
        handles.amp_map{i} = handles.data.amp_rms_map{i}(:,:);
        aux_map_min(i) = min(min(handles.amp_map{i}));
        aux_map_max(i) = max(max(handles.amp_map{i}));
    end
end

% maps color from all maps
map_min = min(aux_map_min);
map_max = max(aux_map_max);

% maps color from 45 degrees map
% map_min = aux_map_min(2);
% map_max = aux_map_max(2);

set(handles.edit_image_color_min, 'String', num2str(map_min,'%.1f'))
set(handles.edit_image_color_max, 'String', num2str(map_max,'%.1f'))

handles = amplitude_map(handles);

% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in radiobutton_map_rms.
function radiobutton_map_rms_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_map_rms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_map_rms


if isfield(handles,'hmap_amp')
    for i = 1:handles.signal.n_conditions
        if ishandle(handles.hmap_amp{i})
            delete(handles.hmap_amp{i})
        end
    end
end

if isfield(handles,'img_color_bar')
    if ishandle(handles.img_color_bar)
        delete(handles.img_color_bar)
    end
end

handles.map_rms_value = get(handles.radiobutton_map_rms,'Value');
aux_map_min = zeros(1, handles.signal.n_conditions);
aux_map_max = zeros(1, handles.signal.n_conditions);

if handles.map_rms_value
    set(handles.radiobutton_map_pp,'Value', 0);
    for i = 1:handles.signal.n_conditions
        handles.amp_map{i} = handles.data.amp_rms_map{i}(:,:);
        aux_map_min(i) = min(min(handles.amp_map{i}));
        aux_map_max(i) = max(max(handles.amp_map{i}));
    end
else
    set(handles.radiobutton_map_pp,'Value', 1);
    for i = 1:handles.signal.n_conditions
        handles.amp_map{i} = handles.data.amp_pp_map{i}(:,:);
        aux_map_min(i) = min(min(handles.amp_map{i}));
        aux_map_max(i) = max(max(handles.amp_map{i}));
    end
end

map_min = min(aux_map_min);
map_max = max(aux_map_max);
set(handles.edit_image_color_min, 'String', num2str(map_min,'%.1f'))
set(handles.edit_image_color_max, 'String', num2str(map_max,'%.1f'))

handles = amplitude_map(handles);

% Update handles structure
guidata(hObject, handles);


function edit_image_color_max_Callback(hObject, eventdata, handles)
% hObject    handle to edit_image_color_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_image_color_max as text
%        str2double(get(hObject,'String')) returns contents of edit_image_color_max as a double


% --- Executes during object creation, after setting all properties.
function edit_image_color_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_image_color_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_image_color_min_Callback(hObject, eventdata, handles)
% hObject    handle to edit_image_color_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_image_color_min as text
%        str2double(get(hObject,'String')) returns contents of edit_image_color_min as a double


% --- Executes during object creation, after setting all properties.
function edit_image_color_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_image_color_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_iz_row_Callback(hObject, eventdata, handles)
% hObject    handle to edit_iz_row (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_iz_row as text
%        str2double(get(hObject,'String')) returns contents of edit_iz_row as a double

handles.data.iz_row = str2double(get(handles.edit_iz_row,'String'));
handles.output = handles.data.iz_row;

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_iz_row_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_iz_row (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function out = amplitude_map(handles)

color_min = str2double(get(handles.edit_image_color_min, 'String'));
color_max = str2double(get(handles.edit_image_color_max, 'String'));

for i = 1:handles.signal.n_conditions
    axes(eval(strcat('handles.axes',num2str(i))));
    handles.hmap_amp{i} = image(interp2(handles.amp_map{i}, 8), 'CDataMapping','scaled', 'parent', gca);
    x_lims = get(gca, 'xlim');
    x_step = abs((x_lims(2)-x_lims(1))/5);
    x_tick =  [x_lims(1)+x_step/2:x_step:x_lims(2)+x_step/2];
    y_lims = get(gca, 'ylim');
    y_step = abs((y_lims(2)-y_lims(1))/13);
    y_tick =  [y_lims(1)+y_step/2:y_step:y_lims(2)+y_step/2];
    set(gca,'xtick',x_tick, 'ytick',y_tick, 'xticklabel',[], 'yticklabel', [],...
        'CLimMode', 'manual','CLim', [color_min color_max])
    if i == 1 || i == 5
        y_ticklabel = {'1';'2';'3';'4';'5';'6';'7';'8';'9';'10';'11';'12';'13'};
        if i == 5
            x_ticklabel = {'1';'2';'3';'4';'5'};
            set(gca,'xticklabel',x_ticklabel, 'yticklabel', y_ticklabel)
        else
            set(gca,'yticklabel', y_ticklabel)
        end
    end
    if i > 5
        x_ticklabel = {'1';'2';'3';'4';'5'};
        set(gca,'xtick',x_tick, 'xticklabel',x_ticklabel)
    end
    colorbar('hide', 'peer', gca);
end

handles.img_color_bar = image(interp2(handles.amp_map{1}, 8),...
    'Visible','off', 'CDataMapping','scaled', 'parent', handles.axes9);
set(handles.axes9,'Visible','off',...
    'CLimMode', 'manual','CLim', [color_min color_max])
colorbar('peer', handles.axes9);

out = handles;
