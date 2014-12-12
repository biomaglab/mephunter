function varargout = MAPVisualization(varargin)
% MAPVISUALIZATION MATLAB code for MAPVisualization.fig
%      MAPVISUALIZATION, by itself, creates a new MAPVISUALIZATION or raises the existing
%      singleton*.
%
%      H = MAPVISUALIZATION returns the handle to a new MAPVISUALIZATION or the handle to
%      the existing singleton*.
%
%      MAPVISUALIZATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAPVISUALIZATION.M with the given input arguments.
%
%      MAPVISUALIZATION('Property','Value',...) creates a new MAPVISUALIZATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MAPVisualization_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MAPVisualization_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MAPVisualization

% Last Modified by GUIDE v2.5 19-Jun-2013 12:38:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MAPVisualization_OpeningFcn, ...
                   'gui_OutputFcn',  @MAPVisualization_OutputFcn, ...
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


% --- Executes just before MAPVisualization is made visible.
function MAPVisualization_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MAPVisualization (see VARARGIN)

if ~isempty (varargin)
    handles.data = varargin{1};
    handles.channel = varargin{2};
end

handles.count = 0;

if size(handles.data.signal.emg_map{1}, 2) == 1
    handles.map_channel = 1;
else
    handles.map_channel = handles.data.channel;
end

handles.amp_map = handles.data.amp_rms_map{handles.map_channel}(:,:);
handles.fmed_map = handles.data.amp_fmed_map{handles.map_channel}(:,:);

% Rearranging variables
averaged_meps = nan(size(handles.data.mepmean{handles.data.channel}{2}, 1), length(handles.data.signal.emg_map));

for i = 1:length(handles.data.signal.emg_map)
    if isempty(handles.data.mepmean{handles.data.channel}{i}) || isempty(handles.data.mepmax{handles.data.channel}{i})
        averaged_meps(:,i) = nan(size(handles.data.mepmean{handles.data.channel}{2}, 1),1);
    else
        averaged_meps(:,i) = handles.data.mepmean{handles.data.channel}{i};
    end
end

% Plotting meps' line
for cols = 1:5
    line([1:size(averaged_meps,1)] + (cols-1)*size(averaged_meps,1),...
        averaged_meps(:,[1:13] + (cols-1)*13)/max(abs(averaged_meps(:))) + repmat(13:-1:1,size(averaged_meps,1),1),...
        'parent', handles.axes1,'color', 'b', 'LineWidth', 2)
end

% Plotting meps' RMS image interpolated
handles.map_min = min(min(handles.amp_map));
handles.map_max = max(max(handles.amp_map));
set(handles.edit_image_color_min, 'String', num2str(handles.map_min,3))
set(handles.edit_image_color_max, 'String', num2str(handles.map_max,3))
amp_map_interp = interp2(handles.amp_map, 8);

handles.hmap_amp = image(amp_map_interp, 'CDataMapping','scaled', 'parent', handles.axes2);
colorbar('peer', handles.axes2);

% Plotting meps' Fmed image interpolated
handles.fmed_map_min = min(min(handles.fmed_map));
handles.fmed_map_max = max(max(handles.fmed_map));
% set(handles.edit_image_color_min, 'String', num2str(handles.map_min,3))
% set(handles.edit_image_color_max, 'String', num2str(handles.map_max,3))
fmed_map_interp = interp2(handles.fmed_map, 8);

handles.hfmed_amp = image(fmed_map_interp, 'CDataMapping','scaled', 'parent', handles.axes4);
colorbar('peer', handles.axes4);

% Setting text and limits values
set(handles.axes1,'xlim',[0 (cols)*size(averaged_meps,1)],'ylim',[0 14],...
    'xtick', size(averaged_meps,1)*([1:cols] - 1/2),'xticklabel',1:cols,...
    'ytick',1:13,'yticklabel', 13:-1:1,'ydir','normal','box','on','xgrid','on','ygrid','on')
set(handles.axes2,'xtick',1:cols,'ytick',1:13, 'ydir','reverse',...
    'CLimMode', 'manual','CLim', [handles.map_min handles.map_max])
set(handles.axes3,'xlim', [0 6],'ylim',[0 14],...
    'xtick',1:cols,'ytick',1:13,'ydir','reverse','box','on',...
    'xgrid','on','ygrid','on')
set(handles.axes4,'xtick',1:cols,'ytick',1:13, 'ydir','reverse',...
    'CLimMode', 'manual','CLim', [handles.fmed_map_min handles.fmed_map_max])

x_label_str = strcat('Columns (', get(handles.edit_ied,'String'), ' mm IED)');
y_label_str = strcat('Rows (', get(handles.edit_ied,'String'), ' mm IED)');

tx(1) = ylabel(handles.axes1, y_label_str, 'fontsize', 11);

set(handles.text_title_line, 'String', sprintf('Raw MEPs (%3.2f microV/Div)',max(abs(averaged_meps(:)))));
set(handles.text_x_label_line, 'String', x_label_str);
set(handles.text_title_map, 'String', 'Map of MEP Amplitude');
set(handles.text_x_label_map, 'String', x_label_str);
set(handles.text_title_fmed_map, 'String', 'Map of Median Frequency');
set(handles.text_x_label_fmed_map, 'String', x_label_str);
set(handles.text_x_label_cluster, 'String', x_label_str);
set(handles.text_title_cluster, 'String', 'Channels detecting greatest MEPs');

% Choose default command line output for MAPVisualization
handles.output = handles.data;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MAPVisualization wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MAPVisualization_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiwait(hObject);
handles = guidata(hObject);

delete(hObject);

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in button_cluster_threshold.
function button_cluster_threshold_Callback(hObject, ~, handles)
% hObject    handle to button_cluster_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'hmarkers')
    if ishandle(handles.hmarkers)
        delete(handles.hmarkers)
        delete(handles.hcog)
    end
end

% Identifying clusters of activity from MEPs' maps
aux_map = handles.amp_map;
aux_map(isnan(aux_map))=0;
Local_MEP = WSCalculation(aux_map,1);

for clusters = 1:max(Local_MEP(:))
    % Mean MEP activity per cluster
    cmean(clusters) = mean(aux_map(Local_MEP==clusters));
    % number of channels per cluster
    cnumberofchannels(clusters) = numel(Local_MEP(Local_MEP==clusters));
    if isnan(cmean(clusters)), cmean(clusters)=0; cnumberofchannels(clusters) = inf; end
end    
cmean(cnumberofchannels<5) = []; cnumberofchannels(cnumberofchannels<5) = [];
% cpos == cluster of channels detecting greatest MEPs
[junk,cpos] = max(cmean);
% indices pointing to channels detecting greatest MEPs
cind = find(Local_MEP == cpos);

threshold = str2double(get(handles.edit_amp_threshold,'String'))/100*max(aux_map(cind));
    
% Identifying channels detecting MEPs over threshold
MEP_channels = aux_map(cind) >= threshold;
MEP_channels = cind(MEP_channels);
[MEP_rows,MEP_columns] = ind2sub(size(aux_map),MEP_channels);
handles.data.mepcluster{handles.data.channel} = [MEP_rows, MEP_columns];

% plotting channels detecting greatest MEPs
axes(handles.axes3);
hold on
% center of gravity
masc = zeros(13,5);
masc(MEP_rows, MEP_columns) = 1;
masc(1, [1 5]) = 0;
masc(13, [1 5]) = 0;

% TODO: it will not work if start mep hunter with hdsemg mono - fix it
if handles.channel == 2 || size(handles.data.signal.emg_map{1}, 2) == 1
    masc(13, :) = 0;
    masc(12, [1 5]) = 0;
end
mapa = mat2gray(handles.amp_map);
mapa_cinza = uint8(round(mapa*255));
mapa_prop = regionprops(masc, mapa_cinza, 'WeightedCentroid');
cog = round(cat(1, mapa_prop.WeightedCentroid));

handles.hmarkers = plot(MEP_columns,MEP_rows,'o', 'markersize', 16,'markerfacecolor','y');
handles.hcog = plot(cog(:, 1),cog(:, 2),'o', 'markersize', 16, 'markerfacecolor', 'g');
hold off

clear cpos threshold aux_map cmean cnumberofchannels

if get(handles.checkbox_map_rms,'Value')
    handles.data.mep_cluster_diff_rms = [MEP_rows, MEP_columns];
    handles.data.size_cluster_diff_rms = size(handles.data.mep_cluster_diff_rms, 1);
    handles.data.cog_cluster_diff_rms = [cog(1,2) cog(1,1)];
    handles.data.cog_cluster_diff_pos_rms = cog(1,2) + (cog(1,1)-1)*13;
else
    handles.data.mep_cluster_diff_pp = [MEP_rows, MEP_columns];
    handles.data.size_cluster_diff_pp = size(handles.data.mep_cluster_diff_pp, 1);
    handles.data.cog_cluster_diff_pp = [cog(1,2) cog(1,1)];
    handles.data.cog_cluster_diff_pos_pp = cog(1,2) + (cog(1,1)-1)*13;
end

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in button_cond_velocity.
function button_cond_velocity_Callback(hObject, eventdata, handles)
% hObject    handle to button_cond_velocity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% CV estimation makes sense only for differential EMGs
% to calculate for each stimulus

axes(handles.axes1);
hold on

handles.count = handles.count + 1;
IED = str2double(get(handles.edit_ied,'String'))/1000;
[x y] = getpts(handles.axes1);
aux_cv_channels = [x y];

if ~isempty(aux_cv_channels)
    
    % Storage of cv_channels and for plotting, were made some calculus
    % to mirror coordinates and fit on axis 1 graph
    aux_cv_channels = [round(1/2+aux_cv_channels(:,1)/size(handles.data.mepmean{handles.channel}{2}, 1)) (round(aux_cv_channels(:,2))*(-1)+14)];
    handles.hcv_channels{handles.count} = plot(((aux_cv_channels(:,1)-1/2)*size(handles.data.mepmean{handles.channel}{2}, 1)),...
        (aux_cv_channels(:,2)*(-1)+14),'ro');
    handles.cv_channels = [13*(aux_cv_channels(:, 1)-1)+(aux_cv_channels(:, 2))]';
    
    % Calculus of conduction velocity and correlation coeficient for the
    % mean signal of the ten stimulus
    signal = zeros(length(handles.cv_channels), size(handles.data.mepmean{handles.data.channel}{2}, 1));
    for i = handles.cv_channels
        signal(i-(handles.cv_channels(1)-1), :) = handles.data.mepmean{handles.data.channel}{i}';
    end
    [aux_cvelocity, aux_correl_coef] = CVCalculation(signal,IED,handles.data.configuration.fsample);
else
    aux_cvelocity = nan;
    aux_correl_coef = nan;   
end

if ~isfield(handles,'cvelocity')
    handles.cvelocity = aux_cvelocity;
    handles.correl_coef = aux_correl_coef;
else
    handles.cvelocity = [handles.cvelocity; aux_cvelocity];
    handles.correl_coef = [handles.correl_coef; aux_correl_coef];
end

handles.cvelocity = mean(handles.cvelocity);
handles.correl_coef = mean(handles.correl_coef);

set(handles.text_cond_velocity,'String', num2str(aux_cvelocity, 4));
set(handles.text_count_number,'String', strcat('#',num2str(handles.count)));

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in button_default_values.
function button_default_values_Callback(hObject, eventdata, handles)
% hObject    handle to button_default_values (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.count = 0;

if isfield(handles,'hcv_channels')
    for i=1:length(handles.hcv_channels)
        if ishandle(handles.hcv_channels{i})
            delete(handles.hcv_channels{i})
        end
    end
end

if isfield(handles,'hmap_amp')
    if ishandle(handles.hmap_amp)
        delete(handles.hmap_amp)
    end
end

if isfield(handles,'hmarkers')
    if ishandle(handles.hmarkers)
        delete(handles.hmarkers)
        delete(handles.hcog)
    end
end

if isfield(handles,'cvelocity')
    handles = rmfield(handles, 'cvelocity');
    handles = rmfield(handles, 'correl_coef');
end

if isfield(handles.data,'MEP_cluster')
    handles.data = rmfield(handles.data, 'MEP_cluster');
    handles.data = rmfield(handles.data, 'COG_cluster');
end

set(handles.checkbox_map_rms,'Value', 1);
checkbox_map_rms_Callback(hObject, eventdata, handles)

set(handles.text_cond_velocity,'String', '0');
set(handles.text_count_number,'String', strcat('#',num2str(handles.count)));
set(handles.edit_amp_threshold,'String', '70');

handles.hmap_amp = image(handles.amp_map, 'CDataMapping','scaled',...
    'parent', handles.axes2);
colorbar('peer', handles.axes2);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in button_mep_image.
function button_mep_image_Callback(hObject, eventdata, handles)
% hObject    handle to button_mep_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'hmap_amp')
    if ishandle(handles.hmap_amp)
        delete(handles.hmap_amp)
    end
end

color_min = str2double(get(handles.edit_image_color_min, 'String'));
color_max = str2double(get(handles.edit_image_color_max, 'String'));

amp_map_interp = interp2(handles.amp_map, 8);

handles.hmap_amp = image(amp_map_interp, 'CDataMapping','scaled', 'parent', handles.axes2);
colorbar('peer', handles.axes2)

set(handles.axes2,'xtick',1:5,'ytick',1:13, 'ydir','reverse',...
    'CLimMode', 'manual','CLim', [color_min color_max])

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in button_ok.
function button_ok_Callback(hObject, eventdata, handles)
% hObject    handle to button_ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'cv_channels');
    handles.data.cv_channels = handles.cv_channels;
end
if isfield(handles,'cvelocity');
    handles.data.cvelocity = handles.cvelocity;
end
if isfield(handles,'correl_coef');
    handles.data.correl_coef = handles.correl_coef;
end

handles.data.iz_row = str2double(get(handles.edit_iz_row,'String'));
handles.output = handles.data;

uiresume

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in button_save_figure.
function button_save_figure_Callback(hObject, eventdata, handles)
% hObject    handle to button_save_figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure_name_default = [num2str(handles.data.configuration.patient_id) '_' handles.data.configuration.emg_side '_' num2str(handles.data.configuration.angle_stim)];

[figure_name, figure_path] = uiputfile({'*.emf','Enhanced metafile (*.emf)';...
    '*.bmp','Windows bitmap (*.bmp)'},...
    'Save figure as...', figure_name_default);

saveas(handles.figure1, [figure_path figure_name])

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in checkbox_map_pp.
function checkbox_map_pp_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_map_pp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_map_pp

handles.map_pp_value = get(handles.checkbox_map_pp,'Value');

if handles.map_pp_value
    set(handles.checkbox_map_rms,'Value', 0);
    handles.amp_map = handles.data.amp_pp_map{handles.map_channel}(:,:);
else
    set(handles.checkbox_map_rms,'Value', 1);
    handles.amp_map = handles.data.amp_rms_map{handles.map_channel}(:,:);
end

handles.map_min = min(min(handles.amp_map));
handles.map_max = max(max(handles.amp_map));
set(handles.edit_image_color_min, 'String', num2str(handles.map_min,3))
set(handles.edit_image_color_max, 'String', num2str(handles.map_max,3))

button_mep_image_Callback(hObject, eventdata, handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in checkbox_map_rms.
function checkbox_map_rms_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_map_rms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_map_rms

handles.map_rms_value = get(handles.checkbox_map_rms,'Value');

if handles.map_rms_value
    set(handles.checkbox_map_pp,'Value', 0);
    handles.amp_map = handles.data.amp_rms_map{handles.map_channel}(:,:);
else
    set(handles.checkbox_map_pp,'Value', 1);
    handles.amp_map = handles.data.amp_pp_map{handles.map_channel}(:,:);
end

handles.map_min = min(min(handles.amp_map));
handles.map_max = max(max(handles.amp_map));
set(handles.edit_image_color_min, 'String', num2str(handles.map_min,3))
set(handles.edit_image_color_max, 'String', num2str(handles.map_max,3))

button_mep_image_Callback(hObject, eventdata, handles)

% Update handles structure
guidata(hObject, handles);


function edit_amp_threshold_Callback(hObject, eventdata, handles)
% hObject    handle to edit_amp_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_amp_threshold as text
%        str2double(get(hObject,'String')) returns contents of edit_amp_threshold as a double

function edit_ied_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ied (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ied as text
%        str2double(get(hObject,'String')) returns contents of edit_ied as a double

function edit_image_color_max_Callback(hObject, eventdata, handles)
% hObject    handle to edit_image_color_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_image_color_max as text
%        str2double(get(hObject,'String')) returns contents of edit_image_color_max as a double

function edit_image_color_min_Callback(hObject, eventdata, handles)
% hObject    handle to edit_image_color_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_image_color_min as text
%        str2double(get(hObject,'String')) returns contents of edit_image_color_min as a double


% --- Executes during object creation, after setting all properties.
function edit_amp_threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_amp_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

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

% --- Executes during object creation, after setting all properties.
function edit_ied_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ied (see GCBO)
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
