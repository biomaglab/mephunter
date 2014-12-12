function varargout = MAPClusters(varargin)
% MAPCLUSTERS MATLAB code for MAPClusters.fig
%      MAPCLUSTERS, by itself, creates a new MAPClUSTERS or raises the existing
%      singleton*.
%
%      H = MAPClUSTERS returns the handle to a new MAPClUSTERS or the handle to
%      the existing singleton*.
%
%      MAPClUSTERS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAPClUSTERS.M with the given input arguments.
%
%      MAPClUSTERS('Property','Value',...) creates a new MAPClUSTERS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MAPClusters_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MAPClusters_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MAPClusters

% Last Modified by GUIDE v2.5 14-Oct-2013 17:45:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MAPClusters_OpeningFcn, ...
                   'gui_OutputFcn',  @MAPClusters_OutputFcn, ...
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


% --- Executes just before MAPClusters is made visible.
function MAPClusters_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MAPClusters (see VARARGIN)

if ~isempty (varargin)
    handles.data = varargin{1};
    handles.signal = varargin{2};
end

% center the figure window on the screen
movegui(hObject, 'center');

if isfield(handles.data, 'iz_row')
    set(handles.edit_iz_row,'String', num2str(handles.data.iz_row));
else
    handles.data.iz_row = 13;
    set(handles.edit_iz_row,'String', num2str(handles.data.iz_row));
end
    

handles.amp_map = handles.data.amp_pp_map;
handles.data.cluster_indices_rms = cell(handles.signal.n_conditions, 1);
handles.data.cluster_indices_pp = cell(handles.signal.n_conditions, 1);
handles.data.cluster_union_rms = cell(handles.signal.n_conditions, 1);
handles.data.cluster_union_pp = cell(handles.signal.n_conditions, 1);
handles.data.abs_cluster_size_rms = cell(handles.signal.n_conditions, 1);
handles.data.abs_cluster_size_pp = cell(handles.signal.n_conditions, 1);

handles = cluster_coordinates(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MAPClusters wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MAPClusters_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiwait(hObject);
handles = guidata(hObject);

delete(hObject);

% Get default command line output from handles structure
varargout{1} = handles.data;


% --- Executes on button press in button_cluster_threshold.
function button_cluster_threshold_Callback(hObject, eventdata, handles)
% hObject    handle to button_cluster_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'hmap_amp')
    for i = 1:handles.signal.n_conditions
        if ishandle(handles.hmap_amp{i})
            delete(handles.hmap_amp{i})
            delete(handles.hcog{i})
        end
    end
end

if isfield(handles,'hmarkers_union')
    if ishandle(handles.hmarkers_union)
        delete(handles.hmarkers_union)
        delete(handles.hcog_union)
    end
end

handles = cluster_coordinates(handles);


% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in button_default_values.
function button_default_values_Callback(hObject, eventdata, handles)
% hObject    handle to button_default_values (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'hmarkers')
    for i = 1:handles.signal.n_conditions
        if ishandle(handles.hmarkers{i})
            delete(handles.hmarkers{i})
        end
    end
end

if isfield(handles,'hmarkers_union')
    if ishandle(handles.hmarkers_union)
        delete(handles.hmarkers_union)
        delete(handles.hcog_union)
    end
end

set(handles.edit_cluster_threshold,'String', '70');

set(handles.radiobutton_map_rms,'Value', 0);
set(handles.radiobutton_map_pp,'Value', 1);

handles.amp_map = handles.data.amp_pp_map;
button_cluster_threshold_Callback(hObject, eventdata, handles)

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

figure_name_default = [num2str(handles.signal.patient_id) '_' handles.signal.emg_side '_cluster'];

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

if isfield(handles,'hmarkers')
    for i = 1:handles.signal.n_conditions
        if ishandle(handles.hmarkers{i})
            delete(handles.hmarkers{i})
            delete(handles.hcog{i})
        end
    end
end

if isfield(handles,'hmarkers_union')
    if ishandle(handles.hmarkers_union)
        delete(handles.hmarkers_union)
        delete(handles.hcog_union)
    end
end

handles.map_pp_value = get(handles.radiobutton_map_pp,'Value');

if handles.map_pp_value
    set(handles.radiobutton_map_rms,'Value', 0);
    handles.amp_map = handles.data.amp_pp_map;
else
    set(handles.radiobutton_map_rms,'Value', 1);
    handles.amp_map = handles.data.amp_rms_map;
end

handles = cluster_coordinates(handles);

% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in radiobutton_map_rms.
function radiobutton_map_rms_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_map_rms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_map_rms

if isfield(handles,'hmarkers')
    for i = 1:handles.signal.n_conditions
        if ishandle(handles.hmarkers{i})
            delete(handles.hmarkers{i})
            delete(handles.hcog{i})
        end
    end
end

if isfield(handles,'hmarkers_union')
    if ishandle(handles.hmarkers_union)
        delete(handles.hmarkers_union)
        delete(handles.hcog_union)
    end
end

handles.map_rms_value = get(handles.radiobutton_map_rms,'Value');

if handles.map_rms_value
    set(handles.radiobutton_map_pp,'Value', 0);
    handles.amp_map = handles.data.amp_rms_map;
else
    set(handles.radiobutton_map_pp,'Value', 1);
    handles.amp_map = handles.data.amp_pp_map;
end

handles = cluster_coordinates(handles);

% Update handles structure
guidata(hObject, handles);


function edit_cluster_threshold_Callback(hObject, eventdata, handles)
% hObject    handle to edit_cluster_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_cluster_threshold as text
%        str2double(get(hObject,'String')) returns contents of edit_cluster_threshold as a double


% --- Executes during object creation, after setting all properties.
function edit_cluster_threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_cluster_threshold (see GCBO)
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

if isfield(handles,'hmarkers')
    for i = 1:handles.signal.n_conditions
        if ishandle(handles.hmarkers{i})
            delete(handles.hmarkers{i})
            delete(handles.hcog{i})
        end
    end
end

if isfield(handles,'hmarkers_union')
    if ishandle(handles.hmarkers_union)
        delete(handles.hmarkers_union)
        delete(handles.hcog_union)
    end
end

handles = cluster_coordinates(handles);

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

function out = cluster_coordinates(handles)

cog = cell(handles.signal.n_conditions,1);
MEP_rows = cell(handles.signal.n_conditions,1);
MEP_columns = cell(handles.signal.n_conditions,1);
cluster_indices = cell(handles.signal.n_conditions, 1);
cluster_union = [];
total_map = zeros(13,5);

% coordinates of electrodes
[x, y] = meshgrid(1:5, 1:13);
% reshaped coordinates map
xr = reshape(x,65,1);
yr = reshape(y,65,1);

for i = 1:handles.signal.n_conditions
    % Identifying clusters of activity from MEPs' maps
    % Peak-to-peak values independ the signal - Victor Hugo Souza (24/10/2013)
    aux_map = abs(handles.amp_map{i});
    % Limit the domain of cluster calculation to the region above the
    % innervation zone - Victor Hugo Souza (24/10/2013)
    % Don not use this for monopolar MEP
%     map_cut = aux_map(handles.data.iz_row:end,:);
%     aux_map = map(handles.data.iz_row:end,:);
    aux_map(1:handles.data.iz_row,:) = 0;
    % ---
    aux_map(isnan(aux_map))=0;
    cmean = [];
    cnumberofchannels = [];
    
    % the WSCalculation does not accept the aux_map matrix with all zeros
    % this try/catch statement treats this error
    try
        Local_MEP = WSCalculation(aux_map,1);
        
        % Exclude channels from the border that are being used in this map
        % Victor Hugo Souza (24/10/2013)
        Local_MEP(1,1) = 0;
        Local_MEP(1,end) = 0;
%         Local_MEP(end,1) = 0;
%         Local_MEP(end,end) = 0;
        % ---
        
        for clusters = 1:max(Local_MEP(:))
            % mean MEP activity per cluster
            cmean(clusters) = mean(aux_map(Local_MEP==clusters));
            % number of channels per cluster
            cnumberofchannels(clusters) = numel(Local_MEP(Local_MEP==clusters));
            if isnan(cmean(clusters))
                cmean(clusters)=0;
                cnumberofchannels(clusters) = inf;
            end
        end
        cmean(cnumberofchannels<5) = [];
        cnumberofchannels(cnumberofchannels<5) = [];
        % cpos == cluster of channels detecting greatest MEPs
        [~,cpos] = max(cmean);
        % indices pointing to channels detecting greatest MEPs
        cind = find(Local_MEP == cpos);
        
        threshold = str2double(get(handles.edit_cluster_threshold,'String'))/100*max(aux_map(cind));
        
%         % By Victor Hugo to fix reduction of aux_map size
%         aux_map = cat(1,zeros(handles.data.iz_row-1, size(map,2)),aux_map);
        
        % identifying channels detecting MEPs over threshold
        MEP_channels = aux_map(cind) >= threshold;
        MEP_channels = cind(MEP_channels);
        [MEP_rows{i},MEP_columns{i}] = ind2sub(size(aux_map),MEP_channels);
        %     handles.data.mepcluster{handles.data.channel} = [MEP_rows, MEP_columns];
        
        % creating binary mask for center of gravity
        mask1 = zeros(size(aux_map));
        for w = 1:length(MEP_rows{i})
            mask1(MEP_rows{i}(w), MEP_columns{i}(w)) = 1;
        end
        mask1(1, [1 end]) = 0;
        mask1(end, [1 end]) = 0;
        
        % this block is to apply clusterdata with spearman
%         % amplitudes for electrodes on cluster from watershed
%         d1 = mask1.*aux_map;
%         % reshaped amplitude and coordinates map for 0 degrees
%         d1r = reshape(d1,size(aux_map,1)*size(aux_map,2),1);
%         % data set prepared for hierarchical clustering (X, Y, AMP)
%         dataset = horzcat(xr, yr, d1r);
%         % cls1r = clusterdata(dataset, 'maxclust', 2);
%         cls1r = clusterdata(dataset, 'criterion', 'distance',...
%             'distance', 'spearman', 'maxclust', 3);
%         % find max peak-peak amplitudes
%         ind = d1r == max(d1r);
%         clmax = (cls1r == cls1r(ind)).*cls1r;
%         mask2 = 1.0.*((reshape(clmax,size(aux_map,1),size(aux_map,2)))~=0);
%         [rows, cols] = find(reshape(clmax,13,5));
%         MEP_rows{i} = rows;
%         MEP_columns{i} = cols;
        % end of spearman block
        
        mask = clusterfilt(mask1);
        [MEP_rows{i}, MEP_columns{i}] = find(mask);
        % excluding channels under the innervation zone
        % Do not use for monopolar MEP
%         rmv_rows = find(MEP_rows{i}>=handles.data.iz_row);
%         MEP_rows{i}(rmv_rows) = [];
%         MEP_columns{i}(rmv_rows) = [];
%         mask(handles.data.iz_row:end,:) = 0;
        
        % finding cluster indices and creating the union of all
        % orientations - possible area for MEPs
        cluster_indices{i} = find(mask);
        cluster_union = union(cluster_union, cluster_indices{i});
        
        % calculating the cog by imaging methods
        map_gray = mat2gray(aux_map);
        map_gray = uint8(round(map_gray*255));
        map_prop = regionprops(mask, map_gray, 'WeightedCentroid');
        cog{i} = round(cat(1, map_prop.WeightedCentroid));
        
        if get(handles.radiobutton_map_rms,'Value')
            handles.data.rms_entropy(i) = entropy(map_gray);
        else
            handles.data.pp_entropy(i) = entropy(map_gray);
        end
        total_map = total_map + handles.amp_map{i};
    catch
        cluster_indices{i} = [];
    end
   
    clear cpos threshold aux_map cmean cnumberofchannels

end

[cl_rows,cl_cols] = ind2sub(size(mask),cluster_union);
% cl_cols = ceil(cluster_union/13);
% cl_rows = mod(cluster_union, 13);
% cl_rows(cl_rows==0)=13;
mask_union = zeros(size(mask));
for w = 1:length(cl_cols)
    mask_union(cl_rows(w), cl_cols(w)) = 1;
end
total_map = total_map/handles.signal.n_conditions;
map_gray2 = mat2gray(total_map);
map_gray2 = uint8(round(map_gray2*255));
map_prop2 = regionprops(mask_union, map_gray2, 'WeightedCentroid');
cog_union = round(cat(1, map_prop2.WeightedCentroid));

for i = 1:handles.signal.n_conditions
    axes(eval(strcat('handles.axes',num2str(i))));
    set(gca,'xlim',[0 6],'ylim',[0 14], 'xtick', 1:5,'xticklabel',[],...
        'ytick',1:13,'yticklabel', [],'ydir','reverse','box','on','xgrid','on','ygrid','on')
    if i == 1 || i == 5
        if i == 5
            set(gca,'xticklabel',1:5,'yticklabel', 1:13)
        else
            set(gca, 'yticklabel', 1:13)
        end
    end
    if i > 5
        set(gca,'xticklabel',1:5)
    end
    set(eval(strcat('handles.text_angle',num2str(i))), 'String',...
        [num2str(handles.signal.angle_stim(i)) 'º'])
    hold on
    
    handles.hmarkers{i} = plot(MEP_columns{i},MEP_rows{i},'o', 'markersize', 10,'markerfacecolor','y');
    
    if sum(size(cog{i})) ~= 0
        handles.data.cog_col(i) = cog{i}(:, 1);
        handles.data.cog_row(i) = cog{i}(:, 2);
        handles.hcog{i} = plot(cog{i}(:, 1),cog{i}(:, 2),'o', 'markersize', 10, 'markerfacecolor', 'g');
    else
        handles.data.cog_col(i) = 0;
        handles.data.cog_row(i) = 0;
        handles.hcog{i} = [];
    end
    hold off
end

set(handles.axes9,'xlim', [0 6],'ylim',[0 14],...
    'xtick',1:5,'ytick',1:13,'ydir','reverse','box','on',...
    'xgrid','on','ygrid','on')

% saving the cluster results for RMS and P-P maps
if get(handles.radiobutton_map_rms,'Value')
    for i = 1:handles.signal.n_conditions
        handles.data.cluster_size_rms{i} = length(cluster_indices{i})/length(cluster_union);
        handles.data.abs_cluster_size_rms{i} = length(cluster_indices{i});
    end
    handles.data.cluster_indices_rms = cluster_indices;
    handles.data.cluster_union_rms = cluster_union;
else
    for i = 1:handles.signal.n_conditions
        handles.data.cluster_size_pp{i} = length(cluster_indices{i})/length(cluster_union);
        handles.data.abs_cluster_size_pp{i} = length(cluster_indices{i});
    end
    handles.data.cluster_indices_pp = cluster_indices;
    handles.data.cluster_union_pp = cluster_union;
end

% cl_cols = ceil(cluster_union/13);
% cl_rows = mod(cluster_union, 13);

axes(handles.axes9);
hold on
handles.hmarkers_union = plot(cl_cols, cl_rows,'o', 'markersize', 10,'markerfacecolor','y');
if ~isempty(cog_union)
    handles.hcog_union = plot(cog_union(:, 1),cog_union(:, 2),'o', 'markersize', 10, 'markerfacecolor', 'g');
end
hold off

out = handles;
