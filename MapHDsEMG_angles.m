function varargout = MapHDsEMG_angles(varargin)
% MAPHDSEMG_ANGLES MATLAB code for MapHDsEMG_angles.fig
%      MAPHDSEMG_ANGLES, by itself, creates a new MAPHDSEMG_ANGLES or raises the existing
%      singleton*.
%
%      H = MAPHDSEMG_ANGLES returns the handle to a new MAPHDSEMG_ANGLES or the handle to
%      the existing singleton*.
%
%      MAPHDSEMG_ANGLES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAPHDSEMG_ANGLES.M with the given input arguments.
%
%      MAPHDSEMG_ANGLES('Property','Value',...) creates a new MAPHDSEMG_ANGLES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MapHDsEMG_angles_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MapHDsEMG_angles_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MapHDsEMG_angles

% Last Modified by GUIDE v2.5 01-Oct-2013 15:52:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MapHDsEMG_angles_OpeningFcn, ...
                   'gui_OutputFcn',  @MapHDsEMG_angles_OutputFcn, ...
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


% --- Executes just before MapHDsEMG_angles is made visible.
function MapHDsEMG_angles_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MapHDsEMG_angles (see VARARGIN)

if ~isempty (varargin)
    handles.data = varargin{1};
end

% loading signal and configuration data
[signal_name, signal_path, ~] = uigetfile({'*.mat','MAT-files (*.mat)'},...
    'Select the signal file');

% [config_name, config_path, ~] = uigetfile({'*.mat','MAT-files (*.mat)'},...
%     'Select the configuration file');

tic
handles.signal_file = [signal_path signal_name];
% handles.config_file = [config_path config_name];
signal_vars = {'signal_path', 'n_conditions', 'angle_stim', 'fsample',...
    'n_channels', 'patient_id','emg_side', 'signal_gain', 'ad_bits', 'emg_diff'};
handles.signal = load(handles.signal_file, signal_vars{:});

% handles.signal = load(handles.signal_file, 'emg_diff');
% handles.configuration = load(handles.signal_file, config_vars{:});
toc

% flag for angle being visualized
handles.id_angle = 1;

% flags for checkbox visulizaton
handles.signalvalue = 0;
handles.mepmeanvalue = 0;
handles.mepsvalue = 0;
handles.linevalue = 0;
handles.triggervalue = 0;
handles.mepminmaxvalue = 0;
handles.ampthresholdvalue = 0;
handles.latencyvalue = 0;

% MEP windowing and trigger
handles.data.slope = str2double(get(handles.edit_Threshold,'String'));
handles.MEPStart = round(str2double(get(handles.edit_MEPStart,'String')))/1000;
handles.data.s0 = zeros(1, handles.signal.n_conditions);
handles.data.s1 = zeros(1, handles.signal.n_conditions);

% cell variables declaration
handles.data.outliers_mep = cell(1, handles.signal.n_conditions);
handles.data.ch_excluded = cell(1, handles.signal.n_conditions);
handles.data.latency = cell(1, handles.signal.n_conditions);
handles.data.duration = cell(1, handles.signal.n_conditions);
handles.data.amp_pp_map = cell(1, handles.signal.n_conditions);
handles.data.amp_rms_map = cell(1, handles.signal.n_conditions);
handles.data.offset_rms_map = cell(1, handles.signal.n_conditions);
handles.data.amp_fmed_map = cell(1, handles.signal.n_conditions);
handles.data.cluster_amp_pp = cell(1, handles.signal.n_conditions);
handles.data.cluster_amp_rms = cell(1, handles.signal.n_conditions);
handles.data.fmed = cell(1, handles.signal.n_conditions);
handles.data.amp_rms = cell(1, handles.signal.n_conditions);
handles.data.offset_fmed = cell(1, handles.signal.n_conditions);
handles.data.offset_rms = cell(1, handles.signal.n_conditions);
handles.data.trigger =  cell(1, handles.signal.n_conditions);
handles.data.cluster_size_rms = cell(handles.signal.n_conditions, 1);
handles.data.cluster_size_pp = cell(handles.signal.n_conditions, 1);
handles.data.abs_cluster_size_pp = cell(handles.signal.n_conditions, 1);
handles.data.abs_cluster_size_rms = cell(handles.signal.n_conditions, 1);
handles.data.iz_row = [];

% array variables declaration
handles.data.mean_total_rms = zeros(1, handles.signal.n_conditions);
handles.data.mean_total_pp = zeros(1, handles.signal.n_conditions);
handles.data.mean_offset_rms = zeros(1, handles.signal.n_conditions);
handles.data.rms_entropy = zeros(handles.signal.n_conditions, 1);
handles.data.pp_entropy = zeros(handles.signal.n_conditions, 1);

angles_names = '';

% create latency and duration cells, x sampled, initial (s0) and final (s1)
% samples
for i = 1:handles.signal.n_conditions
    angles_names = strcat(angles_names,num2str(i),'-',...
        num2str(handles.signal.angle_stim(i)),'_');
    handles.data.s0(i) = round(str2double(get(handles.edit_MEPStart,'String'))*handles.signal.fsample(i)/1000);
    handles.data.s1(i) = round(str2double(get(handles.edit_MEPEnd,'String'))*handles.signal.fsample(i)/1000);
    handles.data.latency{i} = cell(1, handles.signal.n_channels);
    handles.data.duration{i} = cell(1, handles.signal.n_channels);
    for j = 1:handles.signal.n_channels
        handles.hlatencystart{i, j} = [];
        handles.hlatencystop{i, j} = [];
        handles.data.xs{i} = (1:length(handles.signal.emg_diff{i}(:,j)))/handles.signal.fsample(i);
    end
end

% figure name
name = get(handles.figure1, 'Name');
fig_name = [name ' ' handles.signal.patient_id '_' handles.signal.emg_side];
set(handles.figure1,'Name',fig_name)

% text editors values
set(handles.text_ChannelLegend, 'String',angles_names)
set(handles.text_Channel, 'String',num2str(handles.signal.angle_stim(handles.id_angle)))
set(handles.edit_id_number, 'String',num2str(handles.signal.patient_id))
set(handles.edit_hemisphere_side, 'String',num2str(handles.signal.emg_side))
set(handles.edit_number_angle_stim, 'String',num2str(handles.signal.angle_stim(handles.id_angle)))

% handles = visibleoff(handles);

% trigger and meps calculations
handles = trigger_calculation(handles);
handles = meps_calculation(handles);

% updating checkbox values for MEP visualization
set(handles.checkbox_ViewMEPsMean,'Value',1)
set(handles.checkbox_ViewMinMax,'Value',1)
set(handles.checkbox_ViewLatency,'Value',1)
checkbox_ViewMEPsMean_Callback(hObject, eventdata, handles)
checkbox_ViewMinMax_Callback(hObject, eventdata, handles)
checkbox_ViewLatency_Callback(hObject, eventdata, handles)

% Choose default command line output for MapHDsEMG_angles
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MapHDsEMG_angles wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MapHDsEMG_angles_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in button_Trigger.
function button_Trigger_Callback(hObject, eventdata, handles)
% hObject    handle to button_Trigger (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% set current channel plots as invisible
handles = visibleoff(handles);

handles.data.slope = str2double(get(handles.edit_Threshold,'String'));

handles = trigger_calculation(handles);

set(handles.checkbox_ViewSignal,'Value',1)
set(handles.checkbox_ViewTrigger,'Value',1)
set(handles.checkbox_ViewLine,'Value',1)

checkbox_ViewSignal_Callback(hObject, eventdata, handles)
checkbox_ViewTrigger_Callback(hObject, eventdata, handles)
checkbox_ViewLine_Callback(hObject, eventdata, handles)

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in button_Backward.
function button_Backward_Callback(hObject, eventdata, handles)
% hObject    handle to button_Backward (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% current values of each checkbox
a1 = get(handles.checkbox_ViewSignal,'Value');
a2 = get(handles.checkbox_ViewMEPsMean,'Value');
a3 = get(handles.checkbox_ViewMEPsMult,'Value');
a4 = get(handles.checkbox_ViewLine,'Value');
a5 = get(handles.checkbox_ViewTrigger,'Value');
a6 = get(handles.checkbox_AmplitudeThreshold,'Value');
a7 = get(handles.checkbox_ViewMinMax,'Value');
a8 = get(handles.checkbox_ViewLatency,'Value');

% set current channel plots as invisible
handles = visibleoff(handles);

% restore checkbox values
handles.signalvalue = a1;
handles.mepmeanvalue = a2;
handles.mepsvalue = a3;
handles.linevalue = a4;
handles.triggervalue = a5;
handles.ampthresholdvalue = a6;
handles.mepminmaxvalue = a7;
handles.latencyvalue = a8;

set(handles.checkbox_ViewSignal,'Value',handles.signalvalue);
set(handles.checkbox_ViewMEPsMean,'Value',handles.mepmeanvalue);
set(handles.checkbox_ViewMEPsMult,'Value',handles.mepsvalue);
set(handles.checkbox_ViewLine,'Value',handles.linevalue);
set(handles.checkbox_ViewTrigger,'Value',handles.triggervalue);
set(handles.checkbox_AmplitudeThreshold,'Value',handles.ampthresholdvalue);
set(handles.checkbox_ViewMinMax,'Value',handles.mepminmaxvalue);
set(handles.checkbox_ViewLatency,'Value',handles.latencyvalue);
% ------------------------

% change text channel
if handles.id_angle <= 1
    handles.id_angle = handles.signal.n_conditions;
    set(handles.edit_number_angle_stim, 'String',...
        num2str(handles.signal.angle_stim(handles.id_angle)))
    set(handles.text_Channel, 'String',...
        num2str(handles.signal.angle_stim(handles.id_angle)))
else
    handles.id_angle = handles.id_angle-1;
    set(handles.edit_number_angle_stim, 'String',...
        num2str(handles.signal.angle_stim(handles.id_angle)))
    set(handles.text_Channel, 'String',...
        num2str(handles.signal.angle_stim(handles.id_angle)))
end

% execute callback for each checkbox with current value
checkbox_ViewSignal_Callback(hObject, eventdata, handles)
checkbox_ViewMEPsMean_Callback(hObject, eventdata, handles)
checkbox_ViewMEPsMult_Callback(hObject, eventdata, handles)
checkbox_ViewLine_Callback(hObject, eventdata, handles)
checkbox_ViewTrigger_Callback(hObject, eventdata, handles)
checkbox_AmplitudeThreshold_Callback(hObject, eventdata, handles)
checkbox_ViewMinMax_Callback(hObject, eventdata, handles)
checkbox_ViewLatency_Callback(hObject, eventdata, handles)

% execute callback for edit_angles_stim
edit_number_angle_stim_Callback(hObject, eventdata, handles)

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in button_Forward.
function button_Forward_Callback(hObject, eventdata, handles)
% hObject    handle to button_Forward (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% current values of each checkbox
a1 = get(handles.checkbox_ViewSignal,'Value');
a2 = get(handles.checkbox_ViewMEPsMean,'Value');
a3 = get(handles.checkbox_ViewMEPsMult,'Value');
a4 = get(handles.checkbox_ViewLine,'Value');
a5 = get(handles.checkbox_ViewTrigger,'Value');
a6 = get(handles.checkbox_AmplitudeThreshold,'Value');
a7 = get(handles.checkbox_ViewMinMax,'Value');
a8 = get(handles.checkbox_ViewLatency,'Value');

% set current channel plots as invisible
handles = visibleoff(handles);

% restore checkbox values
handles.signalvalue = a1;
handles.mepmeanvalue = a2;
handles.mepsvalue = a3;
handles.linevalue = a4;
handles.triggervalue = a5;
handles.ampthresholdvalue = a6;
handles.mepminmaxvalue = a7;
handles.latencyvalue = a8;

set(handles.checkbox_ViewSignal,'Value',handles.signalvalue);
set(handles.checkbox_ViewMEPsMean,'Value',handles.mepmeanvalue);
set(handles.checkbox_ViewMEPsMult,'Value',handles.mepsvalue);
set(handles.checkbox_ViewLine,'Value',handles.linevalue);
set(handles.checkbox_ViewTrigger,'Value',handles.triggervalue);
set(handles.checkbox_AmplitudeThreshold,'Value',handles.ampthresholdvalue);
set(handles.checkbox_ViewMinMax,'Value',handles.mepminmaxvalue);
set(handles.checkbox_ViewLatency,'Value',handles.mepminmaxvalue);
% ------------------------

% change text channel
if handles.id_angle >= handles.signal.n_conditions    
   handles.id_angle = 1; 
   set(handles.edit_number_angle_stim, 'String',...
       num2str(handles.signal.angle_stim(handles.id_angle)))
   set(handles.text_Channel, 'String',...
       num2str(handles.signal.angle_stim(handles.id_angle)))
else
    handles.id_angle = handles.id_angle+1;        
    set(handles.edit_number_angle_stim, 'String',...
       num2str(handles.signal.angle_stim(handles.id_angle)))
   set(handles.text_Channel, 'String',...
       num2str(handles.signal.angle_stim(handles.id_angle)))
end

% execute callback for each checkbox with current value
checkbox_ViewSignal_Callback(hObject, eventdata, handles)
checkbox_ViewMEPsMean_Callback(hObject, eventdata, handles)
checkbox_ViewMEPsMult_Callback(hObject, eventdata, handles)
checkbox_ViewLine_Callback(hObject, eventdata, handles)
checkbox_ViewTrigger_Callback(hObject, eventdata, handles)
checkbox_AmplitudeThreshold_Callback(hObject, eventdata, handles)
checkbox_ViewMinMax_Callback(hObject, eventdata, handles)
checkbox_ViewLatency_Callback(hObject, eventdata, handles)

% execute callback for edit_angles_stim
edit_number_angle_stim_Callback(hObject, eventdata, handles)

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in button_MEPs.
function button_MEPs_Callback(hObject, eventdata, handles)
% hObject    handle to button_MEPs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% change checkbox values for MEP visualization
handles = visibleoff(handles);

handles = meps_calculation(handles);

% updating checkbox values for MEP visualization
set(handles.checkbox_ViewMEPsMean,'Value',1)
set(handles.checkbox_ViewMinMax,'Value',1)
set(handles.checkbox_ViewLatency,'Value',1)

checkbox_ViewMEPsMean_Callback(hObject, eventdata, handles)
checkbox_ViewMinMax_Callback(hObject, eventdata, handles)
checkbox_ViewLatency_Callback(hObject, eventdata, handles)

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in button_Export.
function button_Export_Callback(hObject, eventdata, handles)
% hObject    handle to button_Export (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tic
previous_data = [];

[filename, pathname, filterindex] = uiputfile({'*.xls;*.xlsx','MS Excel Files (*.xls,*.xlsx)';...
    '*.txt', 'ASCII format (*.txt)'}, 'Export data', 'processed_data.xlsx');

export_data = [handles.signal.signal_path' num2cell(repmat(handles.signal.patient_id,[8 1]))...
    num2cell(repmat(handles.signal.emg_side,[8 1])) num2cell(handles.signal.angle_stim')...
    num2cell(repmat(handles.data.iz_row,[8 1])) num2cell(handles.data.pp_entropy) num2cell(handles.data.rms_entropy)...
    num2cell(handles.data.mean_total_pp') num2cell(handles.data.mean_total_rms') num2cell(handles.data.mean_offset_rms')...
    handles.data.cluster_amp_pp' handles.data.abs_cluster_size_pp handles.data.cluster_amp_rms'...
    handles.data.abs_cluster_size_rms handles.data.cluster_size_pp handles.data.cluster_size_rms...
    num2cell(handles.data.dantec_1_pp) num2cell(handles.data.dantec_2_pp) num2cell(handles.data.seniam_2_pp)...
    num2cell(handles.data.dantec_3_pp) num2cell(handles.data.seniam_3_pp) num2cell(handles.data.dantec_1_fmed)...
    num2cell(handles.data.dantec_2_fmed) num2cell(handles.data.seniam_2_fmed) num2cell(handles.data.dantec_3_fmed)...
    num2cell(handles.data.seniam_3_fmed) num2cell(handles.data.dantec_1_rms) num2cell(handles.data.dantec_2_rms)...
    num2cell(handles.data.seniam_2_rms) num2cell(handles.data.dantec_3_rms) num2cell(handles.data.seniam_3_rms)];

headers = [{'file_name'} {'subject'} {'hemisphere'} {'angle(degrees)'}...
    {'iz(row)'} {'entropy_pp'} {'entropy_rms'} {'mean_total_pp(uV)'}...
    {'mean_total_rms(uV)'} {'mean_offset_rms(uV)'} {'cluster_amp_pp(uV)'}...
    {'abs_cluster_size_pp'} {'cluster_amp_rms(uV)'} {'abs_cluster_size_rms'}...
    {'cluster_size_pp'} {'cluster_size_rms'}...
    {'dantec_1_pp(uV)'} {'dantec_2_pp(uV)'} {'seniam_2_pp(uV)'}...
    {'dantec_3_pp(uV)'} {'seniam_3_pp(uV)'} {'dantec_1_freq(Hz)'}...
    {'dantec_2_freq(Hz)'} {'seniam_2_freq(Hz)'} {'dantec_3_freq(Hz)'}...
    {'seniam_3_freq(Hz)'} {'dantec_1_rms(uV)'} {'dantec_2_rms(uV)'}...
    {'seniam_2_rms(uV)'} {'dantec_3_rms(uV)'} {'seniam_2_rms(uV)'}];

switch filterindex
    case 1
        try
            [~, ~, previous_data] = xlsread([pathname filename]);
        end
        if isempty(previous_data)
            xlswrite([pathname filename], [headers; export_data])
        else
            xlswrite([pathname filename], [previous_data; export_data])
        end
        
    case 2
        fid = fopen([pathname filename]);
        try
            previous_data = fgets(fid);
        end
        the_format = '\n%s %d %s %d %.4f %.4f %.4f %.4f %.4f %.4f %s %d %.4f %.4f %.4f';
        if isempty(previous_data)
            fid = fopen([pathname filename], 'w');
            fprintf(fid, '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s', headers{1,:});
            fprintf(fid, the_format, export_data{1,:});
            fclose(fid);
        else
            fid = fopen([pathname filename], 'a');
            fprintf(fid, the_format, export_data{1,:})
            fclose(fid);
        end
end
toc

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in button_Save.
function button_Save_Callback(hObject, eventdata, handles)
% hObject    handle to button_Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tic
file_name_default = [num2str(handles.signal.patient_id) '_' handles.signal.emg_side '_' 'data'];
file_name_default = strcat(file_name_default,'.mat');

update_waitbar(handles,0.5)

[file_name, file_path, index] = uiputfile({'*.mat','MAT-files (*.mat)'},...
    'Save data as...', file_name_default);

if index == 1
    data = handles.data;
    save([file_path file_name], 'data')
end

update_waitbar(handles,1)
toc

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in button_Reset.
function button_Reset_Callback(hObject, eventdata, handles)
% hObject    handle to button_Reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% set current channel plots as invisible - increase reset velocity
handles = visibleoff(handles);

% delete all plots handles
for j = 1:handles.signal.n_channels
    for i = 1:handles.signal.n_conditions
        
        if isfield(handles,'hsignal')
            if ishandle(handles.hsignal{i, j})
                delete(handles.hsignal{i, j})
            end
        end
        
        if isfield(handles,'htrigger')
            if ishandle(handles.htrigger{i, j})
                delete(handles.htrigger{i, j})
            end
        end
        
        if isfield(handles,'hline')
            if ishandle(handles.hline{i, j})
                delete(handles.hline{i, j})
            end
        end
        
        if isfield(handles,'hampthreshold')
            if ishandle(handles.hampthreshold{i, j})
                delete(handles.hampthreshold{i, j})
            end
        end
        
        if isfield(handles,'hmeps')
            if ishandle(handles.hmeps{i, j})
                delete(handles.hmeps{i, j})
            end
        end
        if isfield(handles,'hmepmean')
            if ishandle(handles.hmepmean{i, j})
                delete(handles.hmepmean{i, j})
            end
        end
        
        if isfield(handles,'hmepmax')
            if ishandle(handles.hmepmax{i, j})
                delete(handles.hmepmax{i, j})
            end
        end
        
        if isfield(handles,'hmepmin')
            if ishandle(handles.hmepmin{i, j})
                delete(handles.hmepmin{i, j})
            end
            
        end
        if isfield(handles,'hlatencystart')
            if ishandle(handles.hlatencystart{i, j})
                delete(handles.hlatencystart{i, j})
            end
        end
        
        if isfield(handles,'hlatencystop')
            if ishandle(handles.hlatencystop{i, j})
                delete(handles.hlatencystop{i, j})
            end
        end
    end
    
    update_waitbar(handles,j/handles.signal.n_channels)
    axes(eval(strcat('handles.axes',num2str(j))))
    cla
end

% update checkbox values
set(handles.checkbox_ViewSignal,'Value',0);
set(handles.checkbox_ViewMEPsMean,'Value',0);
set(handles.checkbox_ViewMEPsMult,'Value',0);
set(handles.checkbox_ViewLine,'Value',0);
set(handles.checkbox_ViewTrigger,'Value',0);
set(handles.checkbox_AmplitudeThreshold,'Value',0);
set(handles.checkbox_ViewMinMax,'Value',0);
set(handles.checkbox_ViewLatency,'Value',0);

% clear plots handles variables
clear handles.hsignal handles.htrigger handles.hline handles.hampthreshold...
    handles.hmeps handles.hmepmean handles.hmepmax handles.hmepmin...
    handles.hlatencystart handles.hlatencystop

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in checkbox_AmplitudeThreshold.
function checkbox_AmplitudeThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_AmplitudeThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_AmplitudeThreshold

% update amplitude threshold lines visualization according to checkbox value
ampthresholdvalue = get(handles.checkbox_AmplitudeThreshold,'Value');
handles.ampthresholdvalue = ampthresholdvalue;

if isfield(handles,'hampthreshold')
    for i = 1:size(handles.hampthreshold, 2)
        if ishandle(handles.hampthreshold{handles.id_angle, i})
            if ampthresholdvalue == 1
                set(handles.hampthreshold{handles.id_angle, i},'Visible','on')
            else
                set(handles.hampthreshold{handles.id_angle, i},'Visible','off')
            end
        end
    end
end

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in checkbox_ViewLatency.
function checkbox_ViewLatency_Callback(hObject, ~, handles)
% hObject    handle to checkbox_ViewLatency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_ViewLatency

% update latency visualization according to checkbox value
latencyvalue = get(handles.checkbox_ViewLatency,'Value');
handles.latencyvalue = latencyvalue;

if isfield(handles,'hlatencystart')
    for j = 1:handles.signal.n_channels
        if ishandle(handles.hlatencystart{handles.id_angle, j})
            if latencyvalue
                set(handles.hlatencystart{handles.id_angle, j},'Visible','on')
                set(handles.hlatencystop{handles.id_angle, j},'Visible','on')
            else
                set(handles.hlatencystart{handles.id_angle, j},'Visible','off')
                set(handles.hlatencystop{handles.id_angle, j},'Visible','off')
            end
        end
    end
end

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in checkbox_ViewSignal.
function checkbox_ViewSignal_Callback(hObject, ~, handles)
% hObject    handle to checkbox_ViewSignal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_ViewSignal

% update signal visualization according to checkbox value
signalvalue = get(handles.checkbox_ViewSignal,'Value');
handles.signalvalue = signalvalue;

if isfield(handles,'hsignal')
    for j = 1:handles.signal.n_channels
        if ishandle(handles.hsignal{handles.id_angle, j})
            if signalvalue
                set(handles.hsignal{handles.id_angle, j},'Visible','on')
            else
                set(handles.hsignal{handles.id_angle, j},'Visible','off')
            end
        end
    end
end

% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in checkbox_ViewMEPsMean.
function checkbox_ViewMEPsMean_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_ViewMEPsMean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_ViewMEPsMean

% update mep mean signal visualization according to checkbox value
mepmeanvalue = get(handles.checkbox_ViewMEPsMean,'Value');
handles.mepmeanvalue = mepmeanvalue;

if isfield(handles,'hmepmean')
    for j = 1:handles.signal.n_channels
        if ishandle(handles.hmepmean{handles.id_angle, j})
            if mepmeanvalue
                set(handles.hmepmean{handles.id_angle, j},'Visible','on')
            else
                set(handles.hmepmean{handles.id_angle, j},'Visible','off')
            end
        end
    end
end

% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in checkbox_ViewMEPsMult.
function checkbox_ViewMEPsMult_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_ViewMEPsMult (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_ViewMEPsMult

% update multiple mep signal visualization according to checkbox value
mepsvalue = get(handles.checkbox_ViewMEPsMult,'Value');
handles.mepsvalue = mepsvalue;

if isfield(handles,'hmeps')
    for j = 1:handles.signal.n_channels
        if ishandle(handles.hmeps{handles.id_angle, j})
            if mepsvalue
                set(handles.hmeps{handles.id_angle, j},'Visible','on')
            else
                set(handles.hmeps{handles.id_angle, j},'Visible','off')
            end
        end
    end
end

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in checkbox_ViewLine.
function checkbox_ViewLine_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_ViewLine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_ViewLine

% update trigger green line visualization according to checkbox value
linevalue = get(handles.checkbox_ViewLine,'Value');
handles.linevalue = linevalue;

if isfield(handles,'hline')
    for j = 1:handles.signal.n_channels
        if ishandle(handles.hline{handles.id_angle, j})
            if linevalue == 1
                set(handles.hline{handles.id_angle, j},'Visible','on')
            else
                set(handles.hline{handles.id_angle, j},'Visible','off')
            end
        end
    end
end

% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in checkbox_ViewTrigger.
function checkbox_ViewTrigger_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_ViewTrigger (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_ViewTrigger

% update trigger red circulus visualization according to checkbox value
triggervalue = get(handles.checkbox_ViewTrigger,'Value');
handles.triggervalue = triggervalue;

if isfield(handles,'htrigger')
    for j = 1:handles.signal.n_channels
        if ishandle(handles.htrigger{handles.id_angle, j})
            if triggervalue
                set(handles.htrigger{handles.id_angle, j},'Visible','on')
            else
                set(handles.htrigger{handles.id_angle, j},'Visible','off')
            end
        end
    end
end

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in checkbox_ViewMinMax.
function checkbox_ViewMinMax_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_ViewMinMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_ViewMinMax

%  update min and max red cross visualization according to checkbox value
mepminmaxvalue = get(handles.checkbox_ViewMinMax,'Value');
handles.mepminmaxvalue = mepminmaxvalue;

if isfield(handles,'hmepmax')
    for j = 1:handles.signal.n_channels
        if ishandle(handles.hmepmax{handles.id_angle, j})
            if mepminmaxvalue
                set(handles.hmepmax{handles.id_angle, j},'Visible','on')
                set(handles.hmepmin{handles.id_angle, j},'Visible','on')
            else
                set(handles.hmepmax{handles.id_angle, j},'Visible','off')
                set(handles.hmepmin{handles.id_angle, j},'Visible','off')
            end
        end
    end
end

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in button_AmplitudeThreshold.
function button_AmplitudeThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to button_AmplitudeThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'hampthreshold')
    for j= 1:handles.signal.n_channels
        if ishandle(handles.hampthreshold{handles.id_angle, j})
            delete(handles.hampthreshold{handles.id_angle, j})
        end
    end
end

% set current channel plots as invisible
handles = visibleoff(handles);

amp_mono = get(handles.edit_AmplitudeThresholdMono,'String');
handles.data.amp_mono = str2double(amp_mono)/2.0;
amp_diff = get(handles.edit_AmplitudeThresholdDiff,'String');
handles.data.amp_diff = str2double(amp_diff)/2.0;

for j = 1:handles.signal.n_channels
    axes(eval(strcat('handles.axes',num2str(j))))
    for i = 1:handles.signal.n_conditions
        
        if ~isempty(handles.data.mepmax{i}{j})
            shift_amp = (handles.data.mepmax{i}{j} + handles.data.mepmin{i}{j})/2;
            handles.data.amp_pp{i}{j} = handles.data.mepmax{i}{j} - handles.data.mepmin{i}{j};
            if handles.data.amp_pp{i}{j} == 0 || isnan(handles.data.amp_pp{i}{j})
                handles.data.amp_pp{i}{j} = NaN;
            end
        else
            handles.data.amp_pp{i}{j} = NaN;
            shift_amp = 0;
        end
        if abs(handles.data.amp_pp{i}{j}) < 2*handles.data.amp_mono
            if isfield(handles,'hmepmax')
                if ishandle(handles.hmepmax{i, j})
                    delete(handles.hmepmax{i, j})
                    handles.data.mepmax{i}{j} = [];
                end
            end
            if isfield(handles,'hmepmin')
                if ishandle(handles.hmepmin{i, j})
                    delete(handles.hmepmin{i, j})
                    handles.data.mepmin{i}{j} = [];
                end
            end
        end
        if sum(handles.data.trigger{i}(:,j)) ~= 0
            a = line([handles.MEPStart,handles.data.xs{i}(handles.data.s1(i)-handles.data.s0(i))+handles.MEPStart],...
                [handles.data.amp_mono+shift_amp,handles.data.amp_mono+shift_amp],'Color','g','Visible','off');
            b = line([handles.MEPStart,handles.data.xs{i}(handles.data.s1(i)-handles.data.s0(i))+handles.MEPStart],...
                [-handles.data.amp_mono+shift_amp,-handles.data.amp_mono+shift_amp],'Color','g','Visible','off');
        else
            a = NaN;
            b = NaN;
        end
        handles.hampthreshold{i, j} = [a b];
    end
    update_waitbar(handles,j/handles.signal.n_channels)
end

% updating checkbox values for MEP visualization
set(handles.checkbox_AmplitudeThreshold,'Value',1)
set(handles.checkbox_ViewMinMax,'Value',1)
set(handles.checkbox_ViewMEPsMean,'Value',1)
set(handles.checkbox_ViewLatency,'Value',1)

checkbox_ViewMinMax_Callback(hObject, eventdata, handles)
checkbox_AmplitudeThreshold_Callback(hObject, eventdata, handles)
checkbox_ViewMEPsMean_Callback(hObject, eventdata, handles)
checkbox_ViewLatency_Callback(hObject, eventdata, handles)

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in button_Load.
function button_Load_Callback(hObject, eventdata, handles)
% hObject    handle to button_Load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tic
button_Reset_Callback(hObject, eventdata, handles)

[data_name, data_path, index] = uigetfile({'*.mat','MAT-files (*.mat)'},...
    'Select the data processed file');

if index ~= 0
    data_file = strcat(data_path,data_name);
else
    return
end
        
loadfile = load (data_file);
handles.data = loadfile.data;
handles.data.amp_mono = 10;
handles.data.amp_diff = 10;

% Trick to round decimal numbers
handles.MEPStart = round((handles.data.s0(handles.id_angle))/handles.signal.fsample(handles.id_angle)*1000)/1000;

angles_names = '';

for i = 1:handles.signal.n_conditions
    angles_names = strcat(angles_names,num2str(i),'-',...
        num2str(handles.signal.angle_stim(i)),'_');
end

set(handles.text_ChannelLegend, 'String',angles_names)
set(handles.text_Channel,'String',num2str(handles.id_angle))
set(handles.edit_id_number, 'String',num2str(handles.signal.patient_id))
set(handles.edit_hemisphere_side, 'String',num2str(handles.signal.emg_side))
set(handles.edit_number_angle_stim, 'String',num2str(handles.signal.angle_stim(handles.id_angle)))
set(handles.edit_Threshold,'String',num2str(handles.data.slope))
set(handles.edit_AmplitudeThresholdDiff,'String',num2str(2*handles.data.amp_diff))
set(handles.edit_AmplitudeThresholdMono,'String',num2str(2*handles.data.amp_mono))
set(handles.edit_MEPStart,'String', num2str(round((handles.data.s0(handles.id_angle)*1000)/handles.signal.fsample(handles.id_angle))))
set(handles.edit_MEPEnd,'String', num2str(round((handles.data.s1(handles.id_angle)*1000)/handles.signal.fsample(handles.id_angle))))

aux_colors = [0.0117 0.0343 0.894; 0.0117 0.6078 0.894;...
    0.0117 0.89411 0.69411; 0.0117 0.8902 0.3255;...
    0.0196 0.8901 0.0196; 0.6274 0.8902 0.3255;...
    0.9804 0.9568 0.1137; 0.9804 0.6941 0.1137;...
    0.9804 0.3294 0.1137; 0.8627 0.0196 0.0196;...
    0.0 0.0 0.0; 0.0 0.0 0.0;...
    0.0 0.0 0.0; 0.0 0.0 0.0;...
    0.0 0.0 0.0; 0.0 0.0 0.0];

% plots
for j = 1:handles.signal.n_channels
    for  i = 1:handles.signal.n_conditions
        
        axes(eval(strcat('handles.axes',num2str(j))))
        hold on
        
        handles.hsignal{i, j} = plot(handles.data.xs{i}, handles.signal.emg_diff{i}(:,j),...
            'Visible','off');
        
        colors = aux_colors;
        if sum(handles.data.trigger{i}(:,j)) ~= 0
            aux_trigger = handles.data.trigger{i}(:,j);
            if ~isempty(handles.data.mepmax{i}{j})
                shift_amp = (handles.data.mepmax{i}{j} + handles.data.mepmin{i}{j})/2;
                handles.hmepmax{i, j} = plot(handles.data.xs{i}(handles.data.pos_max{i}{1, j})+handles.MEPStart,...
                    handles.data.mepmax{i}{j},'+r','Visible','off');
                handles.hmepmin{i, j} = plot(handles.data.xs{i}(handles.data.pos_min{i}{1, j})+handles.MEPStart,...
                    handles.data.mepmin{i}{j},'+r','Visible','off');
            else
                shift_amp = 0;
                handles.hmepmax{i, j} = nan;
                handles.hmepmin{i, j} = nan;
            end
            
            for k = 1:length(aux_trigger) % Number of stimulus
                handles.hmeps{i, j}(k) = plot(handles.data.xs{i}(1:handles.data.s1(i)-handles.data.s0(i))+handles.MEPStart,...
                    handles.data.meps_diff{i}{j}{k},'-','color',colors(k,:),'Visible', 'off');
            end
            handles.hline{i, j} = line([handles.data.xs{i}(handles.data.trigger{i}(:,j))',...
                handles.data.xs{i}(handles.data.trigger{i}(:,j))'],...
                [min(handles.signal.emg_diff{i}(:,j)),...
                max(handles.signal.emg_diff{i}(:,j))],...
                'Color','g','Visible','off');
            handles.hmepmean{i, j} = plot(handles.data.xs{i}(1:handles.data.s1(i)-handles.data.s0(i))+handles.MEPStart,...
                handles.data.mepmean{i}{j},'Visible','off');
            
            if i == 1
                b = line([handles.MEPStart,handles.data.xs{i}(handles.data.s1(i)-handles.data.s0(i))+handles.MEPStart],...
                    [handles.data.amp_mono+shift_amp,handles.data.amp_mono+shift_amp],'Color','g','Visible','off');
                c = line([handles.MEPStart,handles.data.xs{i}(handles.data.s1(i)-handles.data.s0(i))+handles.MEPStart],...
                    [-handles.data.amp_mono+shift_amp,-handles.data.amp_mono+shift_amp],'Color','g','Visible','off');
            else
                b = line([handles.MEPStart,handles.data.xs{i}(handles.data.s1(i)-handles.data.s0(i))+handles.MEPStart],...
                    [handles.data.amp_diff+shift_amp,handles.data.amp_diff+shift_amp],'Color','g','Visible','off');
                c = line([handles.MEPStart,handles.data.xs{i}(handles.data.s1(i)-handles.data.s0(i))+handles.MEPStart],...
                    [-handles.data.amp_diff+shift_amp,-handles.data.amp_diff+shift_amp],'Color','g','Visible','off');
            end
            
        else
            handles.hmeps{i, j} = nan;
            handles.hmepmean{i, j} = nan;
            handles.hline{i, j} = nan;
            b = nan;
            c = nan;
        end
        
        handles.hampthreshold{i, j} = [b c];
        
        if ~isempty(handles.data.latency{i}{j})
            handles.hlatencystart{i, j} = plot(handles.data.latency{i}{j}(1,1),...
                handles.data.latency{i}{j}(1,2),'gv','Visible','off');
            handles.hlatencystop{i, j} = plot(handles.data.latency{i}{j}(2,1),...
                handles.data.latency{i}{j}(2,2),'r^','Visible','off');
        else
            handles.hlatencystart{i, j} = nan;
            handles.hlatencystop{i, j} = nan;
        end
    end
    update_waitbar(handles,j/handles.signal.n_channels)
end

% updating checkbox values for MEP visualization
set(handles.checkbox_ViewMEPsMean,'Value',1)
set(handles.checkbox_ViewMinMax,'Value',1)
set(handles.checkbox_ViewLatency,'Value',1)

checkbox_ViewMEPsMean_Callback(hObject, eventdata, handles)
checkbox_ViewMinMax_Callback(hObject, eventdata, handles)
checkbox_ViewLatency_Callback(hObject, eventdata, handles)

toc

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in button_Map.
function button_Map_Callback(hObject, eventdata, handles)
% hObject    handle to button_Map (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

aux_amp_pp_map = zeros(13,5);
aux_amp_pp_map_2 = zeros(13,5);
aux_amp_rms_map = zeros(13,5);
aux_amp_rms_map_2 = zeros(13,5);
aux_offset_rms_map = zeros(13,5);
aux_offset_rms_map_2 = zeros(13,5);
aux_fmed_map = zeros(13,5);
aux_fmed_map_2 = zeros(13,5);


% refreshing peak to peak, rms and median frequency maps
for j = 1:handles.signal.n_channels
    for i = 1:handles.signal.n_conditions
        if ~isempty(handles.data.mepmax{i}{j})
            handles.data.amp_pp{i}{j} = handles.data.mepmax{i}{j} - handles.data.mepmin{i}{j};
            if handles.data.amp_pp{i}{j} == 0 || isnan(handles.data.amp_pp{i}{j})
                handles.data.amp_pp{i}{j} = NaN;
                handles.data.amp_rms{i}{j} = NaN;
                handles.data.offset_rms{i}{j} = NaN;
                handles.data.fmed{i}{j} = NaN;
            end
        else
            handles.data.amp_pp{i}{j} = NaN;
            handles.data.amp_rms{i}{j} = NaN;
            handles.data.offset_rms{i}{j} = NaN;
            handles.data.fmed{i}{j} = NaN;
        end
    end
end

% Peak-peak amplitude map
for i = 1:length(handles.data.amp_pp)
    for j = 1:length(handles.data.amp_pp{handles.id_angle})
        aux_amp_pp_map(j) = handles.data.amp_pp{i}{j};
    end
    aux_amp_pp_map_2 = cat(3,aux_amp_pp_map_2, aux_amp_pp_map);
end


% RMS and offset RMS amplitude maps
for i = 1:length(handles.data.amp_rms)
    for j = 1: length(handles.data.amp_rms{handles.id_angle})
        aux_amp_rms_map(j) = handles.data.amp_rms{i}{j};
        aux_offset_rms_map(j) = handles.data.offset_rms{i}{j};
    end
    aux_amp_rms_map_2 = cat(3,aux_amp_rms_map_2, aux_amp_rms_map);
    aux_offset_rms_map_2 = cat(3, aux_offset_rms_map_2, aux_offset_rms_map);
end

% Median frequency map
for i = 1:length(handles.data.fmed)
    for j = 1:length(handles.data.fmed{handles.id_angle})
        aux_fmed_map(j) = handles.data.fmed{i}{j};
    end
    aux_fmed_map_2 = cat(3,aux_fmed_map_2, aux_fmed_map);
end

aux_amp_pp_map_2(:,:,1)=[];
aux_amp_rms_map_2(:,:,1)=[];
aux_offset_rms_map_2(:,:,1)=[];
aux_fmed_map_2(:,:,1)=[];


for i = 1:handles.signal.n_conditions
    handles.data.amp_pp_map{i} = aux_amp_pp_map_2(:,:,i);
    handles.data.amp_rms_map{i} = aux_amp_rms_map_2(:,:,i);
    handles.data.offset_rms_map{i} = aux_offset_rms_map_2(:,:,i);
    handles.data.amp_fmed_map{i} = aux_fmed_map_2(:,:,i);
end

% Visualization of lines plots of raw MEPs
handles.data.iz_row = MAPLines_raw(handles.data, handles.signal);

% Visualization of lines plots of MEPs
handles.data.iz_row = MAPLines(handles.data, handles.signal);

% Visualization of images of differential amplitude MEPs
handles.data.iz_row = MAPAmplitudes(handles.data, handles.signal);

% Visualization of cluster plots of MEPs
handles.data = MAPClusters(handles.data, handles.signal);

% Calculus of avarege amplitude rms for all the used electrodes
for i = 1:handles.signal.n_conditions
    handles.data.mean_total_rms(i) = nanmean(nanmean(handles.data.amp_rms_map{i}));
    handles.data.mean_total_pp(i) = nanmean(nanmean(handles.data.amp_pp_map{i}));
    handles.data.mean_offset_rms(i) = nanmean(nanmean(handles.data.offset_rms_map{i}));
end

% Calculus of avarege amplitude rms for the cluster electrodes
total_amp_rms = zeros(1, handles.signal.n_conditions);
total_amp_pp = zeros(1, handles.signal.n_conditions);

for i = 1:handles.signal.n_conditions
    for j = 1:size(handles.data.cluster_indices_rms{i},1)
        if ~isempty(handles.data.cluster_indices_rms{i})
            if isnan(handles.data.amp_rms_map{i}(mod(handles.data.cluster_indices_rms{i}(j), 13), ceil(handles.data.cluster_indices_rms{i}(j)/13)))
                total_amp_rms(i) = total_amp_rms(i);
            else
                total_amp_rms(i) = total_amp_rms(i) + handles.data.amp_rms_map{i}(mod(handles.data.cluster_indices_rms{i}(j), 13),...
                    ceil(handles.data.cluster_indices_rms{i}(j)/13));
            end
        else
            total_amp_rms(i) = 0;
        end
    end
    
    for k = 1:size(handles.data.cluster_indices_pp{i},1)
        if ~isempty(handles.data.cluster_indices_pp{i})
            if isnan(handles.data.amp_pp_map{i}(mod(handles.data.cluster_indices_pp{i}(k), 13),ceil(handles.data.cluster_indices_pp{i}(k)/13)))
                total_amp_rms(i) = total_amp_rms(i);
            else
                total_amp_pp(i) = total_amp_pp(i) + handles.data.amp_pp_map{i}(mod(handles.data.cluster_indices_pp{i}(k), 13),...
                    ceil(handles.data.cluster_indices_pp{i}(k)/13));
            end
        else
            total_amp_pp(i) = 0;
        end
    end
    
    handles.data.cluster_amp_rms{i} = total_amp_rms(i);
    handles.data.cluster_amp_pp{i} = total_amp_pp(i);
end

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in button_Sim_MEP.
function button_Sim_MEP_Callback(hObject, eventdata, handles)
% hObject    handle to button_Sim_MEP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Simulated MEP for 3 different electrode positionings and 2 types of
% electrodes
%
% Type 1 - DANTEC 13k60 - size 6 x 12 x 1,5 mm with IED = 9,6 mm
% Type 2 - SENIAM for dynamic contractionc - 2 x 9,6 mm with IED = 4,8 mm
% IED is the distance from center to center of each electrode surface
%
% Configuration 1 - Monopolar electrode over the central line of the matrix
% Configuration 2 - Differential electrodes far away from the innervation
% zone row
% Configuration 3 - Differential electrodes far away from the central row 
% 
% The signal is simply the sum of the signal in each electrode inside the 
% virtual surface with dimensions of each type
%
% m is the iterator that runs through the columns of the matrix

tic
signal_load = load(handles.signal_file, 'emg_mono');
emg_mono = signal_load.emg_mono;
clear signal_load

% initializing variables
meps_mono = cell(1, handles.signal.n_conditions);
mep_mean_mono = cell(1, handles.signal.n_conditions);

handles.data.dantec_1_pp = zeros(handles.signal.n_conditions,1);
handles.data.dantec_2_pp = zeros(handles.signal.n_conditions,1);
handles.data.seniam_2_pp = zeros(handles.signal.n_conditions,1);
handles.data.dantec_3_pp = zeros(handles.signal.n_conditions,1);
handles.data.seniam_3_pp = zeros(handles.signal.n_conditions,1);

handles.data.dantec_1_fmed = zeros(handles.signal.n_conditions,1);
handles.data.dantec_2_fmed = zeros(handles.signal.n_conditions,1);
handles.data.seniam_2_fmed = zeros(handles.signal.n_conditions,1);
handles.data.dantec_3_fmed = zeros(handles.signal.n_conditions,1);
handles.data.seniam_3_fmed = zeros(handles.signal.n_conditions,1);

handles.data.dantec_1_rms = zeros(handles.signal.n_conditions,1);
handles.data.dantec_2_rms = zeros(handles.signal.n_conditions,1);
handles.data.seniam_2_rms = zeros(handles.signal.n_conditions,1);
handles.data.dantec_3_rms = zeros(handles.signal.n_conditions,1);
handles.data.seniam_3_rms = zeros(handles.signal.n_conditions,1);

sim_dantec_1 = zeros(1, handles.signal.n_conditions);
sim_dantec_2 = zeros(1, handles.signal.n_conditions);
sim_seniam_2 = zeros(1, handles.signal.n_conditions);
sim_dantec_3 = zeros(1, handles.signal.n_conditions);
sim_seniam_3 = zeros(1, handles.signal.n_conditions);
seniam_bipolar =  zeros(1, handles.signal.n_conditions);
seniam_monopolar =  zeros(1, handles.signal.n_conditions);

for i = 1:handles.signal.n_conditions % orientations

    sim_dantec_1(1:handles.data.s1(i)-handles.data.s0(i),i) = 0;
    sim_dantec_2(1:handles.data.s1(i)-handles.data.s0(i),i) = 0;
    sim_seniam_2(1:handles.data.s1(i)-handles.data.s0(i),i) = 0;
    sim_dantec_3(1:handles.data.s1(i)-handles.data.s0(i),i) = 0;
    sim_seniam_3(1:handles.data.s1(i)-handles.data.s0(i),i) = 0;
    seniam_bipolar(1:handles.data.s1(i)-handles.data.s0(i),i) = 0;
    seniam_monopolar(1:handles.data.s1(i)-handles.data.s0(i),i) = 0;

    for j = 1:handles.signal.n_channels % electrodes
        aux_channel = emg_mono{i}(:,j);
        aux_meps_mono = [];
        
        if sum(handles.data.trigger{i}(:,j)) ~= 0
            aux_trigger = handles.data.trigger{i}(:,j);
            for k = 1:length(aux_trigger) % number of stimulus
                meps_mono{i}{j}(:,k) = aux_channel(aux_trigger(k)+handles.data.s0(i):aux_trigger(k)+handles.data.s1(i)-1);
                if k == 1
                    aux_meps_mono = meps_mono{i}{j}(:,k);
                else
                    aux_meps_mono = cat(2, aux_meps_mono, meps_mono{i}{j}(:,k));
                end
            end
        else
            meps_mono{i}{j} = zeros(handles.data.s1(i)-handles.data.s0(i), 10);
        end
         mep_mean_mono{i}(:,j) = mean(meps_mono{i}{j},2);  
    end
    
    % sim 1 - monopolar dantec electrode over the central line
    aux_elec_0 = zeros(handles.data.s1(i)-handles.data.s0(i),1);
    for m = 0:4
        aux_elec_0 = aux_elec_0 + sum(mep_mean_mono{i}(:,(6 + 13*m:8 + 13*m)),2);
    end
    sim_dantec_1(:,i) = aux_elec_0;
    
    % sim 2 - bipolar electrodes positioned far away from the inervation
    % zone row
    % sim 2 - dantec
    aux_elec_1 = zeros(handles.data.s1(i)-handles.data.s0(i),1);
    aux_elec_2 = zeros(handles.data.s1(i)-handles.data.s0(i),1);
    for m = 0:4
        aux_elec_1 = aux_elec_1 + sum(mep_mean_mono{i}(:,(handles.data.iz_row-3 + 13*m:handles.data.iz_row-1 + 13*m)),2);
        aux_elec_2 = aux_elec_2 + sum(mep_mean_mono{i}(:,(handles.data.iz_row+1 + 13*m:handles.data.iz_row+3 + 13*m)),2);
    end
    sim_dantec_2(:,i) = aux_elec_1 - aux_elec_2;
    
    % sim 2 - seniam
    aux_elec_3 = zeros(handles.data.s1(i)-handles.data.s0(i),1);
    aux_elec_4 = zeros(handles.data.s1(i)-handles.data.s0(i),1);
    for m = 0:4
        aux_elec_3 = aux_elec_3 + sum(mep_mean_mono{i}(:,handles.data.iz_row-1 + 13*m),2);
        aux_elec_4 = aux_elec_4 + sum(mep_mean_mono{i}(:,handles.data.iz_row+1 + 13*m),2);
    end
    sim_seniam_2(:,i) = aux_elec_3 - aux_elec_4;
    
    clear aux_elec_0 aux_elec_1 aux_elec_2 aux_elec_3 aux_elec_4
    
    % sim 3 - bipolar electrodes positioned far away from the central row
    % sim 3 - dantec
    aux_elec_1 = zeros(handles.data.s1(i)-handles.data.s0(i),1);
    aux_elec_2 = zeros(handles.data.s1(i)-handles.data.s0(i),1);
    for m = 0:4
        aux_elec_1 = aux_elec_1 + sum(mep_mean_mono{i}(:,(4 + 13*m:6 + 13*m)),2);
        aux_elec_2 = aux_elec_2 + sum(mep_mean_mono{i}(:,(8 + 13*m:10 + 13*m)),2);
    end
    sim_dantec_3(:,i) = aux_elec_1 - aux_elec_2;
    
    % sim 3 - seniam
    aux_elec_3 = zeros(handles.data.s1(i)-handles.data.s0(i),1);
    aux_elec_4 = zeros(handles.data.s1(i)-handles.data.s0(i),1);
    for m = 0:4
        aux_elec_3 = aux_elec_3 + sum(mep_mean_mono{i}(:,6 + 13*m),2);
        aux_elec_4 = aux_elec_4 + sum(mep_mean_mono{i}(:,8 + 13*m),2);
    end
    sim_seniam_3(:,i) = aux_elec_3 - aux_elec_4;
    
    clear aux_elec_1 aux_elec_2 aux_elec_3 aux_elec_4
    
    % sim 4 - square electrodes 5 mm size and 9.6 mm IED
    % seniam_bipolar and seniam_monopolar were used to paper
    % seniam_monopolar is equal to sim_dantec_1
    aux_elec_5 = zeros(handles.data.s1(i)-handles.data.s0(i),1);
    aux_elec_6 = zeros(handles.data.s1(i)-handles.data.s0(i),1);
    for m = 1:3
        aux_elec_5 = aux_elec_5 + sum(mep_mean_mono{i}(:,(4 + 13*m:6 + 13*m)),2);
        aux_elec_6 = aux_elec_6 + sum(mep_mean_mono{i}(:,(8 + 13*m:10 + 13*m)),2);
    end
    seniam_bipolar(:,i) = aux_elec_5 - aux_elec_6;
    
    aux_elec_7 = zeros(handles.data.s1(i)-handles.data.s0(i),1);
    for n = 0:4
        aux_elec_7 = aux_elec_7 + sum(mep_mean_mono{i}(:,(6 + 13*n:8 + 13*n)),2);
    end
    seniam_monopolar(:,i) = aux_elec_7;
    
    clear aux_elec_5 aux_elec_6 aux_elec_7
    
    % signal peak to peak amplitude
    handles.data.dantec_1_pp(i) = max(sim_dantec_1(:,i)) - min(sim_dantec_1(:,i));
    handles.data.dantec_2_pp(i) = max(sim_dantec_2(:,i)) - min(sim_dantec_2(:,i));
    handles.data.seniam_2_pp(i) = max(sim_seniam_2(:,i)) - min(sim_seniam_2(:,i));
    handles.data.dantec_3_pp(i) = max(sim_dantec_3(:,i)) - min(sim_dantec_3(:,i));
    handles.data.seniam_3_pp(i) = max(sim_seniam_3(:,i)) - min(sim_seniam_3(:,i));
    
    % signal rms amplitude and median frequency
    [handles.data.dantec_1_fmed(i), handles.data.dantec_1_rms(i), ~] = Fmed3cla(sim_dantec_1(:,i),...
        handles.signal.fsample(i), handles.data.s1(i)-handles.data.s0(i));
    [handles.data.dantec_2_fmed(i), handles.data.dantec_2_rms(i), ~] = Fmed3cla(sim_dantec_2(:,i),...
        handles.signal.fsample(i), handles.data.s1(i)-handles.data.s0(i));
    [handles.data.seniam_2_fmed(i), handles.data.seniam_2_rms(i), ~] = Fmed3cla(sim_seniam_2(:,i),...
        handles.signal.fsample(i), handles.data.s1(i)-handles.data.s0(i));
    [handles.data.dantec_3_fmed(i), handles.data.dantec_3_rms(i), ~] = Fmed3cla(sim_dantec_3(:,i),...
        handles.signal.fsample(i), handles.data.s1(i)-handles.data.s0(i));
    [handles.data.seniam_3_fmed(i), handles.data.seniam_3_rms(i), ~] = Fmed3cla(sim_seniam_3(:,i),...
        handles.signal.fsample(i), handles.data.s1(i)-handles.data.s0(i));    
end

toc
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in button_ViewClear.
function button_ViewClear_Callback(hObject, eventdata, handles)
% hObject    handle to button_ViewClear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% TODO: This function does not clear the threshold lines. Fix it!
% % % this was the past function of clear, now it deletes all channels in one angle
% % % set current channel plots as invisible
% % handles = visibleoff(handles);

% set current channel plots as invisible
handles = visibleoff(handles);

for j = 1:handles.signal.n_channels
    axes(eval(strcat('handles.axes',num2str(j))))
    
    if ~isempty(handles.data.mepmax{handles.id_angle}{j})
        handles.data.amp_pp{handles.id_angle}{j} = handles.data.mepmax{handles.id_angle}{j} - handles.data.mepmin{handles.id_angle}{j};
        if handles.data.amp_pp{handles.id_angle}{j} == 0 || isnan(handles.data.amp_pp{handles.id_angle}{j})
            handles.data.amp_pp{handles.id_angle}{j} = NaN;
        end
    else
        handles.data.amp_pp{handles.id_angle}{j} = NaN;
    end
    if isfield(handles,'hmepmax')
        if ishandle(handles.hmepmax{handles.id_angle, j})
            delete(handles.hmepmax{handles.id_angle, j})
            handles.data.mepmax{handles.id_angle}{j} = [];
        end
    end
    if isfield(handles,'hmepmin')
        if ishandle(handles.hmepmin{handles.id_angle, j})
            delete(handles.hmepmin{handles.id_angle, j})
            handles.data.mepmin{handles.id_angle}{j} = [];
        end
    end
update_waitbar(handles,j/handles.signal.n_channels)
end

% updating checkbox values for MEP visualization
set(handles.checkbox_ViewMinMax,'Value',1)
set(handles.checkbox_ViewMEPsMean,'Value',1)

checkbox_ViewMinMax_Callback(hObject, eventdata, handles)
checkbox_ViewMEPsMean_Callback(hObject, eventdata, handles)

% Update handles structure
guidata(hObject, handles);


function out = refreshaxes(handles)

% position of electrode on array
pos = get(gca,'Tag');
pos(1:4)=[];
pos = str2num(pos);


if get(handles.checkbox_ViewSignal,'Value')
    % Open trigger manipulation window
    % I did no fix the Trigger Detect to this new version of MEP Hunter -
    % with handles.id_angle and etc
    [trigger handles.data] = TriggerDetect(handles.data,handles.data.trigger{handles.id_angle}(:,pos),handles.data.xs{handles.id_angle}, pos); 
    
    if ishandle(handles.htrigger{handles.id_angle, pos})
        delete(handles.htrigger{handles.id_angle, pos});
        delete(handles.hsignal{handles.id_angle, pos});
        delete(handles.hline{handles.id_angle, pos});
    end
    
    handles.data.trigger{handles.id_angle}(:,pos) = trigger;
    
    if sum(handles.data.trigger{handles.id_angle}(:,pos)) ~= 0
        axes(eval(strcat('handles.axes',num2str(pos))));
        handles.hsignal{handles.id_angle}(pos) = plot(handles.data.xs{handles.id_angle},handles.signal.emg_diff{handles.id_angle}(:,pos));
        % red circulus
%         if sum(handles.data.trigger(:,handles.id_angle)) ~= 0
%              handles.htrigger{handles.id_angle, pos} = plot(handles.data.xs{handles.id_angle}(handles.data.trigger{handles.id_angle}{pos}(:,1)),...
%                  handles.data.trigger{handles.id_angle}{pos}(:,2),'ro');
%         else
%             handles.htrigger{handles.id_angle, pos} = nan;
%         end
        % green line
        if sum(handles.data.trigger{handles.id_angle}(:,pos)) ~= 0
            handles.hline{handles.id_angle, pos} = line([handles.data.xs{handles.id_angle}(handles.data.trigger{handles.id_angle}(:,pos))',...
                handles.data.xs{handles.id_angle}(handles.data.trigger{handles.id_angle}(:,pos))'],...
                [min(handles.signal.emg_diff{handles.id_angle}(:,pos)),...
                max(handles.signal.emg_diff{handles.id_angle}(:,pos))],...
                'Color','g');
        else
            handles.hline{handles.id_angle, pos} = nan;
        end
    end
    
    out = handles;
    
else
    %Open MEP manipulation window to set min, max, latency and coordinate on matrix
       
    [handles.data handles.signal minmax latency] = MEPDetect_angles(handles.data, handles.signal,...
        handles.MEPStart, handles.signal_file,...
        handles.data.xs{handles.id_angle}(1:handles.data.s1(handles.id_angle)-handles.data.s0(handles.id_angle)),...
        [pos handles.id_angle]);
    
    handles.data.latency{handles.id_angle}{pos} = latency;
    
    % Plots of latencies
    if isfield(handles,'hlatencystart')
        if ishandle(handles.hlatencystart{handles.id_angle, pos})
            delete(handles.hlatencystart{handles.id_angle, pos})
        end
    end

    if isfield(handles,'hlatencystop')
        if ishandle(handles.hlatencystop{handles.id_angle, pos})
            delete(handles.hlatencystop{handles.id_angle, pos})
        end
    end

    axes(gca)
    
    % Plot of MEPs Mult
    colors = [0.0117 0.0343 0.894; 0.0117 0.6078 0.894;...
            0.0117 0.89411 0.69411; 0.0117 0.8902 0.3255;...
            0.0196 0.8901 0.0196; 0.6274 0.8902 0.3255;...
            0.9804 0.9568 0.1137; 0.9804 0.6941 0.1137;...
            0.9804 0.3294 0.1137; 0.8627 0.0196 0.0196;...
            0.0 0.0 0.0; 0.0 0.0 0.0;...
            0.0 0.0 0.0; 0.0 0.0 0.0;...
            0.0 0.0 0.0; 0.0 0.0 0.0];
        
    if isfield(handles,'hmeps')
        if ishandle(handles.hmeps{handles.id_angle, pos})
            delete(handles.hmeps{handles.id_angle, pos})
        end
    end   
    if sum(handles.data.trigger{handles.id_angle}(:,pos)) ~= 0
        aux_trigger = handles.data.trigger{handles.id_angle}(:,pos);
    else
        aux_trigger = [];
        colors = [];
    end
    
    % Number of stimulus
    for k = 1:length(aux_trigger)
        handles.hmeps{handles.id_angle, pos}(k) = plot(handles.data.xs{handles.id_angle}(1:handles.data.s1(handles.id_angle)-handles.data.s0(handles.id_angle))+handles.MEPStart,...
            handles.data.meps_diff{handles.id_angle}{pos}{k},'-','color',colors(k,:),'Visible','off');
    end

    mepsvalue = get(handles.checkbox_ViewMEPsMult,'Value');
    if mepsvalue == 1
        if isfield(handles,'hmeps') == 1
            if ishandle(handles.hmeps{handles.id_angle, pos})
                set(handles.hmeps{handles.id_angle, pos},'Visible','on')
            end
        end
    end
    
    % Plot of MEP Mean
    if isfield(handles,'hmepmean')
        if ishandle(handles.hmepmean{handles.id_angle, pos})
            delete(handles.hmepmean{handles.id_angle, pos})
        end
    end
    mepmeanvalue = get(handles.checkbox_ViewMEPsMean,'Value');
    handles.hmepmean{handles.id_angle, pos} = plot(handles.data.xs{handles.id_angle}(1:handles.data.s1(handles.id_angle)-handles.data.s0(handles.id_angle))+handles.MEPStart,...
        handles.data.mepmean{handles.id_angle}{pos}, 'Visible','off');
    if mepmeanvalue == 1
        if isfield(handles,'hmepmean')
            if ishandle(handles.hmepmean{handles.id_angle, pos})
                set(handles.hmepmean{handles.id_angle, pos},'Visible','on')
            end
        end
    end
      
    % Plots of minmax
    if isfield(handles,'hmepmax')
        if ishandle(handles.hmepmax{handles.id_angle, pos})
            delete(handles.hmepmax{handles.id_angle, pos})
        end
    end
    
    if isfield(handles,'hmepmin')
        if ishandle(handles.hmepmin{handles.id_angle, pos})
            delete(handles.hmepmin{handles.id_angle, pos})
        end
    end
    
    if isempty(minmax)
        handles.data.mepmax{handles.id_angle}{pos} = [];
        handles.data.mepmin{handles.id_angle}{pos} = [];
        
    else
        if minmax ~= 0
            handles.data.mepmax{handles.id_angle}{pos} = max(minmax(:,2));
            handles.data.mepmin{handles.id_angle}{pos} = min(minmax(:,2));
            
            aux_pos_max = minmax(:,2) == handles.data.mepmax{handles.id_angle}{pos};
            pos_max = minmax(aux_pos_max,1);
            
            aux_pos_min = minmax(:,2) == handles.data.mepmin{handles.id_angle}{pos};
            pos_min = minmax(aux_pos_min,1);
            
            handles.hmepmax{handles.id_angle, pos} = plot(pos_max + handles.MEPStart,...
                handles.data.mepmax{handles.id_angle}{pos},'+r','Visible','off');
            
            handles.hmepmin{handles.id_angle, pos} = plot(pos_min + handles.MEPStart,...
                handles.data.mepmin{handles.id_angle}{pos},'+r','Visible','off');
            
            mepminmaxvalue = get(handles.checkbox_ViewMinMax,'Value');
            if mepminmaxvalue
                if isfield(handles,'hmepmax')
                    if ishandle(handles.hmepmax{handles.id_angle, pos}) && ishandle(handles.hmepmin{handles.id_angle, pos})
                        set(handles.hmepmax{handles.id_angle, pos},'Visible','on')
                        set(handles.hmepmin{handles.id_angle, pos},'Visible','on')
                    end
                end
            end
        end
    end
    
    if latency ~= 0
        handles.data.duration{handles.id_angle}{pos} = handles.data.latency{handles.id_angle}{pos}(2,1) - handles.data.latency{handles.id_angle}{pos}(1,1);
        handles.hlatencystart{handles.id_angle, pos} = plot(handles.data.latency{handles.id_angle}{pos}(1,1),...
            handles.data.latency{handles.id_angle}{pos}(1,2),'gv','Visible','off');
        handles.hlatencystop{handles.id_angle, pos} = plot(handles.data.latency{handles.id_angle}{pos}(2,1),...
            handles.data.latency{handles.id_angle}{pos}(2,2),'r^','Visible','off');
        latencyvalue = get(handles.checkbox_ViewLatency,'Value');
        if latencyvalue
            if isfield(handles,'hlatencystart')
                if ishandle(handles.hlatencystart{handles.id_angle, pos})
                    set(handles.hlatencystart{handles.id_angle, pos},'Visible','on')
                    set(handles.hlatencystop{handles.id_angle, pos},'Visible','on')
                end
            end
        end
    end
    
    if isempty(aux_trigger)
        handles.hmeps{handles.id_angle, pos} = nan;
        handles.hmepmean{handles.id_angle, pos} = nan;
        handles.hmepmax{handles.id_angle, pos} = nan;
        handles.hmepmin{handles.id_angle, pos} = nan;
    end
       
    % refreshing peak to peak map
    if ~isempty(handles.data.mepmax{handles.id_angle}{pos})
        handles.data.amp_pp{handles.id_angle}{pos} = handles.data.mepmax{handles.id_angle}{pos} - handles.data.mepmin{handles.id_angle}{pos};
        if handles.data.amp_pp{handles.id_angle}{pos} == 0 || isnan(handles.data.amp_pp{handles.id_angle}{pos})
            handles.data.amp_pp{handles.id_angle}{pos} = NaN;
            handles.data.amp_rms{handles.id_angle}{pos} = NaN;
            handles.data.offset_rms{handles.id_angle}{pos} = NaN;
            handles.data.fmed{handles.id_angle}{pos} = NaN;
        end
    else
        handles.data.amp_pp{handles.id_angle}{pos} = NaN;
        handles.data.amp_rms{handles.id_angle}{pos} = NaN;
        handles.data.offset_rms{handles.id_angle}{pos} = NaN;
        handles.data.fmed{handles.id_angle}{pos} = NaN;
    end
    
    out = handles;
end


function update_waitbar(handles,value)
% Update waitbar
h=handles.axes_waitbar;
if ~ishandle(h)
    return
end
set(h,'Visible','On');
axes(h);
cla;
patch([0,value,value,0],[0,0,1,1],'b');
axis([0,1,0,1]);
axis off;


function out = visibleoff(handles)
% set current channel plots as invisible

for i = 1:handles.signal.n_conditions % EMG mode
    for j = 1:handles.signal.n_channels % Electrodes
        
        if isfield(handles,'hampthreshold')
            if ishandle(handles.hampthreshold{i, j})
                set(handles.hampthreshold{i, j},'Visible','off')
            end
        end
        
        if isfield(handles,'hsignal')
            if ishandle(handles.hsignal{i, j})
                set(handles.hsignal{i, j},'Visible','off')
            end
        end
        
        if isfield(handles,'htrigger')
            if ishandle(handles.htrigger{i, j})
                set(handles.htrigger{i, j},'Visible','off')
            end
        end
        
        if isfield(handles,'hline')
            if ishandle(handles.hline{i, j})
                set(handles.hline{i, j},'Visible','off')
            end
        end
        
        if isfield(handles,'hmeps')
            if ishandle(handles.hmeps{i, j})
                set(handles.hmeps{i, j},'Visible','off')
            end
        end
        
        if isfield(handles,'hmepmean')
            if ishandle(handles.hmepmean{i, j})
                set(handles.hmepmean{i, j},'Visible','off')
            end
        end
        
        if isfield(handles,'hmepmax')
            if ishandle(handles.hmepmax{i, j})
                set(handles.hmepmax{i, j},'Visible','off')
            end
        end
        
        if isfield(handles,'hmepmin')
            if ishandle(handles.hmepmin{i, j})
                set(handles.hmepmin{i, j},'Visible','off')
            end
        end
        
        if isfield(handles,'hlatencystart')
            if ishandle(handles.hlatencystart{i, j})
                set(handles.hlatencystart{i, j},'Visible','off')
            end
        end
        
        if isfield(handles,'hlatencystop')
            if ishandle(handles.hlatencystop{i, j})
                set(handles.hlatencystop{i, j},'Visible','off')
            end
        end
    end
    update_waitbar(handles,j/(handles.signal.n_channels))
end

set(handles.checkbox_AmplitudeThreshold,'Value',0)
set(handles.checkbox_ViewLatency,'Value',0)
set(handles.checkbox_ViewLine,'Value',0)
set(handles.checkbox_ViewMEPsMean,'Value',0)
set(handles.checkbox_ViewMEPsMult,'Value',0)
set(handles.checkbox_ViewMinMax,'Value',0)
set(handles.checkbox_ViewSignal,'Value',0)
set(handles.checkbox_ViewTrigger,'Value',0)

handles.signalvalue = 0;
handles.mepmeanvalue = 0;
handles.mepsvalue = 0;
handles.linevalue = 0;
handles.triggervalue = 0;
handles.mepminmaxvalue = 0;
handles.ampthresholdvalue = 0;
handles.latencyvalue = 0;

out = handles;

function out = meps_calculation(handles)

% this is a map of colors similar to the jet
% the blue is referring to the first stimuli and the
% red is used for the last stimuli.
aux_colors = [0.0117 0.0343 0.894; 0.0117 0.6078 0.894;...
    0.0117 0.89411 0.69411; 0.0117 0.8902 0.3255;...
    0.0196 0.8901 0.0196; 0.6274 0.8902 0.3255;...
    0.9804 0.9568 0.1137; 0.9804 0.6941 0.1137;...
    0.9804 0.3294 0.1137; 0.8627 0.0196 0.0196;...
    0.0 0.0 0.0; 0.0 0.0 0.0;...
    0.0 0.0 0.0; 0.0 0.0 0.0;...
    0.0 0.0 0.0; 0.0 0.0 0.0];

for j = 1:handles.signal.n_channels % Electrodes
    axes(eval(strcat('handles.axes',num2str(j))))
    hold on
    for i = 1:handles.signal.n_conditions % Angles
        colors = aux_colors;
        if sum(handles.data.trigger{i}(:,j)) ~= 0
            aux_trigger = handles.data.trigger{i}(:,j);
        else
            colors = [];
            aux_trigger = [];
        end
        
        aux_channel = handles.signal.emg_diff{i}(:,j);
        aux_meps = [];
        aux_meps_column = [];
        aux_offset_signal_columns = [];
        aux_offset_signal = [];
        
        for k = 1:length(aux_trigger) % Number of stimulus
            handles.data.meps_diff{i}{j}{k} = aux_channel(aux_trigger(k)+handles.data.s0(i):aux_trigger(k)+handles.data.s1(i)-1);
            handles.data.offset_signals{i}{j}{k} = aux_channel(aux_trigger(k)-handles.data.s1(i)-30+1:aux_trigger(k)-handles.data.s0(i)-30);
            handles.hmeps{i, j}(k) = plot(handles.data.xs{i}(1:handles.data.s1(i)-handles.data.s0(i))+handles.MEPStart,...
                handles.data.meps_diff{i}{j}{k},'-','color',colors(k,:),'Visible', 'off');
            if k == 1
                aux_meps_column = handles.data.meps_diff{i}{j}{k};
                aux_meps = handles.data.meps_diff{i}{j}{k};
                aux_offset_signal_columns = handles.data.offset_signals{i}{j}{k};
                aux_offset_signal = handles.data.offset_signals{i}{j}{k};
            else
                aux_meps_column = cat(1, aux_meps_column, handles.data.meps_diff{i}{j}{k});
                aux_meps = cat(2, aux_meps, handles.data.meps_diff{i}{j}{k});
                aux_offset_signal_columns = cat(1, aux_offset_signal_columns,...
                    handles.data.offset_signals{i}{j}{k});
                aux_offset_signal = cat(2, aux_offset_signal,...
                    handles.data.offset_signals{i}{j}{k});
            end
        end
        
        if isempty(aux_trigger)
            handles.data.meps_diff{i}{j} = [];
            handles.data.offset_signals{i}{j} = [];
        end
        
        % median frequency, rms values and mean frequency for meps
        [handles.data.fmed{i}{j}, handles.data.amp_rms{i}{j}, ~] = Fmed3cla(aux_meps_column,...
            handles.signal.fsample(i), handles.data.s1(i)-handles.data.s0(i));
        handles.data.fmed{i}{j} = mean(handles.data.fmed{i}{j}, 1);
        handles.data.amp_rms{i}{j} = mean(handles.data.amp_rms{i}{j}, 1);
        
        % median frequency and rms value for the offset signal
        [handles.data.offset_fmed{i}{j}, handles.data.offset_rms{i}{j}, ~] = Fmed3cla(aux_offset_signal_columns,...
            handles.signal.fsample(i), handles.data.s1(i)-handles.data.s0(i));
        handles.data.offset_rms{i}{j} = mean(handles.data.offset_rms{i}{j}, 1);
        handles.data.offset_fmed{i}{j} = mean(handles.data.offset_fmed{i}{j}, 1);
        
                       
        % mean values, max, min and their backups
%         handles.data.meps_diff_bkp = handles.data.meps_diff;
%         handles.data.mepmean{i}{j} = mean(aux_meps,2);
%         handles.data.mepmean_bkp = handles.data.mepmean;
%         handles.data.mepmax{i}{j} = max(handles.data.mepmean{i}{j});
%         handles.data.mepmax_bkp = handles.data.mepmax;
%         handles.data.mepmin{i}{j} = min(handles.data.mepmean{i}{j});
%         handles.data.mepmin_bkp = handles.data.mepmin;
        handles.data.meps_diff_bkp = handles.data.meps_diff;
        handles.data.mepmean{i}{j} = mean(aux_meps,2);
        handles.data.mepmean_bkp = handles.data.mepmean;
        
        if ~isempty(aux_trigger)
            % find mep peak
            Fs = handles.signal.fsample(i);
            mepwindow = handles.data.mepmean{i}{j}(round((handles.MEPStart + 5)*Fs/1000):round((handles.MEPStart + 25)*Fs/1000));
            L = length(mepwindow);              % Length of signal
            t = handles.data.xs{i}(round((handles.MEPStart + 5)*Fs/1000)) + (0:L-1)/Fs;
            [pks, plocs] = findpeaks(mepwindow);
            %         mep_peak = t(plocs(1)); % instant of mep peak
            pp = pks == max(pks);
            
            % find mep valley
            [vls, vlocs] = findpeaks(-mepwindow);
            %         mep_valley = t(vlocs(1)); % instant of mep peak
            vv = vls == max(vls);
            
            if ~isempty(plocs) && ~isempty(vlocs)
                handles.data.mepmax{i}{j} = mepwindow(plocs(pp));
                handles.data.mepmin{i}{j} = mepwindow(vlocs(vv));
                handles.data.pos_max{i}{j} = find(handles.data.xs{i} == t(plocs(pp)));
                handles.data.pos_min{i}{j} = find(handles.data.xs{i} == t(vlocs(vv)));
            else
                handles.data.mepmax{i}{j} = [];
                handles.data.mepmin{i}{j} = [];
                handles.data.pos_max{i}{j} = [];
                handles.data.pos_min{i}{j} = [];
            end
                        
        else
            handles.data.mepmax{i}{j} = [];
            handles.data.mepmin{i}{j} = [];
            handles.data.pos_max{i}{j} = [];
            handles.data.pos_min{i}{j} = [];
        end
        
        handles.data.mepmax_bkp = handles.data.mepmax;
        handles.data.mepmin_bkp = handles.data.mepmin;
        % mean values, max, min and their backups
        handles.data.fmed_bkp = handles.data.fmed;
        handles.data.amp_rms_bkp = handles.data.amp_rms;
        handles.data.offset_rms_bkp = handles.data.offset_rms;
        handles.data.offset_fmed_bkp = handles.data.offset_fmed;
        
%         % position in x axis of min and max values
%         handles.data.pos_max{i}{j} = find(handles.data.mepmean{i}{j} == handles.data.mepmax{i}{j});
%         handles.data.pos_min{i}{j} = find(handles.data.mepmean{i}{j} == handles.data.mepmin{i}{j});
%         handles.data.pos_max_init = handles.data.pos_max;
%         handles.data.pos_min_init = handles.data.pos_min;
        handles.data.pos_max_init = handles.data.pos_max;
        handles.data.pos_min_init = handles.data.pos_min;
        % plots
        if isempty(aux_trigger)
            handles.hmeps{i, j} = nan;
            handles.hmepmean{i, j} = nan;
            handles.hmepmax{i, j} = nan;
            handles.hmepmin{i, j} = nan;
        else
            handles.hmepmean{i, j} = plot(handles.data.xs{i}(1:handles.data.s1(i)-handles.data.s0(i))+handles.MEPStart,...
                handles.data.mepmean{i}{j},'Visible','off');
            handles.hmepmax{i, j} = plot(handles.data.xs{i}(handles.data.pos_max{i}{1, j})+handles.MEPStart,...
                handles.data.mepmax{i}{j},'+r','Visible','off');
            handles.hmepmin{i, j} = plot(handles.data.xs{i}(handles.data.pos_min{i}{1, j})+handles.MEPStart,...
                handles.data.mepmin{i}{j},'+r','Visible','off');
        end
        
    end
    update_waitbar(handles,j/(handles.signal.n_channels))
end

out = handles;

function out = trigger_calculation(handles)

% ad_range - of the emg equipment
% conv_uv - multiplication factor to see signal in microVolts
% TODO: put a text box in the UI for the user to change the conversion
% factor
ad_range = 5.0;
conv_uv = 1000000.0;

handles.hp_filter = str2double(get(handles.edit_HPF,'String'));
handles.lp_filter = str2double(get(handles.edit_LPF,'String'));

samples_triggeron = cell(1, handles.signal.n_conditions);
% handles.data.trigger = cell(handles.signal.n_conditions);

if handles.hp_filter == 25 && handles.lp_filter == 400
    % finding samples denoting trigger onset
    trigger_load = load(handles.signal_file, 'trigger_aux');
    handles.signal.trigger_aux = trigger_load.trigger_aux;
%     handles.signal = rmfield(handles.signal, 'trigger_aux');
    clear trigger_load
else
    signal_load = load(handles.signal_file, 'raw_data');
    handles.signal.raw_data = signal_load.raw_data;
    clear signal_load
end

emg_data_aux = cell(1, handles.signal.n_conditions);
emg_mono_aux = cell(1, handles.signal.n_conditions);
emg_diff_aux = cell(1, handles.signal.n_conditions);

for i = 1:handles.signal.n_conditions
    
    if handles.hp_filter ~= 25 && handles.lp_filter ~= 400
        % Creating new signal with band pass filter specified by the user
        [b,a] = butter(2,[handles.hp_filter handles.lp_filter]*2/handles.signal.fsample(i));
        
        emg_data_aux{i} = handles.signal.raw_data{i}(:,...
            [8 9 19 20 27 30 40 41 51 52 61 1 7 10 18 21 28 31 39 42 50 53 60 62 2 6 11 16 22 29 32 38 43 49 54 59 63 3 5 12 15 23 26 34 37 44 47 55 58 64 4 13 14 24 25 35 36 45 46 56 57]);
        % Channel 16 does not have signel - this line is an interpolation
        emg_data_aux{i}(:, 16) = mean(emg_data_aux{i}(:,[3 4 5 15 17 27 30 31]),2);
        
        emg_data_aux{i} = filtfilt(b,a,...
            (emg_data_aux{i}*ad_range*conv_uv)/(handles.signal.signal_gain(i)*2^(handles.signal.ad_bits)));
        
        emg_mono_aux{i} = [nan(length(emg_data_aux{i}), 1) emg_data_aux{i}(:, 1:11)...
            nan(length(emg_data_aux{i}), 1) emg_data_aux{i}(:, 12:50) nan(length(emg_data_aux{i}), 1)...
            emg_data_aux{i}(:, 51:end) nan(length(emg_data_aux{i}), 1)];
        
        emg_diff_aux{i} = diff(emg_data_aux{i},1,2);        
        emg_diff_aux{i} = [nan(length(emg_diff_aux{i}), 1) emg_diff_aux{i}(:, 1:10)...
            nan(length(emg_diff_aux{i}), 1) nan(length(emg_diff_aux{i}), 1)...
            emg_diff_aux{i}(:, 12:23) nan(length(emg_diff_aux{i}), 1)...
            emg_diff_aux{i}(:, 25:36) nan(length(emg_diff_aux{i}), 1)...
            emg_diff_aux{i}(:, 38:49) nan(length(emg_diff_aux{i}), 1)...
            nan(length(emg_diff_aux{i}), 1) emg_diff_aux{i}(:, 51:end)...
            nan(length(emg_diff_aux{i}), 1) nan(length(emg_diff_aux{i}), 1)];
                      
        % Triggering EMGs
%         handles.signal.trigger_aux{i} = handles.signal.raw_data{i}(:,65)*5/(2^(handles.signal.ad_bits-1));
    end
    
    % finding samples denoting trigger onset
    samples_triggeron{i} = find(handles.signal.trigger_aux{i}>4);
    samples_triggeron{i} = (samples_triggeron{i}(diff([-inf;samples_triggeron{i}])>1));
    
end

if handles.hp_filter ~= 25 && handles.lp_filter ~= 400
    save('signal_temp2', '-struct', 'signal_aux');
    handles.signal = rmfield(handles.signal, {'trigger_aux',...
        'raw_data', 'emg_mono', 'emg_diff'});
end
clear emg_data_aux emg_mono_aux emg_diff_aux

% trigger computation based no HDsEMG aux channel with trigger signal
if handles.data.slope == 0
    for j = 1:handles.signal.n_channels % Electrodes
        axes(eval(strcat('handles.axes',num2str(j))))
        hold on
        for i = 1:handles.signal.n_conditions % Angles
            
            if sum(isnan(handles.signal.emg_diff{i}(samples_triggeron{i},j)))==0
                handles.data.trigger{i}(:,j) = samples_triggeron{i}(:);
            else
                handles.data.trigger{i}(:,j) = zeros(size(samples_triggeron{i},1),1);
            end
            % signal plot
            handles.hsignal{i, j} = plot(handles.data.xs{i},handles.signal.emg_diff{i}(:,j),...
                'Visible','off');
            
%             % red circulus plot
%             if sum(handles.data.trigger{i}(:,j)) ~= 0
%                 handles.htrigger{i, j} = plot(handles.data.xs{i}(handles.data.trigger{i}{j}(:,1)),...
%                     handles.data.trigger{i}{j}(:,2),'ro', 'Visible', 'off');
%             else
%                 handles.htrigger{i, j} = nan;
%             end

            % green line plot
            if sum(handles.data.trigger{i}(:,j))~=0
                a = line([handles.data.xs{i}(handles.data.trigger{i}(:,j))',...
                    handles.data.xs{i}(handles.data.trigger{i}(:,j))'],...
                    [min(handles.signal.emg_diff{i}(:,j)),...
                    max(handles.signal.emg_diff{i}(:,j))],'Color','g','Visible','off');
                handles.hline{i, j} = a;
            else
                handles.hline{i, j} = nan;
            end
        end        
        set(eval(strcat('handles.axes',num2str(j))),'HitTest','on')
        update_waitbar(handles,j/handles.signal.n_channels)
    end
    
% trigger computation based on derivative properties
% this part is not fixed for the actual code
else
    for i = 1:length(handles.signal.emg_diff) 
        axes(eval(strcat('handles.axes',num2str(i))))           
        handles.data.trigger{i} = BiopacTrigger(handles.signal.emg_diff{i}(:,handles.id_angle),handles.data.slope);
        
        % red circulus  
        if ~isempty(handles.data.trigger{i})
            handles.htrigger(i) = plot(handles.data.xs{handles.id_angle}(handles.data.trigger{i}(:,1)),handles.data.trigger{i}(:,2),'ro', 'Visible', 'off');
        else
            handles.htrigger(i) = nan;
        end
        set(eval(strcat('handles.axes',num2str(i))),'HitTest','on')
        
        % green line
        if ~isempty(handles.data.trigger{i})
            a = line([handles.data.xs{handles.id_angle}{i}(handles.data.trigger(:,i))',handles.data.xs{handles.id_angle}{i}(handles.data.trigger(:,i))'],...
                [min(handles.signal.emg_diff{i}(:,handles.id_angle)),...
                max(handles.signal.emg_diff{i}(:,handles.id_angle))],...
                'Color','g','Visible','off');
            handles.hline{i} = a;
        else
            handles.hline{i} = nan;
        end
        update_waitbar(handles,i/length(handles.signal.emg_diff))
    end
end

out = handles;


% --- Executes during object creation, after setting all properties.
function edit_AmplitudeThresholdDiff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_AmplitudeThresholdDiff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function edit_AmplitudeThresholdMono_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_AmplitudeThresholdMono (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function edit_AmpThre_Diff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_AmplitudeThresholdDiff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function edit_AmpThre_Mono_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_AmplitudeThresholdMono (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function edit_HPF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_HPF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function edit_LPF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_LPF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function edit_MEPEnd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_MEPEnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function edit_MEPStart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_MEPStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function edit_Threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_AmplitudeThresholdDiff_Callback(hObject, eventdata, handles)
% hObject    handle to edit_AmplitudeThresholdDiff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_AmplitudeThresholdDiff as text
%        str2double(get(hObject,'String')) returns contents of edit_AmplitudeThresholdDiff as a double

function edit_AmplitudeThresholdMono_Callback(hObject, eventdata, handles)
% hObject    handle to edit_AmplitudeThresholdMono (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_AmplitudeThresholdMono as text
%        str2double(get(hObject,'String')) returns contents of edit_AmplitudeThresholdMono as a double

function edit_AmpThre_Diff_Callback(hObject, eventdata, handles)
% hObject    handle to edit_AmplitudeThresholdDiff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_AmplitudeThresholdDiff as text
%        str2double(get(hObject,'String')) returns contents of edit_AmplitudeThresholdDiff as a double

function edit_AmpThre_Mono_Callback(hObject, eventdata, handles)
% hObject    handle to edit_AmplitudeThresholdMono (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_AmplitudeThresholdMono as text
%        str2double(get(hObject,'String')) returns contents of edit_AmplitudeThresholdMono as a double

function edit_hemisphere_side_Callback(hObject, eventdata, handles)
% hObject    handle to edit_hemisphere_side (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_hemisphere_side as text
%        str2double(get(hObject,'String')) returns contents of edit_hemisphere_side as a double

handles.signal.emg_side = get(handles.edit_hemisphere_side,'String');

% Update handles structure
guidata(hObject, handles);

function edit_HPF_Callback(hObject, eventdata, handles)
% hObject    handle to edit_HPF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_HPF as text
%        str2double(get(hObject,'String')) returns contents of edit_HPF as a double

function edit_id_number_Callback(hObject, eventdata, handles)
% hObject    handle to edit_id_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_id_number as text
%        str2double(get(hObject,'String')) returns contents of edit_id_number as a double

handles.signal.patient_id = str2double(get(handles.edit_id_number,'String'));

% Update handles structure
guidata(hObject, handles);

function edit_LPF_Callback(hObject, eventdata, handles)
% hObject    handle to edit_LPF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_LPF as text
%        str2double(get(hObject,'String')) returns contents of edit_LPF as a double

function edit_MEPEnd_Callback(hObject, eventdata, handles)
% hObject    handle to edit_MEPEnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_MEPEnd as text
%        str2double(get(hObject,'String')) returns contents of edit_MEPEnd as a double

function edit_MEPStart_Callback(hObject, eventdata, handles)
% hObject    handle to edit_MEPStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_MEPStart as text
% str2double(get(hObject,'String')) returns contents of edit_MEPStart as a double

function edit_number_angle_stim_Callback(hObject, eventdata, handles)
% hObject    handle to edit_number_angle_stim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.signal.angle_stim(handles.id_angle) = str2double(get(handles.edit_number_angle_stim,'String'));

angles_names = '';

% Generate latency cell and initial total amplitude
for i = 1:handles.signal.n_conditions
    angles_names = strcat(angles_names,num2str(i),'-',...
        num2str(handles.signal.angle_stim(i)),'_');
end

set(handles.text_ChannelLegend, 'String',angles_names)
set(handles.text_Channel, 'String',num2str(handles.signal.angle_stim(handles.id_angle)))

% Update handles structure
guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit_number_angle_stim as text
%        str2double(get(hObject,'String')) returns contents of edit_number_angle_stim as a double

function edit_Threshold_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Threshold as text
%        str2double(get(hObject,'String')) returns contents of edit_Threshold as a double


% --- Executes on mouse press over axes background.
function axes1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function axes2_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);



% --- Executes on mouse press over axes background.
function axes3_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);



% --- Executes on mouse press over axes background.
function axes4_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);



% --- Executes on mouse press over axes background.
function axes5_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);




% --- Executes on mouse press over axes background.
function axes6_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);



% --- Executes on mouse press over axes background.
function axes7_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);



% --- Executes on mouse press over axes background.
function axes8_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);



% --- Executes on mouse press over axes background.
function axes9_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);



% --- Executes on mouse press over axes background.
function axes10_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);




% --- Executes on mouse press over axes background.
function axes11_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);



% --- Executes on mouse press over axes background.
function axes12_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);



% --- Executes on mouse press over axes background.
function axes13_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);



% --- Executes on mouse press over axes background.
function axes14_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);



% --- Executes on mouse press over axes background.
function axes15_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);



% --- Executes on mouse press over axes background.
function axes16_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);



% --- Executes on mouse press over axes background.
function axes17_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);



% --- Executes on mouse press over axes background.
function axes18_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);



% --- Executes on mouse press over axes background.
function axes19_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);



% --- Executes on mouse press over axes background.
function axes20_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);



% --- Executes on mouse press over axes background.
function axes21_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);



% --- Executes on mouse press over axes background.
function axes22_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);



% --- Executes on mouse press over axes background.
function axes23_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);



% --- Executes on mouse press over axes background.
function axes24_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);



% --- Executes on mouse press over axes background.
function axes25_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function axes26_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function axes27_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function axes28_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function axes29_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function axes30_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function axes31_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function axes32_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function axes33_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function axes34_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function axes35_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function axes36_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function axes37_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function axes38_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function axes39_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes39 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function axes40_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes40 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function axes41_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes41 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function axes42_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes42 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function axes43_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes43 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function axes44_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes44 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function axes45_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes45 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function axes46_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes46 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function axes47_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes47 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function axes48_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes48 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function axes49_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes49 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function axes50_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes50 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function axes51_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes51 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function axes52_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes52 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function axes53_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes53 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function axes54_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes54 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function axes55_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes55 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function axes56_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes56 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function axes57_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes57 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function axes58_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes58 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function axes59_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes59 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function axes60_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes60 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function axes61_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes61 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function axes62_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes62 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function axes63_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes63 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function axes64_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes64 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes on mouse press over axes background.
function axes65_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes65 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = refreshaxes(handles);

% Update handles structure
guidata(hObject, handles);
