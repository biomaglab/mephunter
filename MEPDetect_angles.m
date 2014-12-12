function varargout = MEPDetect_angles(varargin)
% MEPDETECT_ANGLES MATLAB code for MEPDetect_angles.fig
%      MEPDETECT_ANGLES, by itself, creates a new MEPDETECT_ANGLES or raises the existing
%      singleton*.
%
%      H = MEPDETECT_ANGLES returns the handle to a new MEPDETECT_ANGLES or the handle to
%      the existing singleton*.
%
%      MEPDETECT_ANGLES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MEPDETECT_ANGLES.M with the given input arguments.
%
%      MEPDETECT_ANGLES('Property','Value',...) creates a new MEPDETECT_ANGLES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MEPDetect_angles_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MEPDetect_angles_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MEPDetect_angles

% Last Modified by GUIDE v2.5 01-Oct-2013 17:26:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MEPDetect_angles_OpeningFcn, ...
                   'gui_OutputFcn',  @MEPDetect_angles_OutputFcn, ...
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


% --- Executes just before MEPDetect_angles is made visible.
function MEPDetect_angles_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MEPDetect_angles (see VARARGIN)


if ~isempty (varargin)
    handles.data = varargin{1};
    handles.signal = varargin{2};
    handles.MEPStart = varargin{3};
    handles.signal_file = varargin{4};
    handles.xs = varargin{5};
    handles.pos = varargin{6}(1);
    handles.id_angle = varargin{6}(2);
        
    axes(handles.axes1);
    hold on   
end

% get current values of min max positions and plot them
x_max_init =  handles.data.pos_max{handles.id_angle}{handles.pos}/handles.signal.fsample(handles.id_angle);
x_min_init =  handles.data.pos_min{handles.id_angle}{handles.pos}/handles.signal.fsample(handles.id_angle);

handles.hdata = plot(handles.xs+handles.MEPStart,handles.data.mepmean{handles.id_angle}{handles.pos});
if ~isempty(handles.data.mepmax{handles.id_angle}{handles.pos}) && ~isempty(handles.data.mepmin{handles.id_angle}{handles.pos})
    if handles.data.mepmax{handles.id_angle}{handles.pos} ~= 0 && handles.data.mepmin{handles.id_angle}{handles.pos} ~= 0
        handles.hmepmax_init = plot(x_max_init + handles.MEPStart, handles.data.mepmax{handles.id_angle}{handles.pos},'+g');
        handles.hmepmin_init = plot(x_min_init + handles.MEPStart, handles.data.mepmin{handles.id_angle}{handles.pos},'+g');
    end
end
% -------------------------------------

% Choose default command line output for MEPDetect_angles
handles.minmax = [x_min_init handles.data.mepmin{handles.id_angle}{handles.pos};...
    x_max_init handles.data.mepmax{handles.id_angle}{handles.pos}];

handles.output_data = handles.data;
handles.output_signal = handles.signal;
handles.output_minmax = handles.minmax;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MEPDetect_angles wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MEPDetect_angles_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiwait(hObject);
handles = guidata(hObject);

delete(hObject);

% Get default command line output from handles structure
varargout{1} = handles.output_data;
varargout{2} = handles.output_signal;
varargout{3} = handles.output_minmax;
varargout{4} = handles.output_latency;



% --- Executes on button press in button_initial_values.
function button_initial_values_Callback(hObject, eventdata, handles)
% hObject    handle to button_initial_values (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.data.count(handles.id_angle) = handles.data.count(handles.id_angle) - 1;

handles.data.outliers_mep{handles.id_angle}(handles.data.outliers_mep{handles.id_angle} == handles.pos)=[];
handles.data.ch_excluded{handles.id_angle}(handles.data.ch_excluded{handles.id_angle} == handles.pos)=[];

if isfield(handles,'hdata')
    if ishandle(handles.hdata)
        delete(handles.hdata)
    end
end
if isfield(handles,'hmepmax') && isfield(handles,'hmepmin')
    if ishandle(handles.hmepmax{handles.id_angle, handles.pos}) && ishandle(handles.hmepmin{handles.id_angle, handles.pos})
        delete(handles.hmepmax{handles.id_angle, handles.pos})
        delete(handles.hmepmin{handles.id_angle, handles.pos})
    end
end
if isfield(handles,'hmepmax_init') && isfield(handles,'hmepmin_init')
    if ishandle(handles.hmepmax_init) && ishandle(handles.hmepmin_init)
        delete(handles.hmepmax_init)
        delete(handles.hmepmin_init)
    end
end

handles.data.meps_mono{handles.id_angle}{handles.pos} = handles.data.meps_mono_bkp{handles.id_angle}{handles.pos};
handles.data.mepmean{handles.id_angle}{handles.pos} = handles.data.mepmean_bkp{handles.id_angle}{handles.pos};
handles.data.mepmax{handles.id_angle}{handles.pos} = handles.data.mepmax_bkp{handles.id_angle}{handles.pos};
handles.data.mepmin{handles.id_angle}{handles.pos} = handles.data.mepmin_bkp{handles.id_angle}{handles.pos};

handles.data.fmed{handles.id_angle}{handles.pos} = handles.data.fmed_bkp{handles.id_angle}{handles.pos};
handles.data.amp_rms{handles.id_angle}{handles.pos} = handles.data.amp_rms_bkp{handles.id_angle}{handles.pos};

handles.data.offset_rms{handles.id_angle}{handles.pos} = handles.data.offset_rms_bkp{handles.id_angle}{handles.pos};
handles.data.offset_fmed{handles.id_angle}{handles.pos} = handles.data.offset_fmed_bkp{handles.id_angle}{handles.pos};

handles.data.pos_max{handles.id_angle}{handles.pos} = handles.data.pos_max_init{handles.id_angle}{handles.pos};
handles.data.pos_min{handles.id_angle}{handles.pos} = handles.data.pos_min_init{handles.id_angle}{handles.pos};

x_max = handles.data.pos_max{handles.id_angle}{handles.pos}/handles.signal.fsample(handles.id_angle);
x_min = handles.data.pos_min{handles.id_angle}{handles.pos}/handles.signal.fsample(handles.id_angle);

handles.minmax = [x_min handles.data.mepmin{handles.id_angle}{handles.pos}; x_max handles.data.mepmax{handles.id_angle}{handles.pos}];

handles.hdata = plot(handles.xs+handles.MEPStart,handles.data.mepmean{handles.id_angle}{handles.pos});
handles.hmepmax{handles.id_angle, handles.pos} = plot(x_max+handles.MEPStart,...
    handles.data.mepmax{handles.id_angle}{handles.pos},'+g');
handles.hmepmin{handles.id_angle, handles.pos} = plot(x_min+handles.MEPStart,...
    handles.data.mepmin{handles.id_angle}{handles.pos},'+g');

if ~isfield(handles,'minmax');
    handles.output_minmax = [];
else
    handles.output_minmax = handles.minmax;
end
handles.output_latency = [];
handles.output_data = handles.data;
handles.output_signal = handles.signal;

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in button_ManualSelect.
function button_ManualSelect_Callback(hObject, eventdata, handles)
% hObject    handle to button_ManualSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'minmax')
    if ishandle(handles.minmax)
        delete(handles.minmax)
    end
end
% trick to shift the x coordinate using MEPStart
% necessery to combine the plots on this MEP window and the Main one
[x y] = getpts(handles.axes1);
minmax = [x y];

handles.data.ch_excluded{handles.id_angle}(handles.data.ch_excluded{handles.id_angle} == handles.pos)=[];

if ~isempty(minmax)
    axes(handles.axes1);
    handles.hminmax = plot(minmax(:,1),minmax(:,2),'ro');
end
handles.minmax = [x-handles.MEPStart y];

% Update min and max selected
handles.data.pos_min{handles.id_angle}{handles.pos} = round(handles.minmax(1,1)*handles.signal.fsample(handles.id_angle));
handles.data.pos_max{handles.id_angle}{handles.pos} = round(handles.minmax(2,1)*handles.signal.fsample(handles.id_angle));
handles.data.mepmin{handles.id_angle}{handles.pos} = handles.minmax(1,2);
handles.data.mepmax{handles.id_angle}{handles.pos} = handles.minmax(2,2);

% Update handles structure
guidata(hObject, handles);




% --- Executes on button press in button_OK.
function button_OK_Callback(hObject, eventdata, handles)
% hObject    handle to button_OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles,'minmax');
    handles.output_minmax = [];
else
    handles.output_minmax = handles.minmax;
end

if ~isfield(handles,'latency');
    handles.output_latency = [];
else
    handles.output_latency = handles.latency;
end

handles = rmfield(handles, {'xs', 'pos'});

handles.output_data = handles.data;
handles.output_signal = handles.signal;

uiresume

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

handles.minmax = [0 0; 0 0];

handles.data.outliers_mep{handles.id_angle}(handles.data.outliers_mep{handles.id_angle} == handles.pos)=[];
if ~isempty(handles.data.ch_excluded{handles.id_angle})
    if isempty(handles.data.ch_excluded{handles.id_angle}(handles.data.ch_excluded{handles.id_angle} == handles.pos))
        handles.data.ch_excluded{handles.id_angle} = cat(1, handles.data.ch_excluded{handles.id_angle}, handles.pos);
    end
else
    handles.data.ch_excluded{handles.id_angle} = handles.pos;
end

% Update min and max selected
% handles.data.mepmax{handles.id_angle}{handles.pos} = handles.minmax(1,2);
% handles.data.mepmin{handles.id_angle}{handles.pos} = handles.minmax(2,2);
handles.data.pos_max{handles.id_angle}{handles.pos} = [];
handles.data.pos_min{handles.id_angle}{handles.pos} = [];
handles.data.mepmax{handles.id_angle}{handles.pos} = [];
handles.data.mepmin{handles.id_angle}{handles.pos} = [];
handles.latency = [];

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in button_LatencySelect.
function button_LatencySelect_Callback(hObject, eventdata, handles)
% hObject    handle to button_LatencySelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[xlat ylat] = getpts(handles.axes1);
handles.latency = [xlat ylat];


if ~isempty(handles.latency)
    axes(handles.axes1);
    handles.hlatencystart = plot(handles.latency(1,1),handles.latency(1,2),'gv');
end

if length(handles.latency) > 1
    axes(handles.axes1);
    handles.hlatencystop = plot(handles.latency(2,1),handles.latency(2,2),'r^');
end

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in button_mep_outlier.
function button_mep_outlier_Callback(hObject, eventdata, handles)
% hObject    handle to button_mep_outlier (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% count the total number of interpolated electrodes
handles.data.count(handles.id_angle) = handles.data.count(handles.id_angle) + 1;

% coord computation
x = mod(handles.pos, 13);
if x == 0;
    x = 13;
end
y = ceil(handles.pos/13);

handles.data.ch_excluded{handles.id_angle}(handles.data.ch_excluded{handles.id_angle} == handles.pos)=[];
if ~isempty(handles.data.outliers_mep{handles.id_angle})
    if isempty(handles.data.outliers_mep{handles.id_angle}(handles.data.outliers_mep{handles.id_angle} == handles.pos))
        handles.data.outliers_mep{handles.id_angle} = cat(1, handles.data.outliers_mep{handles.id_angle}, handles.pos);
    end
else
    handles.data.outliers_mep{handles.id_angle} = handles.pos;
end

% Outlier computation
emg_diff_aux = reshape(interpolate_chs_nan(reshape(handles.signal.emg_mono{handles.id_angle}',...
    13,5,size(handles.signal.emg_mono{handles.id_angle},1)),[x y]),65,size(handles.signal.emg_mono{handles.id_angle},1))';

% rewritting the emg_map and plotting its new handles.signal
handles.signal.emg_mono{handles.id_angle} = emg_diff_aux;
% -----------------------

% MEPs and minmax computation
if sum(handles.data.trigger{handles.id_angle}(:,handles.pos)) ~= 0
    aux_trigger = handles.data.trigger{handles.id_angle}(:,handles.pos);
else
    aux_trigger = [];
end

aux_channel = handles.signal.emg_mono{handles.id_angle}(:,handles.pos);
aux_meps = [];
aux_meps_column = [];
aux_offset_signal_columns = [];
aux_offset_signal = [];

for k = 1:length(aux_trigger) % Number of stimulus
    handles.data.meps_mono{handles.id_angle}{handles.pos}{k} = aux_channel(aux_trigger(k)+handles.data.s0(handles.id_angle):aux_trigger(k)+handles.data.s1(handles.id_angle)-1);
    handles.data.offset_signals{handles.id_angle}{handles.pos}{k} = aux_channel(aux_trigger(k)-handles.data.s1(handles.id_angle)-30+1:aux_trigger(k)-handles.data.s0(handles.id_angle)-30);
    if k == 1
        aux_meps_column = handles.data.meps_mono{handles.id_angle}{handles.pos}{k};
        aux_meps = handles.data.meps_mono{handles.id_angle}{handles.pos}{k};
        aux_offset_signal_columns = handles.data.offset_signals{handles.id_angle}{k};
        aux_offset_signal = handles.data.offset_signals{handles.id_angle}{handles.pos}{k};
    else
        aux_meps_column = cat(1, aux_meps_column, handles.data.meps_mono{handles.id_angle}{handles.pos}{k});
        aux_meps = cat(2, aux_meps, handles.data.meps_mono{handles.id_angle}{handles.pos}{k});
        aux_offset_signal_columns = cat(1, aux_offset_signal_columns,...
            handles.data.offset_signals{handles.id_angle}{handles.pos}{k});
        aux_offset_signal = cat(2, aux_offset_signal,...
            handles.data.offset_signals{handles.id_angle}{handles.pos}{k});
    end
end

if isempty(aux_trigger)
    handles.data.meps_mono{handles.id_angle}{handles.pos} = [];
    handles.data.offset_signals{handles.id_angle}{handles.pos} = [];
end

% median frequency, rms values and mean frequency for meps
[handles.data.fmed{handles.id_angle}{handles.pos}, handles.data.amp_rms{handles.id_angle}{handles.pos}, ~] = Fmed3cla(aux_meps_column,...
    handles.signal.fsample(handles.id_angle), handles.data.s1(handles.id_angle)-handles.data.s0(handles.id_angle));
handles.data.fmed{handles.id_angle}{handles.pos} = mean(handles.data.fmed{handles.id_angle}{handles.pos}, 1);
handles.data.amp_rms{handles.id_angle}{handles.pos} = mean(handles.data.amp_rms{handles.id_angle}{handles.pos}, 1);

% median frequency and rms value for the offset signal
[handles.data.offset_fmed{handles.id_angle}{handles.pos}, handles.data.offset_rms{handles.id_angle}{handles.pos}, ~] = Fmed3cla(aux_offset_signal_columns,...
    handles.signal.fsample(handles.id_angle), handles.data.s1(handles.id_angle)-handles.data.s0(handles.id_angle));
handles.data.offset_rms{handles.id_angle}{handles.pos} = mean(handles.data.offset_rms{handles.id_angle}{handles.pos}, 1);
handles.data.offset_fmed{handles.id_angle}{handles.pos} = mean(handles.data.offset_fmed{handles.id_angle}{handles.pos}, 1);

handles.data.mepmean{handles.id_angle}{handles.pos} = mean(aux_meps,2);
handles.data.mepmean_outlined = handles.data.mepmean;
handles.data.mepmax{handles.id_angle}{handles.pos} = max(handles.data.mepmean{handles.id_angle}{handles.pos});
handles.data.mepmin{handles.id_angle}{handles.pos} = min(handles.data.mepmean{handles.id_angle}{handles.pos});

handles.data.pos_max{handles.id_angle}{handles.pos} = find(handles.data.mepmean{handles.id_angle}{handles.pos}...
    == handles.data.mepmax{handles.id_angle}{handles.pos});
handles.data.pos_min{handles.id_angle}{handles.pos} = find(handles.data.mepmean{handles.id_angle}{handles.pos}...
    == handles.data.mepmin{handles.id_angle}{handles.pos});

x_max = handles.data.pos_max{handles.id_angle}{handles.pos}/handles.signal.fsample(handles.id_angle);
x_min = handles.data.pos_min{handles.id_angle}{handles.pos}/handles.signal.fsample(handles.id_angle);

if isfield(handles,'minmax')
    if ishandle(handles.minmax)
        delete(handles.minmax)
    end
end
handles.minmax = [x_min handles.data.mepmin{handles.id_angle}{handles.pos}; x_max handles.data.mepmax{handles.id_angle}{handles.pos}];
% ------------------------

% output atribution
if isfield(handles,'minmax') == 0;
    handles.output_minmax = [];
else
    handles.output_minmax = handles.minmax;
end

if isfield(handles,'latency') == 0;
    handles.output_latency = [];
else
    handles.output_latency = handles.latency;
end

handles.output_data = handles;
% --------------

% plot of data and minmax outlined
if isfield(handles,'hdata')
    if ishandle(handles.hdata)
        delete(handles.hdata)
    end
end
if isfield(handles,'hmepmax_init') && isfield(handles,'hmepmin_init')
    if ishandle(handles.hmepmax_init) && ishandle(handles.hmepmin_init)
        delete(handles.hmepmax_init)
        delete(handles.hmepmin_init)
    end
end
handles.hdata = plot(handles.xs+handles.MEPStart,handles.data.mepmean{handles.id_angle}{handles.pos});
handles.hmepmax{handles.id_angle, handles.pos} = plot(x_max+handles.MEPStart,...
    handles.data.mepmax{handles.id_angle}{handles.pos},'+r');
handles.hmepmin{handles.id_angle, handles.pos} = plot(x_min+handles.MEPStart,...
    handles.data.mepmin{handles.id_angle}{handles.pos},'+r');
% -----------------------

% Update handles structure
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
