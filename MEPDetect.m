function varargout = MEPDetect(varargin)
% MEPDETECT MATLAB code for MEPDetect.fig
%      MEPDETECT, by itself, creates a new MEPDETECT or raises the existing
%      singleton*.
%
%      H = MEPDETECT returns the handle to a new MEPDETECT or the handle to
%      the existing singleton*.
%
%      MEPDETECT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MEPDETECT.M with the given input arguments.
%
%      MEPDETECT('Property','Value',...) creates a new MEPDETECT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MEPDetect_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MEPDetect_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MEPDetect

% Last Modified by GUIDE v2.5 22-May-2013 17:34:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MEPDetect_OpeningFcn, ...
                   'gui_OutputFcn',  @MEPDetect_OutputFcn, ...
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


% --- Executes just before MEPDetect is made visible.
function MEPDetect_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MEPDetect (see VARARGIN)


if ~isempty (varargin)
    handles.data = varargin{1};
    handles.MEPStart = varargin{2};
    handles.xs = varargin{3};
    handles.pos = varargin{4};
        
    axes(handles.axes1);
    hold on   
end

% get current values of min max positions and plot them
x_max_init =  handles.data.pos_max{handles.data.channel, handles.pos}/handles.data.configuration.fsample;
x_min_init =  handles.data.pos_min{handles.data.channel, handles.pos}/handles.data.configuration.fsample;

handles.hdata = plot(handles.xs+handles.MEPStart,handles.data.mepmean{handles.data.channel}{handles.pos});
if ~isempty(handles.data.mepmax{handles.data.channel}{handles.pos}) && ~isempty(handles.data.mepmin{handles.data.channel}{handles.pos})
    if handles.data.mepmax{handles.data.channel}{handles.pos} ~= 0 && handles.data.mepmin{handles.data.channel}{handles.pos} ~= 0
        handles.hmepmax_init = plot(x_max_init + handles.MEPStart, handles.data.mepmax{handles.data.channel}{handles.pos},'+g');
        handles.hmepmin_init = plot(x_min_init + handles.MEPStart, handles.data.mepmin{handles.data.channel}{handles.pos},'+g');
    end
end
% -------------------------------------

% Choose default command line output for MEPDetect
handles.minmax = [x_min_init handles.data.mepmin{handles.data.channel}{handles.pos};...
    x_max_init handles.data.mepmax{handles.data.channel}{handles.pos}];
handles.output_minmax = handles.minmax;
handles.output_data = handles.data;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MEPDetect wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MEPDetect_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiwait(hObject);
handles = guidata(hObject);

delete(hObject);

% Get default command line output from handles structure
varargout{1} = handles.output_minmax;
varargout{2} = handles.output_lattency;
varargout{3} = handles.output_data;


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
if isfield(handles,'hmepmax') && isfield(handles,'hmepmin')
    if ishandle(handles.hmepmax{handles.data.channel, handles.pos}) && ishandle(handles.hmepmin{handles.data.channel, handles.pos})
        delete(handles.hmepmax{handles.data.channel, handles.pos})
        delete(handles.hmepmin{handles.data.channel, handles.pos})
    end
end
if isfield(handles,'hmepmax_init') && isfield(handles,'hmepmin_init')
    if ishandle(handles.hmepmax_init) && ishandle(handles.hmepmin_init)
        delete(handles.hmepmax_init)
        delete(handles.hmepmin_init)
    end
end

handles.data.mepmean{handles.data.channel}{handles.pos} = handles.data.mepmean_bkp{handles.data.channel}{handles.pos};
handles.data.mepmax{handles.data.channel}{handles.pos} = handles.data.mepmax_bkp{handles.data.channel}{handles.pos};
handles.data.mepmin{handles.data.channel}{handles.pos} = handles.data.mepmin_bkp{handles.data.channel}{handles.pos};

handles.data.fmed{handles.data.channel}{handles.pos} = handles.data.fmed_bkp{handles.data.channel}{handles.pos};
handles.data.amp_rms{handles.data.channel}{handles.pos} = handles.data.amp_rms_bkp{handles.data.channel}{handles.pos};
handles.data.fmean{handles.data.channel}{handles.pos} = handles.data.fmean_bkp{handles.data.channel}{handles.pos};

handles.data.offset_rms{handles.data.channel}{handles.pos} = handles.data.offset_rms_bkp{handles.data.channel}{handles.pos};
handles.data.offset_fmed{handles.data.channel}{handles.pos} = handles.data.offset_fmed_bkp{handles.data.channel}{handles.pos};

handles.data.pos_max{handles.data.channel, handles.pos} = handles.data.pos_max_init{handles.data.channel, handles.pos};
handles.data.pos_min{handles.data.channel, handles.pos} = handles.data.pos_min_init{handles.data.channel, handles.pos};

x_max = handles.data.pos_max{handles.data.channel, handles.pos}/handles.data.configuration.fsample;
x_min = handles.data.pos_min{handles.data.channel, handles.pos}/handles.data.configuration.fsample;

handles.minmax = [x_min handles.data.mepmin{handles.data.channel}{handles.pos}; x_max handles.data.mepmax{handles.data.channel}{handles.pos}];

handles.hdata = plot(handles.xs+handles.MEPStart,handles.data.mepmean{handles.data.channel}{handles.pos});
handles.hmepmax{handles.data.channel, handles.pos} = plot(x_max+handles.MEPStart,...
    handles.data.mepmax{handles.data.channel}{handles.pos},'+g');
handles.hmepmin{handles.data.channel, handles.pos} = plot(x_min+handles.MEPStart,...
    handles.data.mepmin{handles.data.channel}{handles.pos},'+g');

if isfield(handles,'minmax') == 0;
    handles.output_minmax = [];
else
    handles.output_minmax = handles.minmax;
end
handles.output_lattency = [];
handles.output_data = handles.data;

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

if ~isempty(minmax)
    axes(handles.axes1);
    handles.hminmax = plot(minmax(:,1),minmax(:,2),'ro');
end
handles.minmax = [x-handles.MEPStart y];

% Update min and max selected
handles.data.pos_min{handles.data.channel, handles.pos} = round(handles.minmax(1,1)*handles.data.configuration.fsample);
handles.data.pos_max{handles.data.channel, handles.pos} = round(handles.minmax(2,1)*handles.data.configuration.fsample);
handles.data.mepmin{handles.data.channel}{handles.pos} = handles.minmax(1,2);
handles.data.mepmax{handles.data.channel}{handles.pos} = handles.minmax(2,2);

% Update handles structure
guidata(hObject, handles);




% --- Executes on button press in button_OK.
function button_OK_Callback(hObject, eventdata, handles)
% hObject    handle to button_OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'minmax') == 0;
    handles.output_minmax = [];
else
    handles.output_minmax = handles.minmax;
end

if isfield(handles,'latency') == 0;
    handles.output_lattency = [];
else
    handles.output_lattency = handles.latency;
end

handles.output_data = handles.data;

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

% Update min and max selected
% handles.data.mepmax{handles.data.channel}{handles.pos} = handles.minmax(1,2);
% handles.data.mepmin{handles.data.channel}{handles.pos} = handles.minmax(2,2);
handles.data.pos_max{handles.data.channel, handles.pos} = [];
handles.data.pos_min{handles.data.channel, handles.pos} = [];
handles.data.mepmax{handles.data.channel}{handles.pos} = [];
handles.data.mepmin{handles.data.channel}{handles.pos} = [];
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

% coord computation
x = mod(handles.pos, 13);
if x == 0;
    x = 13;
end
y = ceil(handles.pos/13);
% -----------------

% Outlier computation
% if isfield(handles.data.signal,'emg_data_mono') && handles.data.channel == 1
%     handles.data.signal.emg_data_mono = reshape(interpolate_chs_nan(reshape(handles.data.signal.emg_data_mono',...
%         13,5,size(handles.data.signal.emg_data_mono,1)),[x y]),65,size(handles.data.signal.emg_data_mono,1))';
%     for i = 1:size(handles.data.signal.emg_data_mono, 2)
%         handles.data.signal.emg_map{i}(:, handles.data.channel) = handles.data.signal.emg_data_mono(:, i);
%     end
% else
%     handles.data.signal.emg_data_diff = reshape(interpolate_chs_nan(reshape(handles.data.signal.emg_data_diff',...
%         13,5,size(handles.data.signal.emg_data_diff,1)),[x y]),65,size(handles.data.signal.emg_data_diff,1))';
%     for i = 1:size(handles.data.signal.emg_data_diff, 2)
%         handles.data.signal.emg_map{i}(:, handles.data.channel) = handles.data.signal.emg_data_diff(:, i);
%     end
% end

% This will only work for hdsemg_diff_bin - because of the channel
% selection
% TODO: Fix it for a general case
if size(handles.data.signal.emg_map{1}, 2) == 1
    handles.data.signal.emg_data_diff = reshape(interpolate_chs_nan(reshape(handles.data.signal.emg_data_diff',...
        13,5,size(handles.data.signal.emg_data_diff,1)),[x y]),65,size(handles.data.signal.emg_data_diff,1))';
    for i = 1:size(handles.data.signal.emg_data_diff, 2)
        handles.data.signal.emg_map{i}(:, handles.data.channel) = handles.data.signal.emg_data_diff(:, i);
    end
else
    handles.data.signal.emg_data_mono = reshape(interpolate_chs_nan(reshape(handles.data.signal.emg_data_mono',...
        13,5,size(handles.data.signal.emg_data_mono,1)),[x y]),65,size(handles.data.signal.emg_data_mono,1))';
    for i = 1:size(handles.data.signal.emg_data_mono, 2)
        handles.data.signal.emg_map{i}(:, handles.data.channel) = handles.data.signal.emg_data_mono(:, i);
    end
end
% -----------------------

% MEPs and minmax computation
off0 = -1500;
off1 = -10;
aux_channel = handles.data.signal.emg_map{handles.pos}(:,handles.data.channel);
aux_meps = zeros(handles.data.s1-handles.data.s0,1);
aux_meps_column = [];
aux_offset_signal_columns = [];
aux_offset_signal = [];

if ~isempty(handles.data.trigger{handles.data.channel, handles.pos})
    aux_trigger = handles.data.trigger{handles.data.channel, handles.pos}(:,1);
else
    aux_trigger = [];
end

for k = 1:length(aux_trigger)
    % offset adjustment
%     aux_offset = aux_channel(aux_trigger(k)+off0:aux_trigger(k)+off1-1);
%     offset{handles.data.channel}{handles.pos}(k) = mean(aux_offset);
    % ----------------
    handles.data.meps{handles.data.channel}{handles.pos}{k} = aux_channel(aux_trigger(k)+...
        handles.data.s0:aux_trigger(k)+handles.data.s1-1);
    handles.data.offset_signals{handles.data.channel}{handles.pos}{k} = aux_channel(aux_trigger(k)-handles.data.s1-30+1:aux_trigger(k)-handles.data.s0-30);
    
    if k == 1
        aux_meps = [aux_meps,handles.data.meps{handles.data.channel}{handles.pos}{k}];
        aux_meps = handles.data.meps{handles.data.channel}{handles.pos}{k};
        aux_offset_signal_columns = handles.data.offset_signals{handles.data.channel}{handles.pos}{k};
        aux_offset_signal = handles.data.offset_signals{handles.data.channel}{handles.pos}{k};
    else
        aux_meps_column = cat(1, aux_meps_column, handles.data.meps{handles.data.channel}{handles.pos}{k});
        aux_meps = cat(2, aux_meps, handles.data.meps{handles.data.channel}{handles.pos}{k});
        aux_offset_signal_columns = cat(1, aux_offset_signal_columns, handles.data.offset_signals{handles.data.channel}{handles.pos}{k});
        aux_offset_signal = cat(2, aux_offset_signal, handles.data.offset_signals{handles.data.channel}{handles.pos}{k});
    end
end

if isempty(aux_trigger)
    handles.data.meps{handles.data.channel}{handles.pos} = [];
end

% median frequency, rms values and mean frequency for meps and offset
[handles.data.fmed{handles.data.channel}{handles.pos}, handles.data.amp_rms{handles.data.channel}{handles.pos}, handles.data.fmean{handles.data.channel}{handles.pos}] = Fmed3cla(aux_meps_column,...
    handles.data.configuration.fsample, handles.data.s1-handles.data.s0);
handles.data.fmed{handles.data.channel}{handles.pos} = mean(handles.data.fmed{handles.data.channel}{handles.pos}, 1);
handles.data.amp_rms{handles.data.channel}{handles.pos} = mean(handles.data.amp_rms{handles.data.channel}{handles.pos}, 1);
handles.data.fmean{handles.data.channel}{handles.pos} = mean(handles.data.fmean{handles.data.channel}{handles.pos}, 1);

[handles.data.offset_fmed{handles.data.channel}{handles.pos}, handles.data.offset_rms{handles.data.channel}{handles.pos}, ~] = Fmed3cla(aux_offset_signal_columns,...
    handles.data.configuration.fsample, handles.data.s1-handles.data.s0);
handles.data.offset_fmed{handles.data.channel}{handles.pos} = mean(handles.data.offset_fmed{handles.data.channel}{handles.pos}, 1);
handles.data.offset_rms{handles.data.channel}{handles.pos} = mean(handles.data.offset_rms{handles.data.channel}{handles.pos}, 1);

handles.data.mepmean{handles.data.channel}{handles.pos} = mean(aux_meps,2);
handles.data.mepmean_outlined = handles.data.mepmean;
handles.data.mepmax{handles.data.channel}{handles.pos} = max(handles.data.mepmean{handles.data.channel}{handles.pos});
handles.data.mepmin{handles.data.channel}{handles.pos} = min(handles.data.mepmean{handles.data.channel}{handles.pos});

handles.data.pos_max{handles.data.channel, handles.pos} = find(handles.data.mepmean{handles.data.channel}{handles.pos}...
    == handles.data.mepmax{handles.data.channel}{handles.pos});
handles.data.pos_min{handles.data.channel, handles.pos} = find(handles.data.mepmean{handles.data.channel}{handles.pos}...
    == handles.data.mepmin{handles.data.channel}{handles.pos});

x_max = handles.data.pos_max{handles.data.channel, handles.pos}/handles.data.configuration.fsample;
x_min = handles.data.pos_min{handles.data.channel, handles.pos}/handles.data.configuration.fsample;

if isfield(handles,'minmax')
    if ishandle(handles.minmax)
        delete(handles.minmax)
    end
end
handles.minmax = [x_min handles.data.mepmin{handles.data.channel}{handles.pos}; x_max handles.data.mepmax{handles.data.channel}{handles.pos}];
% ------------------------

% output atribution
if isfield(handles,'minmax') == 0;
    handles.output_minmax = [];
else
    handles.output_minmax = handles.minmax;
end

if isfield(handles,'latency') == 0;
    handles.output_lattency = [];
else
    handles.output_lattency = handles.latency;
end

handles.output_data = handles.data;
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
handles.hdata = plot(handles.xs+handles.MEPStart,handles.data.mepmean{handles.data.channel}{handles.pos});
handles.hmepmax{handles.data.channel, handles.pos} = plot(x_max+handles.MEPStart,...
    handles.data.mepmax{handles.data.channel}{handles.pos},'+r');
handles.hmepmin{handles.data.channel, handles.pos} = plot(x_min+handles.MEPStart,...
    handles.data.mepmin{handles.data.channel}{handles.pos},'+r');
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

function [fmed, rms, fmean] = Fmed3cla(x,fsamp,epoch)
%[P,f]=psd(x,fsamp*epoch_len,fsamp,ones(fsamp*epoch_len,1),0);

if mod(length(x),epoch), error('vector length should be multiple of the epoch chosen'), return, end
if all(isnan(x)),
    fmed = nan(round(size(x,1)/epoch),1);
    rms = fmed;
    fmean = rms;
    return
end
    
ReshapedX = detrend(reshape(x,epoch,numel(x)/epoch),'linear');
rms = zeros(size(ReshapedX,2), 1);
fmean = zeros(size(ReshapedX,2), 1);
fmed = zeros(size(ReshapedX,2), 1);
for column = 1:size(ReshapedX,2)
    x = ReshapedX(:,column);
    
    rms(column,1)=norm(x)/sqrt(length(x));
    x=x-mean(x);
%     [P,f] = psd(x,fsamp,fsamp,boxcar(length(x)),0);
    % I changed psd to pwelch to supress matlab warning
    % I made tests and pwelch and psd return very similar P and f outputs
    % Doing hamming windowing besides boxcar the difference increases
    [P, f] = pwelch(x, boxcar(length(x)), 0, fsamp, fsamp);

    if sum(P)~=0
        num=sum(f.*P);
        den=sum(P);
        fmean(column,1)=num/den;
        den=sum(P);
        k=1;
        while (sum(P(1:k)))<=den/2
            k=k+1;
        end
        i=k;

        if i<length(P)
            alfa1=sum(P(1:i-1))/den;
            if P(i)<=P(i+1),
                radi=roots([abs(P(i+1)-P(i)) 2*P(i) -2*(0.5-alfa1)*den]);
            else
                radi=roots([abs(P(i+1)-P(i)) P(i)+P(i+1) -2*(0.5-alfa1)*den]);
            end
            if radi(1)>0 && radi(1)<1,
                x=radi(1);
            else
                x=radi(2);
            end
            df=f(2)-f(1);
            fmed(column,1)=f(k)+df*x;
        else
            fmean(column,1)=NaN;
            fmed(column,1)=NaN;
        end
    else
        fmean(column,1)=NaN;
        fmed(column,1)=NaN;
    end
end
