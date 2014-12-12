function varargout = MapHDsEMG(varargin)
% MAPHDSEMG MATLAB code for MapHDsEMG.fig
%      MAPHDSEMG, by itself, creates a new MAPHDSEMG or raises the existing
%      singleton*.
%
%      H = MAPHDSEMG returns the handle to a new MAPHDSEMG or the handle to
%      the existing singleton*.
%
%      MAPHDSEMG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAPHDSEMG.M with the given input arguments.
%
%      MAPHDSEMG('Property','Value',...) creates a new MAPHDSEMG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MapHDsEMG_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MapHDsEMG_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MapHDsEMG

% Last Modified by GUIDE v2.5 15-Aug-2013 12:24:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MapHDsEMG_OpeningFcn, ...
                   'gui_OutputFcn',  @MapHDsEMG_OutputFcn, ...
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


% --- Executes just before MapHDsEMG is made visible.
function MapHDsEMG_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MapHDsEMG (see VARARGIN)


if ~isempty (varargin)
    handles.data = varargin{1};
end

handles.signalvalue = 0;
handles.mepmeanvalue = 0;
handles.mepsvalue = 0;
handles.linevalue = 0;
handles.triggervalue = 0;
handles.mepminmaxvalue = 0;
handles.ampthresholdvalue = 0;
handles.latencyvalue = 0;

% Data to export
handles.data.amp_pp_sim_diff = 0;
handles.data.amp_rms_sim_diff = 0;
handles.data.latency_sim = 0;
handles.data.fmed_sim_diff = 0;
handles.data.cog_cluster_diff_pp = [0 0];
handles.data.size_cluster_diff_pp = 0;
handles.data.cluster_amp_pp_diff = 0;
handles.data.freq_cluster_pp = 0;
handles.data.cog_cluster_diff_rms = [0 0];
handles.data.size_cluster_diff_rms = 0;
handles.data.cluster_amp_rms_diff = 0;
handles.data.freq_cluster_rms = 0;
handles.data.cvelocity = 0;
handles.data.correl_coef = 0;
handles.data.iz_row = 0;
handles.data.mean_total_rms = 0;
handles.data.mean_offset_rms = 0;
handles.data.latency_avarege = 0;

handles.data.outliers_mep = {};
handles.data.outliers_trig = {};

handles.data.channel = 1;

for i = 1:size(handles.data.signal.emg_map{handles.data.channel},2)
    handles.data.outliers_trig{i} = [];
    handles.data.outliers_mep{i} = [];
end

channels_names = '';

for i = 1:length(handles.data.configuration.emg_mode)
    channels_names = strcat(channels_names,num2str(i),'-',handles.data.configuration.emg_mode(i),'__');
end

% figure name
name = get(handles.figure1, 'Name');
fig_name = [name ' ' handles.data.configuration.patient_id '_' handles.data.configuration.emg_side];
set(handles.figure1,'Name',fig_name)

set(handles.text_ChannelLegend, 'String',channels_names)
set(handles.text_Channel, 'String',num2str(handles.data.channel))

set(handles.edit_id_number, 'String',num2str(handles.data.configuration.patient_id))
set(handles.edit_hemisphere_side, 'String',num2str(handles.data.configuration.emg_side))
set(handles.edit_number_angle_stim, 'String',num2str(handles.data.configuration.angle_stim))

% Plot signal on each axes
for i = 1:length(handles.data.signal.emg_map)
    axes(eval(strcat('handles.axes',num2str(i))))
    hold on
    for j = 1:size(handles.data.signal.emg_map{handles.data.channel},2)
        handles.data.xs{j}{i} = (1:length(handles.data.signal.emg_map{i}(:,j)))/handles.data.configuration.fsample;
        handles.hsignal{j}(i) = plot(handles.data.xs{j}{i},handles.data.signal.emg_map{i}(:,j),'Visible','off');
    end
    set(eval(strcat('handles.axes',num2str(i))),'HitTest','on')    
    update_waitbar(handles,i/length(handles.data.signal.emg_map))
end

set(handles.hsignal{handles.data.channel},'Visible','on')
set(handles.checkbox_ViewSignal,'Value',1)


% Generate latency cell and initial total amplitude
for i = 1:size(handles.data.signal.emg_map{handles.data.channel},2)
    for j = 1:length(handles.data.signal.emg_map)
        handles.data.latency{i}{j} = [];
        handles.data.duration{i}{j} = [];
        handles.hlatencystart{i}{j} = [];
        handles.hlatencystop{i}{j} = [];
    end
end

% Choose default command line output for MapHDsEMG
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MapHDsEMG wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MapHDsEMG_OutputFcn(hObject, eventdata, handles) 
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

set(handles.checkbox_ViewSignal,'Value',0)
checkbox_ViewSignal_Callback(hObject, eventdata, handles)

for j = 1:size(handles.data.signal.emg_map{handles.data.channel},2)
    if isfield(handles,'hsignal')
        if ishandle(handles.hsignal{j})
            delete(handles.hsignal{j})
        end
    end
    for i = 1:length(handles.data.signal.emg_map)
        if isfield(handles,'htrigger')
            if ishandle(handles.htrigger{j, i})
                delete(handles.htrigger{j, i})
            end
        end
        
        if isfield(handles,'hline') == 1
            if ishandle(handles.hline{j, i})
                delete(handles.hline{j, i})
            end
        end
    end
end

clear handles.data.slope
handles.data.slope = get(handles.edit_Threshold,'String');
handles.data.slope = str2double(handles.data.slope);

hp_filter = str2double(get(handles.edit_HPF,'String'));
lp_filter = str2double(get(handles.edit_LPF,'String'));

[b,a] = butter(2,[hp_filter lp_filter]*2/handles.data.configuration.fsample);

% ad_range - of the emg equipment
% conv_uv - multiplication factor to see signal in microVolts
% TODO: put a text box in the UI for the user to change the conversion
% factor
ad_range = 5.0;
conv_uv = 1000000.0;

handles.data.signal.emg_data = filtfilt(b,a,(handles.data.signal.emg_data*ad_range*conv_uv)/(handles.data.configuration.signal_gain*2^(handles.data.configuration.ad_bits)));

handles.data.signal.emg_data_mono = [nan(length(handles.data.signal.emg_data), 1) handles.data.signal.emg_data(:, 1:11)...
    nan(length(handles.data.signal.emg_data), 1) handles.data.signal.emg_data(:, 12:50) nan(length(handles.data.signal.emg_data), 1)...
    handles.data.signal.emg_data(:, 51:end) nan(length(handles.data.signal.emg_data), 1)];

handles.data.signal.emg_data_diff = diff(handles.data.signal.emg_data,1,2);
handles.data.signal.emg_data_diff = [nan(length(handles.data.signal.emg_data_diff), 1) handles.data.signal.emg_data_diff(:, 1:10)...
    nan(length(handles.data.signal.emg_data_diff), 1) nan(length(handles.data.signal.emg_data_diff), 1) handles.data.signal.emg_data_diff(:, 12:23)...
    nan(length(handles.data.signal.emg_data_diff), 1) handles.data.signal.emg_data_diff(:, 25:36)...
    nan(length(handles.data.signal.emg_data_diff), 1) handles.data.signal.emg_data_diff(:, 38:49)...
    nan(length(handles.data.signal.emg_data_diff), 1) nan(length(handles.data.signal.emg_data_diff), 1)...
    handles.data.signal.emg_data_diff(:, 51:end) nan(length(handles.data.signal.emg_data_diff), 1) nan(length(handles.data.signal.emg_data_diff), 1)];

handles.data.signal.emg_data_all = {handles.data.signal.emg_data_mono handles.data.signal.emg_data_diff};


% trigger computation based no HDsEMG aux channel with trigger signal
if handles.data.slope == 0
    trigger = handles.data.signal.raw_data(:,65)*5/(2^(handles.data.configuration.ad_bits-1));% Triggering EMGs
    
    % finding samples denoting trigger onset
    samples_triggeron = find(trigger>4);
    samples_triggeron = (samples_triggeron(diff([-inf;samples_triggeron])>1));
    
    % rewritting the emg_map and plotting its new signal beacuse of the
    % butterworth applied just before
    % TODO: This trick is temporary just to use with emg_diff map,
    % fix it to the general case.
    for i = 1:length(handles.data.signal.emg_map)
        if size(handles.data.signal.emg_map{1}, 2) == 1
            handles.data.signal.emg_map{i} = handles.data.signal.emg_data_diff(:, i);
        else
            for j = 1:size(handles.data.signal.emg_map{handles.data.channel},2)
                handles.data.signal.emg_map{i}(:, j) = handles.data.signal.emg_data_all{j}(:, i);
            end
        end
    end

      
    for i = 1:length(handles.data.signal.emg_map)
        axes(eval(strcat('handles.axes',num2str(i))))
        for j = 1:size(handles.data.signal.emg_map{handles.data.channel},2)
            
            handles.hsignal{j}(i) = plot(handles.data.xs{j}{i},handles.data.signal.emg_map{i}(:,j),'Visible','off');
            
            if sum(isnan(handles.data.signal.emg_map{i}(samples_triggeron, j)))==0
                handles.data.trigger{j, i} = [samples_triggeron,...
                    handles.data.signal.emg_map{i}(samples_triggeron, j)];
            else
                handles.data.trigger{j, i} = [];
            end
            % red circulus
            if ~isempty(handles.data.trigger{j, i})
                handles.htrigger{j, i} = plot(handles.data.xs{j}{i}(handles.data.trigger{j,i}(:,1)),...
                    handles.data.trigger{j,i}(:,2),'ro', 'Visible', 'off');
            else
                handles.htrigger{j, i} = nan;
            end
            % green line
            if ~isempty(handles.data.trigger{j, i})
                a = line([handles.data.xs{j}{i}(handles.data.trigger{j, i}(:,1))',...
                    handles.data.xs{j}{i}(handles.data.trigger{j, i}(:,1))'],...
                    [min(handles.data.signal.emg_map{i}(:,j)),...
                    max(handles.data.signal.emg_map{i}(:,j))],...
                    'Color','g','Visible','off');
                handles.hline{j, i} = a;
            else
                handles.hline{j, i} = nan;
            end
        end
        
        set(eval(strcat('handles.axes',num2str(i))),'HitTest','on')
        update_waitbar(handles,i/length(handles.data.signal.emg_map))        
    end
% trigger computation based on derivative properties      
else
    for i = 1:length(handles.data.signal.emg_map) 
        axes(eval(strcat('handles.axes',num2str(i))))           
        handles.data.trigger{i} = BiopacTrigger(handles.data.signal.emg_map{i}(:,handles.data.channel),handles.data.slope);
        
        % red circulus  
        if ~isempty(handles.data.trigger{i})
            handles.htrigger(i) = plot(handles.data.xs{handles.data.channel}{i}(handles.data.trigger{i}(:,1)),handles.data.trigger{i}(:,2),'ro', 'Visible', 'off');
        else
            handles.htrigger(i) = nan;
        end
        set(eval(strcat('handles.axes',num2str(i))),'HitTest','on')
        
        % green line
        if ~isempty(handles.data.trigger{i})
            a = line([handles.data.xs{handles.data.channel}{i}(handles.data.trigger{i}(:,1))',handles.data.xs{handles.data.channel}{i}(handles.data.trigger{i}(:,1))'],...
                [min(handles.data.signal.emg_map{i}(:,handles.data.channel)),...
                max(handles.data.signal.emg_map{i}(:,handles.data.channel))],...
                'Color','g','Visible','off');
            handles.hline{i} = a;
        else
            handles.hline{i} = nan;
        end
        update_waitbar(handles,i/length(handles.data.signal.emg_map))
    end
end

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
if handles.data.channel <= 1
    handles.data.channel = size(handles.data.signal.emg_map{1},2);
    set(handles.text_Channel, 'String',num2str(handles.data.channel))
else
    handles.data.channel = handles.data.channel-1;
    set(handles.text_Channel, 'String',num2str(handles.data.channel))
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
if handles.data.channel >= size(handles.data.signal.emg_map{1},2)    
   handles.data.channel = 1; 
   set(handles.text_Channel, 'String',num2str(handles.data.channel))   
else
    handles.data.channel = handles.data.channel+1;        
    set(handles.text_Channel, 'String',num2str(handles.data.channel))      
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


% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in button_MEPs.
function button_MEPs_Callback(hObject, eventdata, handles)
% hObject    handle to button_MEPs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% change checkbox values for MEP visualization
set(handles.checkbox_ViewSignal,'Value',0)
set(handles.checkbox_ViewTrigger,'Value',0)
set(handles.checkbox_ViewLine,'Value',0)

checkbox_ViewSignal_Callback(hObject, eventdata, handles)
checkbox_ViewTrigger_Callback(hObject, eventdata, handles)
checkbox_ViewLine_Callback(hObject, eventdata, handles)

% MEP windowing and calculation
handles.MEPStart = round(str2double(get(handles.edit_MEPStart,'String')))/1000;
handles.data.s0 = round(str2double(get(handles.edit_MEPStart,'String'))*handles.data.configuration.fsample/1000);
handles.data.s1 = round(str2double(get(handles.edit_MEPEnd,'String'))*handles.data.configuration.fsample/1000);

handles.data.fmed = {};
handles.data.amp_rms = {};
handles.data.offset_fmed = {};
handles.data.offset_rms = {};

off0 = -1500;
off1 = -10;

% aux_colors = rand(10,3);
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

for j = 1:length(handles.data.signal.emg_map) % Electrodes
    axes(eval(strcat('handles.axes',num2str(j))))
    hold on
    for i = 1:size(handles.data.signal.emg_map{handles.data.channel},2) % EMG mode
        
        aux_channel = handles.data.signal.emg_map{j}(:,i);   
        aux_meps = [];
        aux_meps_column = [];
        aux_offset_signal_columns = [];
        aux_offset_signal = [];
        colors = aux_colors;
        if ~isempty(handles.data.trigger{i, j})
            aux_trigger = handles.data.trigger{i, j}(:,1);
        else
            aux_trigger = [];
            colors = [];
        end
                
        for k = 1:length(aux_trigger) % Number of stimulus
            % this offset window the signal before the trigger point to
            %calculus the mean value and then subtract from MEP values
%             aux_offset = aux_channel(aux_trigger(k)+off0:aux_trigger(k)+off1-1);
%             offset{i}{j}(k) = mean(aux_offset);
            
            % offset adjustment - didnt understand why is necessary to
            % subtract the offset
%             handles.data.meps{i}{j}{k} = aux_channel(aux_trigger(k)+handles.data.s0:aux_trigger(k)+handles.data.s1-1) - offset{i}{j}(k);
            handles.data.meps{i}{j}{k} = aux_channel(aux_trigger(k)+handles.data.s0:aux_trigger(k)+handles.data.s1-1);
            handles.data.offset_signals{i}{j}{k} = aux_channel(aux_trigger(k)-handles.data.s1-30+1:aux_trigger(k)-handles.data.s0-30);
            handles.hmeps{i, j}(k) = plot(handles.data.xs{i}{j}(1:handles.data.s1-handles.data.s0)+handles.MEPStart,...
                handles.data.meps{i}{j}{k},'-','color',colors(k,:),'Visible', 'off');
            if k == 1
                aux_meps_column = handles.data.meps{i}{j}{k};
                aux_meps = handles.data.meps{i}{j}{k};
                aux_offset_signal_columns = handles.data.offset_signals{i}{j}{k};
                aux_offset_signal = handles.data.offset_signals{i}{j}{k};
            else
                aux_meps_column = cat(1, aux_meps_column, handles.data.meps{i}{j}{k});
                aux_meps = cat(2, aux_meps, handles.data.meps{i}{j}{k});
                aux_offset_signal_columns = cat(1, aux_offset_signal_columns, handles.data.offset_signals{i}{j}{k});
                aux_offset_signal = cat(2, aux_offset_signal, handles.data.offset_signals{i}{j}{k});
            end
        end
        
        if isempty(aux_trigger)
            handles.data.meps{i}{j} = []; 
        end
        
        % median frequency, rms values and mean frequency for meps
        [handles.data.fmed{i}{j}, handles.data.amp_rms{i}{j}, handles.data.fmean{i}{j}] = Fmed3cla(aux_meps_column,...
            handles.data.configuration.fsample, handles.data.s1-handles.data.s0);
        handles.data.fmed{i}{j} = mean(handles.data.fmed{i}{j}, 1);
        handles.data.amp_rms{i}{j} = mean(handles.data.amp_rms{i}{j}, 1);
        handles.data.fmean{i}{j} = mean(handles.data.fmean{i}{j}, 1);        
        
        % median frequency and rms value for the offset signal
        [handles.data.offset_fmed{i}{j}, handles.data.offset_rms{i}{j}, ~] = Fmed3cla(aux_offset_signal_columns,...
            handles.data.configuration.fsample, handles.data.s1-handles.data.s0);
        handles.data.offset_rms{i}{j} = mean(handles.data.offset_rms{i}{j}, 1);
        handles.data.offset_fmed{i}{j} = mean(handles.data.offset_fmed{i}{j}, 1);
        
        % mean values, max, min and their backups
        handles.data.mepmean{i}{j} = mean(aux_meps,2);
        handles.data.mepmean_bkp = handles.data.mepmean;
        handles.data.mepmax{i}{j} = max(handles.data.mepmean{i}{j});
        handles.data.mepmax_bkp = handles.data.mepmax;
        handles.data.mepmin{i}{j} = min(handles.data.mepmean{i}{j});
        handles.data.mepmin_bkp = handles.data.mepmin;
        handles.data.fmed_bkp = handles.data.fmed;
        handles.data.amp_rms_bkp = handles.data.amp_rms;
        handles.data.fmean_bkp = handles.data.fmean;
        handles.data.offset_rms_bkp = handles.data.offset_rms;
        handles.data.offset_fmed_bkp = handles.data.offset_fmed;
              
        % position in x axis of min and max values
        handles.data.pos_max{i, j} = find(handles.data.mepmean{i}{j} == handles.data.mepmax{i}{j});
        handles.data.pos_min{i, j} = find(handles.data.mepmean{i}{j} == handles.data.mepmin{i}{j});
        handles.data.pos_max_init = handles.data.pos_max;
        handles.data.pos_min_init = handles.data.pos_min;
        
        % plots
        if isempty(aux_trigger)
            handles.hmeps{i, j} = nan;
            handles.hmepmean{i, j} = nan;
            handles.hmepmax{i, j} = nan;
            handles.hmepmin{i, j} = nan;
        else
            handles.hmepmean{i, j} = plot(handles.data.xs{i}{j}(1:handles.data.s1-handles.data.s0)+handles.MEPStart,...
                handles.data.mepmean{i}{j},'Visible','off');
            handles.hmepmax{i, j} = plot(handles.data.xs{i}{j}(handles.data.pos_max{i, j}(1))+handles.MEPStart,...
                handles.data.mepmax{i}{j},'+r','Visible','off');
            handles.hmepmin{i, j} = plot(handles.data.xs{i}{j}(handles.data.pos_min{i, j}(1))+handles.MEPStart,...
                handles.data.mepmin{i}{j},'+r','Visible','off');
        end
        
    end
    update_waitbar(handles,j/(length(handles.data.signal.emg_map)))
end

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

previous_data = [];

[filename, pathname, filterindex] = uiputfile({'*.xls;*.xlsx','MS Excel Files (*.xls,*.xlsx)';...
    '*.txt', 'ASCII format (*.txt)'}, 'Export data', 'processed_data.xlsx');

export_data = [{handles.data.configuration.signal_path handles.data.configuration.patient_id...
    handles.data.configuration.emg_side handles.data.configuration.angle_stim...
    handles.data.amp_pp_sim_diff handles.data.amp_rms_sim_diff handles.data.latency_sim...
    handles.data.fmed_sim_diff num2str(handles.data.cog_cluster_diff_pp)...
    handles.data.size_cluster_diff_pp handles.data.cluster_amp_pp_diff...
    handles.data.freq_cluster_pp num2str(handles.data.cog_cluster_diff_rms)...
    handles.data.size_cluster_diff_rms handles.data.cluster_amp_rms_diff...
    handles.data.freq_cluster_rms handles.data.cvelocity handles.data.correl_coef...
    handles.data.iz_row handles.data.mean_total_rms handles.data.mean_offset_rms...
    handles.data.latency_avarege}];

headers = [{'file_name'} {'subject'} {'hemisphere'} {'angle(degrees)'}...
    {'sim_amp_pp(uV)'} {'sim_amp_rms(uV)'} {'latency_sim(s)'}...
    {'sim_freq_med(Hz)'} {'cog_pp'} {'size_cluster_pp'}...
    {'cluster_amp_pp(uV)'} {'freq_cluster_pp(Hz)'} {'cog_rms'}...
    {'size_cluster_rms'} {'cluster_amp_rms(uV)'} {'freq_cluster_rms(Hz)'}...
    {'cv_central(m/s)'} {'cv_coef_correl'} {'iz(row)'} {'mean_total_rms(uV)'}...
    {'mean_offset_rms(uV)'} {'latency_avarege(ms)'}];

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

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in button_Save.
function button_Save_Callback(hObject, eventdata, handles)
% hObject    handle to button_Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

file_name_default = [num2str(handles.data.configuration.patient_id) '_' handles.data.configuration.emg_side '_' num2str(handles.data.configuration.angle_stim)];
file_name_default = strcat(file_name_default,'.mat');

update_waitbar(handles,0.5)

[file_name, file_path, index] = uiputfile({'*.mat','MAT-files (*.mat)'},...
    'Save data as...', file_name_default);

if index == 1
    data = handles.data;
    save([file_path file_name], 'data')
end

update_waitbar(handles,1)

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
for j = 1:length(handles.data.signal.emg_map)
    for i = 1:size(handles.data.signal.emg_map{handles.data.channel},2)
        
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
            if ishandle(handles.hlatencystart{i}{j})
                delete(handles.hlatencystart{i}{j})
            end
        end
        
        if isfield(handles,'hlatencystop')
            if ishandle(handles.hlatencystop{i}{j})
                delete(handles.hlatencystop{i}{j})
            end
        end
    end
    
    if isfield(handles,'hsignal')
        if ishandle(handles.hsignal{i})
            delete(handles.hsignal{i})
        end
    end
    update_waitbar(handles,j/length(handles.data.signal.emg_map))
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
        if ishandle(handles.hampthreshold{handles.data.channel, i})
            if ampthresholdvalue == 1
                set(handles.hampthreshold{handles.data.channel, i},'Visible','on')
            else
                set(handles.hampthreshold{handles.data.channel, i},'Visible','off')
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
    for i = 1:size(handles.hlatencystart{handles.data.channel}, 2)
        if ishandle(handles.hlatencystart{handles.data.channel}{i})
            if latencyvalue == 1
                set(handles.hlatencystart{handles.data.channel}{i},'Visible','on')
                set(handles.hlatencystop{handles.data.channel}{i},'Visible','on')
            else
                set(handles.hlatencystart{handles.data.channel}{i},'Visible','off')
                set(handles.hlatencystop{handles.data.channel}{i},'Visible','off')
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

if isfield(handles,'hsignal') == 1
    if signalvalue == 1
        set(handles.hsignal{handles.data.channel},'Visible','on')
    else
        set(handles.hsignal{handles.data.channel},'Visible','off')
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
    for i = 1:size(handles.hmepmean, 2)
        if ishandle(handles.hmepmean{handles.data.channel, i})
            if mepmeanvalue == 1
                
                set(handles.hmepmean{handles.data.channel, i},'Visible','on')
            else
                set(handles.hmepmean{handles.data.channel, i},'Visible','off')
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

if isfield(handles,'hmeps') == 1
    for i = 1:size(handles.hmeps, 2)
        if ishandle(handles.hmeps{handles.data.channel, i})
            if mepsvalue == 1
                
                set(handles.hmeps{handles.data.channel, i},'Visible','on')
            else
                set(handles.hmeps{handles.data.channel, i},'Visible','off')
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
    for i = 1:size(handles.hline, 2)
        if ishandle(handles.hline{handles.data.channel, i})
            if linevalue == 1
                set(handles.hline{handles.data.channel, i},'Visible','on')
            else
                set(handles.hline{handles.data.channel, i},'Visible','off')
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
    for i = 1:size(handles.htrigger, 2)
        if ishandle(handles.htrigger{handles.data.channel, i})
            if triggervalue == 1
                set(handles.htrigger{handles.data.channel, i},'Visible','on')
            else
                set(handles.htrigger{handles.data.channel, i},'Visible','off')
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
    for i = 1: length(handles.data.signal.emg_map)
        if ishandle(handles.hmepmax{handles.data.channel, i})
            if mepminmaxvalue == 1
                set(handles.hmepmax{handles.data.channel, i},'Visible','on')
                set(handles.hmepmin{handles.data.channel, i},'Visible','on')
            else
                set(handles.hmepmax{handles.data.channel, i},'Visible','off')
                set(handles.hmepmin{handles.data.channel, i},'Visible','off')
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
    for i = 1:length(handles.data.signal.emg_map)
        if ishandle(handles.hampthreshold{handles.data.channel, i})
            delete(handles.hampthreshold{handles.data.channel, i})
        end
    end
end

% set current channel plots as invisible
handles = visibleoff(handles);

amp_mono = get(handles.edit_AmplitudeThresholdMono,'String');
handles.data.amp_mono = str2double(amp_mono)/2.0;
amp_diff = get(handles.edit_AmplitudeThresholdDiff,'String');
handles.data.amp_diff = str2double(amp_diff)/2.0;

for j = 1:length(handles.data.signal.emg_map)
    axes(eval(strcat('handles.axes',num2str(j))))
    for i = 1:size(handles.data.signal.emg_map{handles.data.channel},2) % EMG mode
        
        if ~isempty(handles.data.mepmax{i}{j})
            shift_amp = (handles.data.mepmax{i}{j} + handles.data.mepmin{i}{j})/2;
            handles.data.amp_pp{i, j} = handles.data.mepmax{i}{j} - handles.data.mepmin{i}{j};
            if handles.data.amp_pp{i, j} == 0 || isnan(handles.data.amp_pp{i, j})
                handles.data.amp_pp{i, j} = NaN;
            end
        else
            handles.data.amp_pp{i, j} = NaN;
            shift_amp = 0;
        end
        
        if i == 1
            if abs(handles.data.amp_pp{i, j}) < 2*handles.data.amp_mono
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
            if ~isempty(handles.data.trigger{i, j})
                a = line([handles.MEPStart,handles.data.xs{i}{j}(handles.data.s1-handles.data.s0)+handles.MEPStart],...
                    [handles.data.amp_mono+shift_amp,handles.data.amp_mono+shift_amp],'Color','g','Visible','off');
                b = line([handles.MEPStart,handles.data.xs{i}{j}(handles.data.s1-handles.data.s0)+handles.MEPStart],...
                    [-handles.data.amp_mono+shift_amp,-handles.data.amp_mono+shift_amp],'Color','g','Visible','off');
            else
                a = NaN;
                b = NaN;
            end
        else
            if abs(handles.data.amp_pp{i, j}) < 2*handles.data.amp_diff
                if isfield(handles,'hmepmax')
                    if ishandle(handles.hmepmax{i, j})
                        delete(handles.hmepmax{i, j})
                        handles.data.mepmax{i}{j} = [];
                        handles.data.pos_max{i, j} = [];
                    end
                end
                if isfield(handles,'hmepmin')
                    if ishandle(handles.hmepmin{i, j})
                        delete(handles.hmepmin{i, j})
                        handles.data.mepmin{i}{j} = [];
                        handles.data.pos_min{i, j} = [];
                    end
                end
            end
            if ~isempty(handles.data.trigger{i, j})
                a = line([handles.MEPStart,handles.data.xs{i}{j}(handles.data.s1-handles.data.s0)+handles.MEPStart],...
                [handles.data.amp_diff+shift_amp,handles.data.amp_diff+shift_amp],'Color','g','Visible','off');
                b = line([handles.MEPStart,handles.data.xs{i}{j}(handles.data.s1-handles.data.s0)+handles.MEPStart],...
                [-handles.data.amp_diff+shift_amp,-handles.data.amp_diff+shift_amp],'Color','g','Visible','off');
            else
                a = NaN;
                b = NaN;
            end         
        end        
        handles.hampthreshold{i, j} = [a b];
    end
    update_waitbar(handles,j/length(handles.data.signal.emg_map))
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

button_Reset_Callback(hObject, eventdata, handles)

[fn, pn, index] = uigetfile({'*.mat'}, 'Selecione o arquivo .mat');

if index ~= 0
    data_path = strcat(pn,fn);
else
    return
end
        
loadfile = load (data_path);
handles.data = loadfile.data;

% Trick to round decimal numbers
handles.MEPStart = round((handles.data.s0)/handles.data.configuration.fsample*1000)/1000;

channels_names = '';

for i = 1:length(handles.data.configuration.emg_mode)
    channels_names = strcat(channels_names,num2str(i),'-',handles.data.configuration.emg_mode(i),'__');
end

set(handles.text_ChannelLegend, 'String',channels_names)
set(handles.text_Channel,'String',num2str(handles.data.channel))
set(handles.edit_id_number, 'String',num2str(handles.data.configuration.patient_id))
set(handles.edit_hemisphere_side, 'String',num2str(handles.data.configuration.emg_side))
set(handles.edit_number_angle_stim, 'String',num2str(handles.data.configuration.angle_stim))
set(handles.edit_Threshold,'String',num2str(handles.data.slope))
set(handles.edit_AmplitudeThresholdDiff,'String',num2str(2*handles.data.amp_diff))
set(handles.edit_AmplitudeThresholdMono,'String',num2str(2*handles.data.amp_mono))
set(handles.edit_MEPStart,'String', num2str(round((handles.data.s0*1000)/handles.data.configuration.fsample)))
set(handles.edit_MEPEnd,'String', num2str(round((handles.data.s1*1000)/handles.data.configuration.fsample)))

for j = 1:length(handles.data.signal.emg_map)
    for  i = 1:size(handles.data.signal.emg_map{handles.data.channel},2)
        
        axes(eval(strcat('handles.axes',num2str(j))))
        hold on
        
        handles.hsignal{i}(j) = plot(handles.data.xs{i}{j}, handles.data.signal.emg_map{j}(:,i),'Visible','off');
        
        if ~isempty(handles.data.trigger{i, j})
            if ~isempty(handles.data.mepmax{i}{j})
                shift_amp = (handles.data.mepmax{i}{j} + handles.data.mepmin{i}{j})/2;
            else
                shift_amp = 0;
            end
            aux_trigger = handles.data.trigger{i, j}(:,1);
            handles.htrigger{i,j} = plot(handles.data.xs{handles.data.channel}{j}(handles.data.trigger{i, j}(:,1)),...
                handles.data.trigger{i, j}(:,2),'ro', 'Visible', 'off');
            handles.hline{i, j} = line([handles.data.xs{i}{j}(handles.data.trigger{i, j}(:,1))',...
                handles.data.xs{i}{j}(handles.data.trigger{i, j}(:,1))'],...
                [min(handles.data.signal.emg_map{j}(:,i)),...
                max(handles.data.signal.emg_map{j}(:,i))],...
                'Color','g','Visible','off');
            cores = rand(length(aux_trigger),3);
            for k = 1:length(aux_trigger)
                handles.hmeps{i, j}(k) = plot(handles.data.xs{i}{j}(1:handles.data.s1-handles.data.s0)+handles.MEPStart,...
                    handles.data.meps{i}{j}{k},'-','color',cores(k,:),'Visible', 'off');
            end
            
            handles.hmepmean{i, j} = plot(handles.data.xs{i}{j}(1:handles.data.s1-handles.data.s0)+handles.MEPStart,...
                handles.data.mepmean{i}{j},'Visible','off');
            if ~isempty(handles.data.mepmax{i}{j})
                handles.hmepmax{i, j} = plot(handles.data.xs{i}{j}(handles.data.pos_max{i, j}(1))+handles.MEPStart,...
                    handles.data.mepmax{i}{j},'+r','Visible','off');
                handles.hmepmin{i, j} = plot(handles.data.xs{i}{j}(handles.data.pos_min{i, j}(1))+handles.MEPStart,...
                    handles.data.mepmin{i}{j},'+r','Visible','off');
            end
            if i == 1
                b = line([handles.MEPStart,handles.data.xs{i}{j}(handles.data.s1-handles.data.s0)+handles.MEPStart],...
                    [handles.data.amp_mono+shift_amp,handles.data.amp_mono+shift_amp],'Color','g','Visible','off');
                c = line([handles.MEPStart,handles.data.xs{i}{j}(handles.data.s1-handles.data.s0)+handles.MEPStart],...
                    [-handles.data.amp_mono+shift_amp,-handles.data.amp_mono+shift_amp],'Color','g','Visible','off');
            else
                b = line([handles.MEPStart,handles.data.xs{i}{j}(handles.data.s1-handles.data.s0)+handles.MEPStart],...
                    [handles.data.amp_diff+shift_amp,handles.data.amp_diff+shift_amp],'Color','g','Visible','off');
                c = line([handles.MEPStart,handles.data.xs{i}{j}(handles.data.s1-handles.data.s0)+handles.MEPStart],...
                    [-handles.data.amp_diff+shift_amp,-handles.data.amp_diff+shift_amp],'Color','g','Visible','off');
            end    
        else
            handles.hmeps{i, j} = nan;
            handles.hmepmean{i, j} = nan;
            handles.hmepmax{i, j} = nan;
            handles.hmepmin{i, j} = nan;
            handles.htrigger{i,j} = nan;
            handles.hline{i, j} = nan;
            b = nan;
            c = nan;
        end
        
        handles.hampthreshold{i, j} = [b c];
        
        if ~isempty(handles.data.latency{i}{j})
            handles.hlatencystart{i}{j} = plot(handles.data.latency{i}{j}(1,1),...
                handles.data.latency{i}{j}(1,2),'gv','Visible','off');
            handles.hlatencystop{i}{j} = plot(handles.data.latency{i}{j}(2,1),...
                handles.data.latency{i}{j}(2,2),'r^','Visible','off');
        else
            handles.hlatencystart{i}{j} = nan;
            handles.hlatencystop{i}{j} = nan;
        end
    end
    update_waitbar(handles,j/(length(handles.data.signal.emg_map)))
end

% updating checkbox values for MEP visualization
set(handles.checkbox_ViewMEPsMean,'Value',1)
set(handles.checkbox_ViewMinMax,'Value',1)
set(handles.checkbox_ViewLatency,'Value',1)

checkbox_ViewMEPsMean_Callback(hObject, eventdata, handles)
checkbox_ViewMinMax_Callback(hObject, eventdata, handles)
checkbox_ViewLatency_Callback(hObject, eventdata, handles)


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
for j = 1:length(handles.data.signal.emg_map)
    for i = 1:size(handles.data.signal.emg_map{handles.data.channel},2) % EMG mode
        if ~isempty(handles.data.mepmax{i}{j})
            handles.data.amp_pp{i, j} = handles.data.mepmax{i}{j} - handles.data.mepmin{i}{j};
            if handles.data.amp_pp{i, j} == 0 || isnan(handles.data.amp_pp{i, j})
                handles.data.amp_pp{i, j} = NaN;
                handles.data.amp_rms{i}{j} = NaN;
                handles.data.offset_rms{i}{j} = NaN;
                handles.data.fmed{i}{j} = NaN;
            end
        else
            handles.data.amp_pp{i, j} = NaN;
            handles.data.amp_rms{i}{j} = NaN;
            handles.data.offset_rms{i}{j} = NaN;
            handles.data.fmed{i}{j} = NaN;
        end
    end
end

% Peak-peak amplitude map
for i = 1:size(handles.data.amp_pp, 1)
    for j = 1:length(handles.data.amp_pp)
        aux_amp_pp_map(j) = handles.data.amp_pp{i,j};
    end
    aux_amp_pp_map_2 = cat(3,aux_amp_pp_map_2, aux_amp_pp_map);
end


% RMS and offset RMS amplitude maps
for i = 1:length(handles.data.amp_rms)
    for j = 1: length(handles.data.amp_rms{handles.data.channel})
        aux_amp_rms_map(j) = handles.data.amp_rms{i}{j};
        aux_offset_rms_map(j) = handles.data.offset_rms{i}{j};
    end
    aux_amp_rms_map_2 = cat(3,aux_amp_rms_map_2, aux_amp_rms_map);
    aux_offset_rms_map_2 = cat(3, aux_offset_rms_map_2, aux_offset_rms_map);
end

% Median frequency map
for i = 1:length(handles.data.fmed)
    for j = 1:length(handles.data.fmed{handles.data.channel})
        aux_fmed_map(j) = handles.data.fmed{i}{j};
    end
    aux_fmed_map_2 = cat(3,aux_fmed_map_2, aux_fmed_map);
end

aux_amp_pp_map_2(:,:,1)=[];
aux_amp_rms_map_2(:,:,1)=[];
aux_offset_rms_map_2(:,:,1)=[];
aux_fmed_map_2(:,:,1)=[];


for i = 1:size(handles.data.signal.emg_map{1}, 2)
    handles.data.amp_pp_map{i} = aux_amp_pp_map_2(:,:,i);
    handles.data.amp_rms_map{i} = aux_amp_rms_map_2(:,:,i);
    handles.data.offset_rms_map{i} = aux_offset_rms_map_2(:,:,i);
    handles.data.amp_fmed_map{i} = aux_fmed_map_2(:,:,i);
end

% MAP Visualization for Monopolar MEPs
% [handles.data] = MAPVisualization(handles.data, 1);
% MAP Visualization for Differential MEPs
[handles.data] = MAPVisualization(handles.data, handles.data.channel);
aux_latency_avarege = [];
% Calculus of avarege latency 2 rows away from the innervation zone

% TODO: this trick is to use hdsemg_diff - fix it for general case. This
% trick is inside MAPVisualization too
if size(handles.data.signal.emg_map{1}, 2) == 1
    latency_channel = 1;
else
    latency_channel = handles.data.channel;
end

for i=1:10
    if i == 1
        if ~isempty(handles.data.latency{latency_channel}{(handles.data.iz_row-2) + (i-1)*13})
            aux_latency_avarege = handles.data.latency{handles.data.channel}{(handles.data.iz_row-2) + (i-1)*13}(1,1);
        end
    elseif i > 1 && i < 6
        if ~isempty(handles.data.latency{latency_channel}{(handles.data.iz_row-2) + (i-1)*13})
            aux_latency_avarege = [aux_latency_avarege handles.data.latency{handles.data.channel}{(handles.data.iz_row-2) + (i-1)*13}(1,1)];
        end
    else
        if ~isempty(handles.data.latency{latency_channel}{(handles.data.iz_row+2) + (i-6)*13});
            aux_latency_avarege = [aux_latency_avarege handles.data.latency{handles.data.channel}{(handles.data.iz_row+2) + (i-6)*13}(1,1)];
        end
    end
end

handles.data.latency_avarege = mean(aux_latency_avarege);

% Calculus of avarege amplitude rms for all the used electrodes
% TODO: Fix it to the general case with various channels (mono and diff)
handles.data.mean_total_rms = nanmean(nanmean(handles.data.amp_rms_map{1}));
handles.data.mean_offset_rms = nanmean(nanmean(handles.data.offset_rms_map{1}));

% Calculus of avarege amplitude rms for the cluster electrodes
total_amp_pp_diff = 0;
total_amp_rms_diff = 0;
total_freq_cluster_pp = 0;
total_freq_cluster_rms = 0;

if size(handles.data.signal.emg_map{1}, 2) == 1
    map_channel = 1;
else
    map_channel = handles.data.channel;
end

for i = 1:size(handles.data.mep_cluster_diff_rms,1)
    total_amp_rms_diff = total_amp_rms_diff + handles.data.amp_rms_map{map_channel}(handles.data.mep_cluster_diff_rms(i,1),handles.data.mep_cluster_diff_rms(i,2));
    total_freq_cluster_rms = total_freq_cluster_rms +handles.data.fmed{map_channel}{handles.data.mep_cluster_diff_rms(i,1) + 13*(handles.data.mep_cluster_diff_rms(i,2)-1)};
end

for j = 1:size(handles.data.mep_cluster_diff_pp,1)
    total_amp_pp_diff = total_amp_pp_diff + handles.data.amp_pp_map{map_channel}(handles.data.mep_cluster_diff_pp(j,1),handles.data.mep_cluster_diff_pp(j,2));
    total_freq_cluster_pp = total_freq_cluster_pp + handles.data.fmed{map_channel}{handles.data.mep_cluster_diff_pp(j,1) + 13*(handles.data.mep_cluster_diff_pp(j,2)-1)};
end

handles.data.cluster_amp_rms_diff = total_amp_rms_diff/i;
handles.data.freq_cluster_rms = total_freq_cluster_rms/i;
handles.data.cluster_amp_pp_diff = total_amp_pp_diff/j;
handles.data.freq_cluster_pp = total_freq_cluster_pp/j;

% count = 0;
% total_amp_pp_diff = 0;
% for i = 1:size(handles.data.amp_pp_map(:,:,2),1)
%     for j = 1:size(handles.data.amp_pp_map(:,:,2),2)
%         if ~isnan(handles.data.amp_pp_map(i,j,2))
%             total_amp_pp_diff = handles.data.amp_pp_map(i,j,2)+total_amp_pp_diff;
%             count = count + 1;
%         end
%     end
% end
% handles.data.avarege_amp_pp_diff = total_amp_pp_diff/count;
% 
% count = 0;
% total_amp_rms_diff = 0;
% for i = 1:size(handles.data.amp_rms_map(:,:,2),1)
%     for j = 1:size(handles.data.amp_rms_map(:,:,2),2)
%         if ~isnan(handles.data.amp_rms_map(i,j,2))
%             total_amp_rms_diff = handles.data.amp_rms_map(i,j,2)+total_amp_rms_diff;
%             count = count + 1;
%         end
%     end
% end
% handles.data.avarege_amp_rms_diff = total_amp_rms_diff/count;

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in button_Sim_MEP.
function button_Sim_MEP_Callback(hObject, eventdata, handles)
% hObject    handle to button_Sim_MEP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Simulated MEP for a 10 mm diameter circular electrode
% The signal is simply the sum of the signal in each electrode inside the 
% virtual circle with 10 mm diameter

% p and aux2 - electrode simulated
% m - runs through one line of the matrix
% n - runs through one column of the matrix

% If I want to calculate 3 electrodes
% change p=0:1 to p=0:2 e aux=aux+4 to aux=aux+3

handles.data.mep_sim = zeros(handles.data.s1-handles.data.s0, 2);
aux = 0;
for p = 0:1
    aux2 = p+1;
    for m = 0:4
        for n = 1:5
            if m == 0 && n ==  1
                handles.data.mep_sim(:, aux2) = handles.data.mep_sim(:, aux2);
            elseif m == 0 && n == 5
                handles.data.mep_sim(:, aux2) = handles.data.mep_sim(:, aux2);
            elseif m == 4 && n == 1
                handles.data.mep_sim(:, aux2) = handles.data.mep_sim(:, aux2);
            elseif m == 4 && n == 5
                handles.data.mep_sim(:, aux2) = handles.data.mep_sim(:, aux2);
            else
                aux3 = n+p+aux + 13*m;
                handles.data.mep_sim(:, aux2) = handles.data.mep_sim(:, aux2) + handles.data.mepmean{1}{aux3}(:);
            end
        end
    end
    aux = aux + 4;
end

% handles.data.mep_sim_amp(1:2) = max(handles.data.mep_sim(:,1:2)) - min(handles.data.mep_sim(:,1:2));

% Differential MEP for the simulated circular electrodes with IED ~ 10.0 mm
% The electrodes are tangential
handles.data.mep_sim_diff = handles.data.mep_sim(:,2) - handles.data.mep_sim(:,1);
handles.data.sim_pos_min = find(handles.data.mep_sim_diff == min(handles.data.mep_sim_diff));
handles.data.sim_pos_max = find(handles.data.mep_sim_diff == max(handles.data.mep_sim_diff));

[sim_minmax sim_latency] = MEPDetect_sim(handles.data, handles.data.mep_sim_diff, handles.MEPStart,...
    handles.data.xs{1}{29}(1:handles.data.s1-handles.data.s0));

if ~isempty(sim_minmax)
    handles.data.amp_pp_sim_diff = sim_minmax(2,2) - sim_minmax(1,2);
else
    handles.data.amp_pp_sim_diff = [];
end

if ~isempty(sim_latency)
    handles.data.duration_sim = sim_latency(1,1) - sim_latency(2,1);
else
    handles.data.duration_sim = [];
end

handles.data.latency_sim = sim_latency(1,1);

[handles.data.fmed_sim_diff, handles.data.amp_rms_sim_diff, handles.data.fmean_sim_diff] = Fmed3cla(handles.data.mep_sim_diff,...
            handles.data.configuration.fsample, length(handles.data.mep_sim_diff));

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in button_ViewClear.
function button_ViewClear_Callback(hObject, eventdata, handles)
% hObject    handle to button_ViewClear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% TODO: This function does not clear the threshold lines. Fix it!

% set current channel plots as invisible
handles = visibleoff(handles);

% Update handles structure
guidata(hObject, handles);


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


function out = refreshaxes(handles)

% position of electrode on array
pos = get(gca,'Tag');
pos(1:4)=[];
pos = str2num(pos);

% channel selected
channel = handles.data.channel;

if get(handles.checkbox_ViewSignal,'Value') == 1
    % Open trigger manipulation window
    [trigger handles.data] = TriggerDetect(handles.data,handles.data.trigger{handles.data.channel, pos},handles.data.xs{channel}{pos}, pos); 
    
    if ishandle(handles.htrigger{handles.data.channel, pos})
        delete(handles.htrigger{handles.data.channel, pos});
        delete(handles.hsignal{handles.data.channel}(pos));
        delete(handles.hline{handles.data.channel, pos});
    end
    
    handles.data.trigger{handles.data.channel, pos} = trigger;
    
    if ~isempty(handles.data.trigger{handles.data.channel, pos})
        axes(eval(strcat('handles.axes',num2str(pos))));
        handles.hsignal{handles.data.channel}(pos) = plot(handles.data.xs{channel}{pos},handles.data.signal.emg_map{pos}(:,handles.data.channel));
        % red circulus
        if ~isempty(handles.data.trigger{handles.data.channel, pos})
            handles.htrigger{handles.data.channel, pos} = plot(handles.data.xs{channel}{pos}(handles.data.trigger{handles.data.channel,pos}(:,1)),...
                handles.data.trigger{handles.data.channel, pos}(:,2),'ro');
        else
            handles.htrigger{handles.data.channel, pos} = nan;
        end
        % green line
        if ~isempty(handles.data.trigger{handles.data.channel, pos})
            handles.hline{handles.data.channel, pos} = line([handles.data.xs{handles.data.channel}{pos}(handles.data.trigger{handles.data.channel, pos}(:,1))',...
                handles.data.xs{handles.data.channel}{pos}(handles.data.trigger{handles.data.channel, pos}(:,1))'],...
                [min(handles.data.signal.emg_map{pos}(:,handles.data.channel)),...
                max(handles.data.signal.emg_map{pos}(:,handles.data.channel))],...
                'Color','g');
        else
            handles.hline{handles.data.channel, pos} = nan;
        end
    end
    
    out = handles;
    
else
    %Open MEP manipulation window to set min, max, latency and coordinate on matrix
       
    [minmax latency handles.data] = MEPDetect(handles.data, handles.MEPStart,...
        handles.data.xs{channel}{pos}(1:handles.data.s1-handles.data.s0), pos);
    
    handles.data.latency{channel}{pos} = latency;
    
    % Plots of latencies
    if isfield(handles,'hlatencystart')
        if ishandle(handles.hlatencystart{channel}{pos})
            delete(handles.hlatencystart{channel}{pos})
        end
    end

    if isfield(handles,'hlatencystop')
        if ishandle(handles.hlatencystop{channel}{pos})
            delete(handles.hlatencystop{channel}{pos})
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
        if ishandle(handles.hmeps{channel, pos})
            delete(handles.hmeps{channel, pos})
        end
    end   
    if ~isempty(handles.data.trigger{channel, pos})
        aux_trigger = handles.data.trigger{channel, pos}(:,1);
    else
        aux_trigger = [];
        colors = [];
    end
    
    % Number of stimulus
    for k = 1:length(aux_trigger)
        handles.hmeps{channel, pos}(k) = plot(handles.data.xs{channel}{pos}(1:handles.data.s1-handles.data.s0)+handles.MEPStart,...
            handles.data.meps{channel}{pos}{k},'-','color',colors(k,:),'Visible','off');
    end

    mepsvalue = get(handles.checkbox_ViewMEPsMult,'Value');
    if mepsvalue == 1
        if isfield(handles,'hmeps') == 1
            if ishandle(handles.hmeps{channel, pos})
                set(handles.hmeps{channel, pos},'Visible','on')
            end
        end
    end
    
    % Plot of MEP Mean
    if isfield(handles,'hmepmean')
        if ishandle(handles.hmepmean{channel, pos})
            delete(handles.hmepmean{channel, pos})
        end
    end
    mepmeanvalue = get(handles.checkbox_ViewMEPsMean,'Value');
    handles.hmepmean{channel, pos} = plot(handles.data.xs{channel}{pos}(1:handles.data.s1-handles.data.s0)+handles.MEPStart,...
        handles.data.mepmean{channel}{pos}, 'Visible','off');
    if mepmeanvalue == 1
        if isfield(handles,'hmepmean')
            if ishandle(handles.hmepmean{channel, pos})
                set(handles.hmepmean{channel, pos},'Visible','on')
            end
        end
    end
      
    % Plots of minmax
    if isfield(handles,'hmepmax')
        if ishandle(handles.hmepmax{channel, pos})
            delete(handles.hmepmax{channel, pos})
        end
    end
    
    if isfield(handles,'hmepmin')
        if ishandle(handles.hmepmin{channel, pos})
            delete(handles.hmepmin{channel, pos})
        end
    end
    
    if isempty(minmax)
        handles.data.mepmax{channel}{pos} = [];
        handles.data.mepmin{channel}{pos} = [];
        
    else
        if minmax ~= 0
            handles.data.mepmax{channel}{pos} = max(minmax(:,2));
            handles.data.mepmin{channel}{pos} = min(minmax(:,2));
            
            aux_pos_max = minmax(:,2) == handles.data.mepmax{channel}{pos};
            pos_max = minmax(aux_pos_max,1);
            
            aux_pos_min = minmax(:,2) == handles.data.mepmin{channel}{pos};
            pos_min = minmax(aux_pos_min,1);
            
            handles.hmepmax{channel, pos} = plot(pos_max + handles.MEPStart,...
                handles.data.mepmax{channel}{pos},'+r','Visible','off');
            
            handles.hmepmin{channel, pos} = plot(pos_min + handles.MEPStart,...
                handles.data.mepmin{channel}{pos},'+r','Visible','off');
            
            mepminmaxvalue = get(handles.checkbox_ViewMinMax,'Value');
            if mepminmaxvalue
                if isfield(handles,'hmepmax')
                    if ishandle(handles.hmepmax{channel, pos}) && ishandle(handles.hmepmin{channel, pos})
                        set(handles.hmepmax{channel, pos},'Visible','on')
                        set(handles.hmepmin{channel, pos},'Visible','on')
                    end
                end
            end
        end
    end
    
    if latency ~= 0
        handles.data.duration{channel}{pos} = handles.data.latency{channel}{pos}(2,1) - handles.data.latency{channel}{pos}(1,1);
        handles.hlatencystart{channel}{pos} = plot(handles.data.latency{channel}{pos}(1,1),...
            handles.data.latency{channel}{pos}(1,2),'gv','Visible','off');
        handles.hlatencystop{channel}{pos} = plot(handles.data.latency{channel}{pos}(2,1),...
            handles.data.latency{channel}{pos}(2,2),'r^','Visible','off');
        latencyvalue = get(handles.checkbox_ViewLatency,'Value');
        if latencyvalue
            if isfield(handles,'hlatencystart')
                if ishandle(handles.hlatencystart{channel}{pos})
                    set(handles.hlatencystart{channel}{pos},'Visible','on')
                    set(handles.hlatencystop{channel}{pos},'Visible','on')
                end
            end
        end
    end
    
    if isempty(aux_trigger)
        handles.hmeps{channel, pos} = nan;
        handles.hmepmean{channel, pos} = nan;
        handles.hmepmax{channel, pos} = nan;
        handles.hmepmin{channel, pos} = nan;
    end
       
    % refreshing peak to peak map
    if ~isempty(handles.data.mepmax{channel}{pos})
        handles.data.amp_pp{channel, pos} = handles.data.mepmax{channel}{pos} - handles.data.mepmin{channel}{pos};
        if handles.data.amp_pp{channel, pos} == 0 || isnan(handles.data.amp_pp{channel, pos})
            handles.data.amp_pp{channel, pos} = NaN;
            handles.data.amp_rms{channel}{pos} = NaN;
            handles.data.offset_rms{channel}{pos} = NaN;
            handles.data.fmed{channel}{pos} = NaN;
        end
    else
        handles.data.amp_pp{channel, pos} = NaN;
        handles.data.amp_rms{channel}{pos} = NaN;
        handles.data.offset_rms{channel}{pos} = NaN;
        handles.data.fmed{channel}{pos} = NaN;
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

for i = 1:size(handles.data.signal.emg_map{handles.data.channel},2) % EMG mode
    for j = 1:length(handles.data.signal.emg_map) % Electrodes

        if isfield(handles,'hampthreshold')
            if ishandle(handles.hampthreshold{i, j})
                set(handles.hampthreshold{i, j},'Visible','off')
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
        
        % last version had an if for each i and j
        if isfield(handles,'hmepmax')
            if ishandle(handles.hmepmax{i, j})
                set(handles.hmepmax{i, j},'Visible','off')
            end
        end
        
        % last version had an if for each i and j
        if isfield(handles,'hmepmin')
            if ishandle(handles.hmepmin{i, j})
                set(handles.hmepmin{i, j},'Visible','off')
            end
        end
        
        if isfield(handles,'hlatencystart')
            if ishandle(handles.hlatencystart{i}{j})
                set(handles.hlatencystart{i}{j},'Visible','off')
            end
        end
        
        if isfield(handles,'hlatencystop')
            if ishandle(handles.hlatencystop{i}{j})
                set(handles.hlatencystop{i}{j},'Visible','off')
            end
        end
    end
    
    if isfield(handles,'hsignal')
        if ishandle(handles.hsignal{i})
            set(handles.hsignal{i},'Visible','off')
        end
    end   
    update_waitbar(handles,i/size(handles.data.signal.emg_map{handles.data.channel},2))
end

set(handles.checkbox_ViewSignal,'Value',0);
set(handles.checkbox_ViewMEPsMean,'Value',0);
set(handles.checkbox_ViewMEPsMult,'Value',0);
set(handles.checkbox_ViewLine,'Value',0);
set(handles.checkbox_ViewTrigger,'Value',0);
set(handles.checkbox_AmplitudeThreshold,'Value',0);
set(handles.checkbox_ViewMinMax,'Value',0);
set(handles.checkbox_ViewLatency,'Value',0);

handles.signalvalue = 0;
handles.mepmeanvalue = 0;
handles.mepsvalue = 0;
handles.linevalue = 0;
handles.triggervalue = 0;
handles.mepminmaxvalue = 0;
handles.ampthresholdvalue = 0;
handles.latencyvalue = 0;

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

handles.data.configuration.emg_side = get(handles.edit_hemisphere_side,'String');

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

handles.data.configuration.patient_id = str2double(get(handles.edit_id_number,'String'));

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

handles.data.configuration.angle_stim = str2double(get(handles.edit_number_angle_stim,'String'));

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
