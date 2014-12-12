function varargout = MAPLines_raw(varargin)
% MAPLINES_RAW MATLAB code for MAPLines_raw.fig
%      MAPLINES_RAW, by itself, creates a new MAPLINES_RAW or raises the existing
%      singleton*.
%
%      H = MAPLINES_RAW returns the handle to a new MAPLINES_RAW or the handle to
%      the existing singleton*.
%
%      MAPLINES_RAW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAPLINES_RAW.M with the given input arguments.
%
%      MAPLINES_RAW('Property','Value',...) creates a new MAPLINES_RAW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MAPLines_raw_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MAPLines_raw_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MAPLines_raw

% Last Modified by GUIDE v2.5 24-Oct-2013 08:53:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MAPLines_raw_OpeningFcn, ...
                   'gui_OutputFcn',  @MAPLines_raw_OutputFcn, ...
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


% --- Executes just before MAPLines_raw is made visible.
function MAPLines_raw_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MAPLines_raw (see VARARGIN)

if ~isempty (varargin)
    handles.data = varargin{1};
    handles.signal = varargin{2};
end

% center the figure window on the screen
movegui(hObject, 'center');

if isfield(handles.data, 'iz_row')
    if isempty(handles.data.iz_row)
        handles.data.iz_row = 7;
        set(handles.edit_iz_row,'String', num2str(handles.data.iz_row));
    else
        set(handles.edit_iz_row,'String', num2str(handles.data.iz_row));
    end
else
    handles.data.iz_row = 7;
    set(handles.edit_iz_row,'String', num2str(handles.data.iz_row));
end

handles.amp_map = cell(1, handles.signal.n_conditions);

% Rearranging variables
averaged_meps = cell(handles.signal.n_conditions, 1);
averaged_meps_raw = cell(handles.signal.n_conditions, 1);
max_amps = zeros(handles.signal.n_conditions,1);

for i = 1:handles.signal.n_conditions
    averaged_meps{i} = nan(size(handles.data.mepmean_bkp{i}{2}, 1), handles.signal.n_channels);
    averaged_meps_raw{i} = nan(size(handles.data.mepmean_bkp{i}{2}, 1), handles.signal.n_channels);
    for j = 1:handles.signal.n_channels
        if isempty(handles.data.mepmean_bkp{i}{j}) || isempty(handles.data.mepmax{i}{j})
            averaged_meps{i}(:,j) = nan(size(handles.data.mepmean_bkp{i}{2}, 1),1);
        else
            averaged_meps{i}(:,j) = handles.data.mepmean_bkp{i}{j};
        end
        
        if isempty(handles.data.mepmean{i}{j})
            averaged_meps_raw{i}(:,j) = nan(size(handles.data.mepmean_bkp{i}{2}, 1),1);
        else
            averaged_meps_raw{i}(:,j) = handles.data.mepmean_bkp{i}{j};
        end
    end
    max_amps(i) = max(abs(averaged_meps{i}(:)));
end

% Plotting meps' line
for i = 1:handles.signal.n_conditions
    axes(eval(strcat('handles.axes',num2str(i))));
    for cols = 1:5
        line([1:size(averaged_meps_raw{i},1)] + (cols-1)*size(averaged_meps_raw{i},1),...
            averaged_meps_raw{i}(:,[1:13] + (cols-1)*13)/max(max_amps) + repmat(13:-1:1,size(averaged_meps_raw{i},1),1),...
            'parent', gca,'color', 'b', 'LineWidth', 2)
    end

    set(gca,'xlim',[0 (cols)*size(averaged_meps_raw{i},1)],'ylim',[0 14],...
        'xtick', size(averaged_meps_raw{i},1)*([1:cols] - 1/2),'xticklabel',[],...
        'ytick',1:13,'yticklabel', [],'ydir','normal','box','on','xgrid','on','ygrid','on')
    if i == 1 || i == 5
        if i == 5
            set(gca,'xticklabel',1:cols,'yticklabel', 13:-1:1)
        else
            set(gca, 'yticklabel', 13:-1:1)
        end
    end
    if i > 5
        set(gca,'xticklabel',1:cols)
    end
    scale = sprintf(' (%.2f microV/Div)', max(max_amps));
    set(eval(strcat('handles.text_angle',num2str(i))), 'String',...
         [num2str(handles.signal.angle_stim(i)) 'º' scale])
end

handles.output = handles.data.iz_row;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MAPLines_raw wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MAPLines_raw_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiwait(hObject);
handles = guidata(hObject);
delete(hObject);

% Get default command line output from handles structure
varargout{1} = handles.output;



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

figure_name_default = [num2str(handles.signal.patient_id) '_' handles.signal.emg_side '_lines'];

[figure_name, figure_path] = uiputfile({'*.emf','Enhanced metafile (*.emf)';...
    '*.bmp','Windows bitmap (*.bmp)'},...
    'Save figure as...', figure_name_default);

saveas(handles.figure1, [figure_path figure_name])

% Update handles structure
guidata(hObject, handles);



function edit_iz_row_Callback(hObject, eventdata, handles)
% hObject    handle to edit_iz_row (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_iz_row as text
%        str2double(get(hObject,'String')) returns contents of edit_iz_row as a double

handles.data.iz_row = str2double(get(handles.edit_iz_row, 'String'));
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
