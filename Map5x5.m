function varargout = Map5x5(varargin)
% MAP5X5 MATLAB code for Map5x5.fig
%      MAP5X5, by itself, creates a new MAP5X5 or raises the existing
%      singleton*.
%
%      H = MAP5X5 returns the handle to a new MAP5X5 or the handle to
%      the existing singleton*.
%
%      MAP5X5('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAP5X5.M with the given input arguments.
%
%      MAP5X5('Property','Value',...) creates a new MAP5X5 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Map5x5_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Map5x5_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Map5x5

% Last Modified by GUIDE v2.5 27-Nov-2011 22:49:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Map5x5_OpeningFcn, ...
                   'gui_OutputFcn',  @Map5x5_OutputFcn, ...
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


% --- Executes just before Map5x5 is made visible.
function Map5x5_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Map5x5 (see VARARGIN)


if length (varargin) ~= 0
    handles.data = varargin{1};
end

%--------------------------------------------------------------------------
%Possivel erro de aquisicao, passa a atribuir a freq de amostragem para
%9600 para todas as aquisicoes
handles.data.Configuracao.FrequenciaAmostragem = 9600;


name = handles.data.Configuracao.ArquivoOriginal;
tam = length(name);
name(tam-3:tam) = []


for i = length(name):-1:1
    aux_name = str2num(name(i))
    
    if length(aux_name) == 0
        break
    else
        name(i) = [];
    end
end

set(handles.figure1,'Name',name)


handles.signalvalue = 0;
handles.mepmeanvalue = 0;
handles.mepsvalue = 0;
handles.linevalue = 0;
handles.triggervalue = 0;
handles.mepminmaxvalue = 0;
handles.ampthresholdvalue = 0;


handles.data.channel = 1;
handles.data.slope = 0.5;
set(handles.edit_Threshold,'String',num2str(handles.data.slope))

nome = '';



for i = 1:length(handles.data.Sinal.Nome)
    nome = strcat(nome,num2str(i),'-',handles.data.Sinal.Nome(i),'__');
end

set(handles.text_ChannelLegend, 'String',nome)
set(handles.text_Channel, 'String',num2str(handles.data.channel))


% hwb = waitbar(0,'Carregando...')

for i = 1:length(handles.data.Sinal.Mapa)
    
    axes(eval(strcat('handles.axes',num2str(i))))
    hold on
    
    for j = 1:size(handles.data.Sinal.Mapa{handles.data.channel},2)
        handles.data.xs{j}{i} = (1:length(handles.data.Sinal.Mapa{i}(:,j)))/handles.data.Configuracao.FrequenciaAmostragem;
        handles.hsignal{j}(i) = plot(handles.data.xs{j}{i},handles.data.Sinal.Mapa{i}(:,j),'Visible','off');
    end
       
    set(eval(strcat('handles.axes',num2str(i))),'HitTest','on')
    
    update_waitbar(handles,i/length(handles.data.Sinal.Mapa))
%     waitbar(i/length(handles.data.Sinal.Mapa),hwb)
end

% close(hwb)
% clear hwb

set(handles.hsignal{1},'Visible','on')
set(handles.checkbox_ViewSignal,'Value',1)


%Gera a cell de latência
for i = 1:25
    for j = 1:5
        handles.data.latency{i}{j} = [];
        handles.hlatencystart{i}{j} = [];
        handles.hlatencystop{i}{j} = [];
    end
end


% Choose default command line output for Map5x5
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Map5x5 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Map5x5_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function update_waitbar(handles,value)
% Update waitbar
h=handles.axes_waitbar;
if ~ishandle(h);return;end;
set(h,'Visible','On');
axes(h);
cla;
patch([0,value,value,0],[0,0,1,1],'b');
axis([0,1,0,1]);
axis off;


function edit_Threshold_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Threshold as text
%        str2double(get(hObject,'String')) returns contents of edit_Threshold as a double


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


% --- Executes on button press in button_Trigger.
function button_Trigger_Callback(hObject, eventdata, handles)
% hObject    handle to button_Trigger (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


set(handles.checkbox_ViewSignal,'Value',0)
checkbox_ViewSignal_Callback(hObject, eventdata, handles)


if isfield(handles,'htrigger')
    for j = 1:length(handles.htrigger)
        if ishandle(handles.htrigger(j)) == 1
            delete(handles.htrigger(j))
        end
    end
end

if isfield(handles,'hline') == 1
    for i = 1:length(handles.hline)
        if ishandle(handles.hline{i}) == 1
            delete(handles.hline{i})
        end
    end
end

clear handles.data.slope
handles.data.slope = get(handles.edit_Threshold,'String');
handles.data.slope = str2num(handles.data.slope);



set(handles.checkbox_ViewTrigger,'Value',1)
for i = 1:length(handles.data.Sinal.Mapa)
    
    axes(eval(strcat('handles.axes',num2str(i))))
    
    handles.data.trigger{i} = BiopacTrigger(handles.data.Sinal.Mapa{i}(:,handles.data.channel),handles.data.slope);
    
    if length(handles.data.trigger{i}) ~= 0
        handles.htrigger(i) = plot(handles.data.xs{handles.data.channel}{i}(handles.data.trigger{i}(:,1)),handles.data.trigger{i}(:,2),'ro', 'Visible', 'off');
    else
        handles.htrigger(i) = nan;
    end
        
    set(eval(strcat('handles.axes',num2str(i))),'HitTest','on')
    
    
    if length(handles.data.trigger{i}) ~= 0
        a = line([handles.data.xs{handles.data.channel}{i}(handles.data.trigger{i}(:,1))',handles.data.xs{handles.data.channel}{i}(handles.data.trigger{i}(:,1))'],...
            [min(handles.data.Sinal.Mapa{i}(:,handles.data.channel)),...
            max(handles.data.Sinal.Mapa{i}(:,handles.data.channel))],...
            'Color','g','Visible','off');
        handles.hline{i} = a;
    end
    
    update_waitbar(handles,i/length(handles.data.Sinal.Mapa))
end


set(handles.checkbox_ViewSignal,'Value',1)
set(handles.checkbox_ViewTrigger,'Value',1)
set(handles.checkbox_ViewLine,'Value',1)

checkbox_ViewSignal_Callback(hObject, eventdata, handles)
checkbox_ViewTrigger_Callback(hObject, eventdata, handles)
checkbox_ViewLine_Callback(hObject, eventdata, handles)

% Update handles structure
guidata(hObject, handles);




function out = refreshaxes(handles)


pos = get(gca,'Tag');
pos(1:4)=[];
pos = str2num(pos);

channel = get(handles.text_Channel,'String');
channel = str2num(channel);


if get(handles.checkbox_ViewSignal,'Value') == 1
    
    trigger = TriggerDetect(handles.data.Sinal.Mapa{pos}(:,channel),handles.data.trigger{pos},handles.data.xs{channel}{pos});
    
    if ishandle(handles.htrigger(pos)) == 1
        delete(handles.htrigger(pos))
    end
    
    handles.data.trigger{pos} = trigger;
    
    if length(handles.data.trigger{pos}) > 0
        axes(eval(strcat('handles.axes',num2str(pos))))
        handles.htrigger(pos) = plot(handles.data.xs{channel}{pos}(handles.data.trigger{pos}(:,1)),handles.data.trigger{pos}(:,2),'ro');
    end
    
    out = handles;
    
    
else
    
    [minmax latency] = MEPDetect(handles.data.mepmean{channel}{pos}, handles.data.xs{channel}{pos}(1:handles.data.s1-handles.data.s0)+handles.MEPStart); % Descobrir o erro!!!!!!
    
    handles.data.latency{pos}{channel} = latency;
    
    if length(latency) == 0
        
         if isfield(handles,'hlatencystart') == 1
            if ishandle(handles.hlatencystart{pos}{channel}) == 1
                delete(handles.hlatencystart{pos}{channel})
            end
        end
        
        if isfield(handles,'hlatencystop') == 1
            if ishandle(handles.hlatencystop{pos}{channel}) == 1
                delete(handles.hlatencystop{pos}{channel})
            end
        end
    end
    
    if length(handles.data.latency{pos}{channel})== 2
        
        if isfield(handles,'hlatencystart') == 1
            if ishandle(handles.hlatencystart{pos}{channel}) == 1
                delete(handles.hlatencystart{pos}{channel})
            end
        end
        
        if isfield(handles,'hlatencystop') == 1
            if ishandle(handles.hlatencystop{pos}{channel}) == 1
                delete(handles.hlatencystop{pos}{channel})
            end
        end
        
        axes(gca)
        handles.hlatencystart{pos}{channel} = plot(handles.data.latency{pos}{channel}(1,1),handles.data.latency{pos}{channel}(1,2),'gv');
        handles.hlatencystop{pos}{channel} = plot(handles.data.latency{pos}{channel}(2,1),handles.data.latency{pos}{channel}(2,2),'r^');
    end
    
    if length(minmax) == 0
        handles.data.mepmax{channel}{pos} = [];
        handles.data.mepmin{channel}{pos} = [];
        
        if isfield(handles,'hmepmax') == 1
            if ishandle(handles.hmepmax{channel}(pos)) == 1
                delete(handles.hmepmax{channel}(pos))
            end
        end
        
        if isfield(handles,'hmepmin') == 1
            if ishandle(handles.hmepmin{channel}(pos)) == 1
                delete(handles.hmepmin{channel}(pos))
            end
        end
        
    else
        
        if minmax ~= 0
            
            if isfield(handles,'hmepmax') == 1
                if ishandle(handles.hmepmax{channel}(pos)) == 1
                    delete(handles.hmepmax{channel}(pos))
                end
            end
            
            if isfield(handles,'hmepmin') == 1
                if ishandle(handles.hmepmin{channel}(pos)) == 1
                    delete(handles.hmepmin{channel}(pos))
                end
            end
            
            
            handles.data.mepmax{channel}{pos} = max(minmax(:,2));
            handles.data.mepmin{channel}{pos} = min(minmax(:,2));
            
            aux_pos_max = find(minmax(:,2) == handles.data.mepmax{channel}{pos});
            pos_max = minmax(aux_pos_max,1);
            
            aux_pos_min = find(minmax(:,2) == handles.data.mepmin{channel}{pos});
            pos_min = minmax(aux_pos_min,1);
            
            
            handles.hmepmax{channel}(pos) = plot(pos_max,...
                handles.data.mepmax{channel}{pos},'+r','Visible','on');
            
            handles.hmepmin{channel}(pos) = plot(pos_min,...
                handles.data.mepmin{channel}{pos},'+r','Visible','on');
            
            set(handles.checkbox_ViewMinMax,'Value',1)
        end
        
        
    end
    
    out = handles;
end







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



% --- Executes on button press in button_Forward.
function button_Forward_Callback(hObject, eventdata, handles)
% hObject    handle to button_Forward (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

a1 = get(handles.checkbox_ViewSignal,'Value');
a2 = get(handles.checkbox_ViewMEPsMean,'Value');
a3 = get(handles.checkbox_ViewMEPsMult,'Value');
a4 = get(handles.checkbox_ViewLine,'Value');
a5 = get(handles.checkbox_ViewTrigger,'Value');
a6 = get(handles.checkbox_AmplitudeThreshold,'Value');
a7 = get(handles.checkbox_ViewMinMax,'Value');

handles = visibleoff(handles);

handles.signalvalue = a1;
handles.mepmeanvalue = a2;
handles.mepsvalue = a3;
handles.linevalue = a4;
handles.triggervalue = a5;
handles.ampthresholdvalue = a6;
handles.mepminmaxvalue = a7;

set(handles.checkbox_ViewSignal,'Value',handles.signalvalue);
set(handles.checkbox_ViewMEPsMean,'Value',handles.mepmeanvalue);
set(handles.checkbox_ViewMEPsMult,'Value',handles.mepsvalue);
set(handles.checkbox_ViewLine,'Value',handles.linevalue);
set(handles.checkbox_ViewTrigger,'Value',handles.triggervalue);
set(handles.checkbox_AmplitudeThreshold,'Value',handles.ampthresholdvalue);
set(handles.checkbox_ViewMinMax,'Value',handles.mepminmaxvalue);

if handles.data.channel >= size(handles.data.Sinal.Mapa{1},2)    
   handles.data.channel = 1; 
   set(handles.text_Channel, 'String',num2str(handles.data.channel))   

else
    handles.data.channel = handles.data.channel+1;        
    set(handles.text_Channel, 'String',num2str(handles.data.channel))      
end

checkbox_ViewSignal_Callback(hObject, eventdata, handles)
checkbox_ViewMEPsMean_Callback(hObject, eventdata, handles)
checkbox_ViewMEPsMult_Callback(hObject, eventdata, handles)
checkbox_ViewLine_Callback(hObject, eventdata, handles)
checkbox_ViewTrigger_Callback(hObject, eventdata, handles)
checkbox_AmplitudeThreshold_Callback(hObject, eventdata, handles)
checkbox_ViewMinMax_Callback(hObject, eventdata, handles)


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Map5x5 wait for user response (see UIRESUME)
% uiwait(handles.figure1);






% --- Executes on button press in button_Backward.
function button_Backward_Callback(hObject, eventdata, handles)
% hObject    handle to button_Backward (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% button_ViewClear_Callback(hObject, eventdata, handles)
% set(handles.checkbox_ViewSignal,'Value',handles.signalvalue);
% set(handles.checkbox_ViewMEPsMean,'Value',handles.mepmeanvalue);
% set(handles.checkbox_ViewMEPsMult,'Value',handles.mepsvalue);
% set(handles.checkbox_ViewLine,'Value',handles.linevalue);
% set(handles.checkbox_ViewTrigger,'Value',handles.triggervalue);

a1 = get(handles.checkbox_ViewSignal,'Value');
a2 = get(handles.checkbox_ViewMEPsMean,'Value');
a3 = get(handles.checkbox_ViewMEPsMult,'Value');
a4 = get(handles.checkbox_ViewLine,'Value');
a5 = get(handles.checkbox_ViewTrigger,'Value');
a6 = get(handles.checkbox_AmplitudeThreshold,'Value');
a7 = get(handles.checkbox_ViewMinMax,'Value');

handles = visibleoff(handles);

handles.signalvalue = a1;
handles.mepmeanvalue = a2;
handles.mepsvalue = a3;
handles.linevalue = a4;
handles.triggervalue = a5;
handles.ampthresholdvalue = a6;
handles.mepminmaxvalue = a7;

set(handles.checkbox_ViewSignal,'Value',handles.signalvalue);
set(handles.checkbox_ViewMEPsMean,'Value',handles.mepmeanvalue);
set(handles.checkbox_ViewMEPsMult,'Value',handles.mepsvalue);
set(handles.checkbox_ViewLine,'Value',handles.linevalue);
set(handles.checkbox_ViewTrigger,'Value',handles.triggervalue);
set(handles.checkbox_AmplitudeThreshold,'Value',handles.ampthresholdvalue);
set(handles.checkbox_ViewMinMax,'Value',handles.mepminmaxvalue);


if handles.data.channel <= 1
    handles.data.channel = size(handles.data.Sinal.Mapa{1},2); 
   set(handles.text_Channel, 'String',num2str(handles.data.channel))   

else
    handles.data.channel = handles.data.channel-1;        
    set(handles.text_Channel, 'String',num2str(handles.data.channel))      
end

checkbox_ViewSignal_Callback(hObject, eventdata, handles)
checkbox_ViewMEPsMean_Callback(hObject, eventdata, handles)
checkbox_ViewMEPsMult_Callback(hObject, eventdata, handles)
checkbox_ViewLine_Callback(hObject, eventdata, handles)
checkbox_ViewTrigger_Callback(hObject, eventdata, handles)
checkbox_AmplitudeThreshold_Callback(hObject, eventdata, handles)
checkbox_ViewMinMax_Callback(hObject, eventdata, handles)


% Update handles structure
guidata(hObject, handles);
    


% --- Executes on button press in button_MEPs.
function button_MEPs_Callback(hObject, eventdata, handles)
% hObject    handle to button_MEPs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.checkbox_ViewSignal,'Value',0)
set(handles.checkbox_ViewTrigger,'Value',0)
set(handles.checkbox_ViewLine,'Value',0)

checkbox_ViewSignal_Callback(hObject, eventdata, handles)
checkbox_ViewTrigger_Callback(hObject, eventdata, handles)
checkbox_ViewLine_Callback(hObject, eventdata, handles)

% hwb = waitbar(0,'Encontrando os MEPs...');
handles.MEPStart = round(str2num(get(handles.edit_MEPStart,'String')))/1000;
handles.data.s0 = round(str2num(get(handles.edit_MEPStart,'String'))*handles.data.Configuracao.FrequenciaAmostragem/1000);
handles.data.s1 = round(str2num(get(handles.edit_MEPEnd,'String'))*handles.data.Configuracao.FrequenciaAmostragem/1000);

off0 = -1500;
off1 = -10;
contwb = 1;

for i = 1:size(handles.data.Sinal.Mapa{handles.data.channel},2) % Channel
    for j = 1:length(handles.data.Sinal.Mapa) % Position of EMT
        
        
        axes(eval(strcat('handles.axes',num2str(j))))
        hold on
        aux_channel = handles.data.Sinal.Mapa{j}(:,i);
        
        aux_meps = zeros(handles.data.s1-handles.data.s0,1);
        
        if length(handles.data.trigger{j}) > 0
            aux_trigger = handles.data.trigger{j}(:,1);
        else
            aux_trigger = []
        end
        cores = rand(length(aux_trigger),3);
     
        
        for k = 1:length(aux_trigger) % Number of stimulus
            
            aux_offset = aux_channel(aux_trigger(k)+off0:aux_trigger(k)+off1-1);
            offset{i}{j}(k) = mean(aux_offset);
            handles.data.meps{i}{j}{k} = aux_channel(aux_trigger(k)+handles.data.s0:aux_trigger(k)+handles.data.s1-1) - offset{i}{j}(k);
            handles.hmeps{i}{j}(k) = plot(handles.data.xs{i}{j}(1:handles.data.s1-handles.data.s0)+handles.MEPStart,handles.data.meps{i}{j}{k},'-','color',cores(k,:),'Visible', 'off');
            aux_meps = [aux_meps,handles.data.meps{i}{j}{k}];
            
        end
        handles.data.mepmean{i}{j} = mean(aux_meps,2);
        handles.data.mepmax{i}{j} = max(handles.data.mepmean{i}{j});
        handles.data.mepmax_bkp = handles.data.mepmax;
        handles.data.mepmin{i}{j} = min(handles.data.mepmean{i}{j});
        handles.data.mepmin_bkp = handles.data.mepmin;
        
        handles.hmepmean{i}(j) = plot(handles.data.xs{i}{j}(1:handles.data.s1-handles.data.s0)+handles.MEPStart,handles.data.mepmean{i}{j},'Visible','off');
        
        pos_max = find(handles.data.mepmean{i}{j} == handles.data.mepmax{i}{j});
        pos_min = find(handles.data.mepmean{i}{j} == handles.data.mepmin{i}{j});
        
        handles.hmepmax{i}(j) = plot(handles.data.xs{i}{j}(pos_max(1)),...            
            handles.data.mepmax{i}{j},'+r','Visible','off');
       
        
        handles.hmepmin{i}(j) = plot(handles.data.xs{i}{j}(pos_min(1)),...
            handles.data.mepmin{i}{j},'+r','Visible','off');
        
        
        update_waitbar(handles,contwb/(5*length(handles.data.Sinal.Mapa)))
        contwb = contwb + 1;
%         waitbar(i/size(handles.data.Sinal.Mapa{handles.data.channel},2),hwb,strcat('Carregando:',num2str(i),'---->',num2str(j),'---->',num2str(k)));
    end
end

% close (hwb)
% clear hwb

set(handles.checkbox_ViewMEPsMean,'Value',1)
set(handles.checkbox_ViewMEPsMult,'Value',1)

checkbox_ViewMEPsMean_Callback(hObject, eventdata, handles)
checkbox_ViewMEPsMult_Callback(hObject, eventdata, handles)

% Update handles structure
guidata(hObject, handles);





% --- Executes on button press in button_Export.
function button_Export_Callback(hObject, eventdata, handles)
% hObject    handle to button_Export (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
name = handles.data.Configuracao.ArquivoOriginal;
tam = length(name);
name(tam-3:tam) = []


for i = length(name):-1:1
    aux_name = str2num(name(i))
    
    if length(aux_name) == 0
        break
    else
        name(i) = [];
    end
end

name = strcat(name,'_MAP')
name_mat = strcat(name,'.mat');


[file,path,index] = uiputfile('*.txt','Selecione o arquivo para os mapas de TMS:',name);

savefile_mat = strcat(path,name_mat);
maps = handles.data.AmplitudeMap;
save(savefile_mat,'maps');

if index == 1
    for i = 1:size(handles.data.AmplitudeMap,3)
        filename = strcat(file(1:length(file)-4),num2str(i),...
            file(length(file)-3:length(file)));
        savefile = strcat(path,filename);
        map = handles.data.AmplitudeMap(:,:,i);
        
        save(savefile,'map','-ASCII')
        update_waitbar(handles,i/size(handles.data.AmplitudeMap,3))
    end
end


% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in button_Save.
function button_Save_Callback(hObject, eventdata, handles)
% hObject    handle to button_Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
name = handles.data.Configuracao.ArquivoOriginal;
tam = length(name);
name(tam-3:tam) = []


for i = length(name):-1:1
    aux_name = str2num(name(i))
    
    if length(aux_name) == 0
        break
    else
        name(i) = [];
    end
end

name = strcat(name,'.mat');

[file,path,index] = uiputfile(name);

if index == 1
    savefile = strcat(path,file);
    data = handles.data;
    update_waitbar(handles,0.5)
    pause(0.1)
    save(savefile,'data')
    update_waitbar(handles,1)
end


% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in button_Reset.
function button_Reset_Callback(hObject, eventdata, handles)
% hObject    handle to button_Reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



for i = 1:size(handles.data.Sinal.Mapa{handles.data.channel},2) % Channel
    for j = 1:length(handles.data.Sinal.Mapa) % Position of EMT
        
        if i == 1
            
            if isfield(handles,'htrigger')
                if ishandle(handles.htrigger(j))
                    delete(handles.htrigger(j))
                end
            end
            
            if isfield(handles,'hline')
                if ishandle(handles.hline{j})
                    delete(handles.hline{j})
                end
            end
            
            if isfield(handles,'hampthreshold')
                if ishandle(handles.hampthreshold{j})
                    delete(handles.hampthreshold{j})
                end
            end            
        end
        
                
        if isfield(handles,'hmeps')
            if ishandle(handles.hmeps{i}{j})
                delete(handles.hmeps{i}{j})
            end
        end
               
    end
    
    if isfield(handles,'hsignal')
        if ishandle(handles.hsignal{i})
            delete(handles.hsignal{i})
        end
    end
    
    if isfield(handles,'hmepmean')
        if ishandle(handles.hmepmean{i})
            delete(handles.hmepmean{i})
        end
    end
    
    if isfield(handles,'hmepmax')
        if ishandle(handles.hmepmax{i})
            delete(handles.hmepmax{i})
        end
    end
    
    if isfield(handles,'hmepmin')
        if ishandle(handles.hmepmin{i})
            delete(handles.hmepmin{i})
        end
    end
    update_waitbar(handles,i/size(handles.data.Sinal.Mapa{handles.data.channel},2))
end

for i = 1:length(handles.data.Sinal.Mapa)
    axes(eval(strcat('handles.axes',num2str(i))))
    cla
    update_waitbar(handles,i/length(handles.data.Sinal.Mapa))
end
    
set(handles.checkbox_ViewSignal,'Value',0);
set(handles.checkbox_ViewMEPsMean,'Value',0);
set(handles.checkbox_ViewMEPsMult,'Value',0);
set(handles.checkbox_ViewLine,'Value',0);
set(handles.checkbox_ViewTrigger,'Value',0);
set(handles.checkbox_AmplitudeThreshold,'Value',0);
set(handles.checkbox_ViewMinMax,'Value',0);


clear handles.hsignal handles.htrigger handles.hline handles.hampthreshold...
    handles.hmeps handles.hmepmean handles.hmepmax handles.hmepmin

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in checkbox_ViewSignal.
function checkbox_ViewSignal_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_ViewSignal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_ViewSignal



signalvalue = get(handles.checkbox_ViewSignal,'Value');
handles.signalvalue = signalvalue;

if signalvalue == 1
    if isfield(handles,'hsignal') == 1
        set(handles.hsignal{handles.data.channel},'Visible','on')
    end
else
    if isfield(handles,'hsignal') == 1
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

%mepmean{i}(j)

mepmeanvalue = get(handles.checkbox_ViewMEPsMean,'Value');
handles.mepmeanvalue = mepmeanvalue;

if mepmeanvalue == 1
    if isfield(handles,'hmepmean') == 1
        set(handles.hmepmean{handles.data.channel},'Visible','on')
    end
else
    if isfield(handles,'hmepmean') == 1
        set(handles.hmepmean{handles.data.channel},'Visible','off')
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


mepsvalue = get(handles.checkbox_ViewMEPsMult,'Value');
handles.mepsvalue = mepsvalue;

if mepsvalue == 1
    if isfield(handles,'hmeps') == 1
        for i = 1:length(handles.hmeps{handles.data.channel})
            set(handles.hmeps{handles.data.channel}{i},'Visible','on')
        end
    end
else
    if isfield(handles,'hmeps') == 1
        for i = 1:length(handles.hmeps{handles.data.channel})
            set(handles.hmeps{handles.data.channel}{i},'Visible','off')
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

linevalue = get(handles.checkbox_ViewLine,'Value');
handles.linevalue = linevalue;

if linevalue == 1
    if isfield(handles,'hline') == 1
        for i = 1:length(handles.hline)
            set(handles.hline{i},'Visible','on')
        end
    end
else
    if isfield(handles,'hline') == 1
        for i = 1:length(handles.hline)
            set(handles.hline{i},'Visible','off')
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

triggervalue = get(handles.checkbox_ViewTrigger,'Value');
handles.triggervalue = triggervalue;

if triggervalue == 1
    if isfield(handles,'htrigger')
        for i = 1:length(handles.htrigger)
            if ishandle(handles.htrigger(i))
                set(handles.htrigger(i),'Visible','on')
            end
        end
    end
else
    if isfield(handles,'htrigger')
        for i = 1:length(handles.htrigger)
            if ishandle(handles.htrigger(i))
                set(handles.htrigger(i),'Visible','off')
            end
        end
    end
end

% Update handles structure
guidata(hObject, handles);



function edit_MEPEnd_Callback(hObject, eventdata, handles)
% hObject    handle to edit_MEPEnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_MEPEnd as text
%        str2double(get(hObject,'String')) returns contents of edit_MEPEnd as a double


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



function edit_MEPStart_Callback(hObject, eventdata, handles)
% hObject    handle to edit_MEPStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_MEPStart as text
%        str2double(get(hObject,'String')) returns contents of edit_MEPStart as a double


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


% --- Executes on button press in button_ViewClear.
function button_ViewClear_Callback(hObject, eventdata, handles)
% hObject    handle to button_ViewClear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = visibleoff(handles);

% Update handles structure
guidata(hObject, handles);


function out = visibleoff(handles)

for i = 1:size(handles.data.Sinal.Mapa{handles.data.channel},2) % Channel
    for j = 1:length(handles.data.Sinal.Mapa) % Position of EMT
        
        if i == 1
            
            if isfield(handles,'htrigger')
                if ishandle(handles.htrigger(j))
                    set(handles.htrigger(j),'Visible','off')
                end
            end
            
            if isfield(handles,'hline')
                if ishandle(handles.hline{j})
                    set(handles.hline{j},'Visible','off')
                end
            end
            
            if isfield(handles,'hampthreshold')
                if ishandle(handles.hampthreshold{j})
                    set(handles.hampthreshold{j},'Visible','off')
                end
            end            
        end
        
                
        if isfield(handles,'hmeps')
            if ishandle(handles.hmeps{i}{j})
                set(handles.hmeps{i}{j},'Visible','off')
            end
        end
        
        if isfield(handles,'hmepmax')
            if i <= length(handles.hmepmax)
                if j <= length(handles.hmepmax{1})
                    if ishandle(handles.hmepmax{i}(j))
                        set(handles.hmepmax{i}(j),'Visible','off')
                    end
                end
            end
        end
        
        if isfield(handles,'hmepmin')
            if i <= length(handles.hmepmin)
                if j <= length(handles.hmepmin{1})
                    if ishandle(handles.hmepmin{i}(j))
                        set(handles.hmepmin{i}(j),'Visible','off')
                    end
                end
            end
        end
        
    end
    
    if isfield(handles,'hsignal')
        if ishandle(handles.hsignal{i})
            set(handles.hsignal{i},'Visible','off')
        end
    end
    
    if isfield(handles,'hmepmean')
        if ishandle(handles.hmepmean{i})
            set(handles.hmepmean{i},'Visible','off')
        end
    end
    
    if isfield(handles,'hlatencystart') == 1
        if ishandle(handles.hlatencystart) == 1
            set(handles.hlatencystart,'Visible','off')
        end
    end
    
    if isfield(handles,'hlatencystop') == 1
        if ishandle(handles.hlatencystop) == 1
            set(handles.hlatencystop,'Visible','off')
        end
    end
    
    
    update_waitbar(handles,i/size(handles.data.Sinal.Mapa{handles.data.channel},2))
end

set(handles.checkbox_ViewSignal,'Value',0);
set(handles.checkbox_ViewMEPsMean,'Value',0);
set(handles.checkbox_ViewMEPsMult,'Value',0);
set(handles.checkbox_ViewLine,'Value',0);
set(handles.checkbox_ViewTrigger,'Value',0);
set(handles.checkbox_AmplitudeThreshold,'Value',0);
set(handles.checkbox_ViewMinMax,'Value',0);

handles.carai = 'carai';

handles.signalvalue = 0;
handles.mepmeanvalue = 0;
handles.mepsvalue = 0;
handles.linevalue = 0;
handles.triggervalue = 0;
handles.mepminmaxvalue = 0;
handles.ampthresholdvalue = 0;

out = handles;

% for i = 1:length(handles.hsignal)
%     set(handles.hsignal{i},'Visible','off')    
%     set(handles.hmepmean{i},'Visible','off')    
%     
%     for j = 1:length(handles.hmeps{i})
%         set(handles.hmeps{i}{j},'Visible','off')
%     end
% end
% 
% for j = 1:length(handles.hmeps{i})
%     set(handles.hline{j},'Visible','off')
% end
% 
% set(handles.htrigger,'Visible','off')


% --- Executes on button press in checkbox_ViewMinMax.
function checkbox_ViewMinMax_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_ViewMinMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_ViewMinMax

mepminmaxvalue = get(handles.checkbox_ViewMinMax,'Value');
handles.mepminmaxvalue = mepminmaxvalue;


if mepminmaxvalue == 1
    for i = 1: length(handles.data.Sinal.Mapa)
        if isfield(handles,'hmepmax') == 1
            if ishandle(handles.hmepmax{handles.data.channel}(i))
                set(handles.hmepmax{handles.data.channel}(i),'Visible','on')
                set(handles.hmepmin{handles.data.channel}(i),'Visible','on')
            end
        end
        
        if isfield(handles,'hlatencystart') == 1
            if ishandle(handles.hlatencystart{i}{handles.data.channel})
                set(handles.hlatencystart{i}{handles.data.channel},'Visible','on')
                set(handles.hlatencystop{i}{handles.data.channel},'Visible','on')
            end
        end       
        
        
    end
else
    for i = 1: length(handles.data.Sinal.Mapa)
        if isfield(handles,'hmepmax') == 1
            if ishandle(handles.hmepmax{handles.data.channel}(i))
                set(handles.hmepmax{handles.data.channel}(i),'Visible','off')
                set(handles.hmepmin{handles.data.channel}(i),'Visible','off')
            end
        end 
        
        if isfield(handles,'hlatencystart') == 1
            if ishandle(handles.hlatencystart{i}{handles.data.channel})
                set(handles.hlatencystart{i}{handles.data.channel},'Visible','off')
                set(handles.hlatencystop{i}{handles.data.channel},'Visible','off')
            end
        end 
        
    end
end



% Update handles structure
guidata(hObject, handles);



function edit_AmplitudeThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to edit_AmplitudeThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_AmplitudeThreshold as text
%        str2double(get(hObject,'String')) returns contents of edit_AmplitudeThreshold as a double


% --- Executes during object creation, after setting all properties.
function edit_AmplitudeThreshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_AmplitudeThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_AmplitudeThresholdMode_Callback(hObject, eventdata, handles)
% hObject    handle to edit_AmplitudeThresholdMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_AmplitudeThresholdMode as text
%        str2double(get(hObject,'String')) returns contents of edit_AmplitudeThresholdMode as a double


% --- Executes during object creation, after setting all properties.
function edit_AmplitudeThresholdMode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_AmplitudeThresholdMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in button_AmplitudeThreshold.
function button_AmplitudeThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to button_AmplitudeThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.data.mepmax{handles.data.channel} = handles.data.mepmax_bkp{handles.data.channel};
handles.data.mepmin{handles.data.channel} = handles.data.mepmin_bkp{handles.data.channel};



if isfield(handles,'hampthreshold') == 1
    for i = 1: length(handles.data.Sinal.Mapa)
        if ishandle(handles.hampthreshold{i}) == 1
            delete(handles.hampthreshold{i})
        end
    end
end

amp = get(handles.edit_AmplitudeThreshold,'String');
amp = str2num(amp);


for i = 1: length(handles.data.Sinal.Mapa)
    
    if isfield(handles,'hmepmax') == 1
        if ishandle(handles.hmepmax{handles.data.channel}(i))
            delete(handles.hmepmax{handles.data.channel}(i))            
        end
    end
    
    if isfield(handles,'hmepmin') == 1
        if ishandle(handles.hmepmin{handles.data.channel}(i))
            delete(handles.hmepmin{handles.data.channel}(i))            
        end
    end
    
    axes(eval(strcat('handles.axes',num2str(i))))
    pos_max = find(handles.data.mepmean{handles.data.channel}{i} == handles.data.mepmax{handles.data.channel}{i});
    pos_min = find(handles.data.mepmean{handles.data.channel}{i} == handles.data.mepmin{handles.data.channel}{i});
    
    handles.hmepmax{handles.data.channel}(i) = plot(handles.data.xs{handles.data.channel}{i}(pos_max(1))+handles.MEPStart,...
        handles.data.mepmax{handles.data.channel}{i},'+r','Visible','on');
    
    
    handles.hmepmin{handles.data.channel}(i) = plot(handles.data.xs{handles.data.channel}{i}(pos_min(1))+handles.MEPStart,...
        handles.data.mepmin{handles.data.channel}{i},'+r','Visible','on');
    
    if handles.data.mepmax{handles.data.channel}{i} < amp|handles.data.mepmin{handles.data.channel}{i} > -amp
        
        if isfield(handles,'hmepmax') == 1
            if ishandle(handles.hmepmax{handles.data.channel}(i))
                delete(handles.hmepmax{handles.data.channel}(i))
                handles.data.mepmax{handles.data.channel}{i} = [];
            end
        end
        
        if isfield(handles,'hmepmin') == 1
            if ishandle(handles.hmepmin{handles.data.channel}(i))
                delete(handles.hmepmin{handles.data.channel}(i))
                handles.data.mepmin{handles.data.channel}{i} = [];
            end
        end
        
    end
    
    
    axes(eval(strcat('handles.axes',num2str(i))))
    a = line([0,handles.data.xs{handles.data.channel}{i}(handles.data.s1-handles.data.s0)],[amp,amp],'Color','g','Visible','on');
    b = line([0,handles.data.xs{handles.data.channel}{i}(handles.data.s1-handles.data.s0)],[-amp,-amp],'Color','g','Visible','on');
    handles.hampthreshold{i} = [a b];
    update_waitbar(handles,i/length(handles.data.Sinal.Mapa))
end

set(handles.checkbox_AmplitudeThreshold,'Value',1)
set(handles.checkbox_ViewMinMax,'Value',1);

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in checkbox_AmplitudeThreshold.
function checkbox_AmplitudeThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_AmplitudeThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_AmplitudeThreshold


ampthresholdvalue = get(handles.checkbox_AmplitudeThreshold,'Value');
handles.ampthresholdvalue = ampthresholdvalue;

if ampthresholdvalue == 1
    if isfield(handles,'hampthreshold') == 1
        for i = 1:length(handles.hampthreshold)
            set(handles.hampthreshold{i},'Visible','on')
        end
    end
else
    if isfield(handles,'hampthreshold') == 1
        for i = 1:length(handles.hampthreshold)
            set(handles.hampthreshold{i},'Visible','off')
        end
    end
end

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in button_Map.
function button_Map_Callback(hObject, eventdata, handles)
% hObject    handle to button_Map (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

aux_AmplitudeMap = zeros(5);
handles.data.AmplitudeMap = zeros(5);

for i = 1 : length(handles.data.mepmax)
    for j = 1: length(handles.data.mepmax{i})
        if length(handles.data.mepmax{i}{j}) ~= 0
            handles.data.Amplitude{i}(j) = handles.data.mepmax{i}{j} - handles.data.mepmin{i}{j};
        else
            handles.data.Amplitude{i}(j) = 0;
        end
        aux_AmplitudeMap(j) = handles.data.Amplitude{i}(j);
        
    end
    handles.data.AmplitudeMap = cat(3,handles.data.AmplitudeMap,aux_AmplitudeMap');
end

handles.data.AmplitudeMap(:,:,1)=[];

figure
for i = 1:4
    subplot(2,2,i)
    imagesc(handles.data.AmplitudeMap(:,:,i))
    colorbar
end

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
update_waitbar(handles,0.5)
pause(0.1)
        
loadfile = load (data_path);
handles.data = loadfile.data;

set(handles.edit_Threshold,'String',num2str(handles.data.slope))
set(handles.text_Channel,'String',num2str(handles.data.channel))

contwb = 1;

for i = 1:size(handles.data.Sinal.Mapa{handles.data.channel},2) % Channel
    for j = 1:length(handles.data.Sinal.Mapa) % Position of EMT
        
        
        axes(eval(strcat('handles.axes',num2str(j))))
        hold on
        
        handles.hsignal{i}(j) = plot(handles.data.xs{i}{j},...
            handles.data.Sinal.Mapa{j}(:,i),'Visible','off');
        
        if i == 1
            handles.htrigger(j) = plot(handles.data.xs{handles.data.channel}{j}(handles.data.trigger{j}(:,1)),handles.data.trigger{j}(:,2),'ro', 'Visible', 'off');
            
            if length(handles.data.trigger{j}) ~= 0
                a = line([handles.data.xs{handles.data.channel}{j}(handles.data.trigger{j}(:,1))',handles.data.xs{handles.data.channel}{j}(handles.data.trigger{j}(:,1))'],...
                    [min(handles.data.Sinal.Mapa{j}(:,handles.data.channel)),...
                    max(handles.data.Sinal.Mapa{j}(:,handles.data.channel))],...
                    'Color','g','Visible','off');
                handles.hline{j} = a;
            end
            
            amp = get(handles.edit_AmplitudeThreshold,'String');
            amp = str2num(amp);
            
            b = line([0,handles.data.xs{handles.data.channel}{j}(handles.data.s1-handles.data.s0)],[amp,amp],'Color','g','Visible','off');
            c = line([0,handles.data.xs{handles.data.channel}{j}(handles.data.s1-handles.data.s0)],[-amp,-amp],'Color','g','Visible','off');
            handles.hampthreshold{j} = [b c];
            
        end
        
        aux_channel = handles.data.Sinal.Mapa{j}(:,i);
        
        aux_meps = zeros(handles.data.s1-handles.data.s0,1);
        aux_trigger = handles.data.trigger{j}(:,1);
        cores = rand(length(aux_trigger),3);
        
        
        for k = 1:length(aux_trigger) % Number of stimulus
            
            handles.hmeps{i}{j}(k) = plot(handles.data.xs{i}{j}(1:handles.data.s1-handles.data.s0),handles.data.meps{i}{j}{k},'-','color',cores(k,:),'Visible', 'off');
            aux_meps = [aux_meps,handles.data.meps{i}{j}{k}];
        end
        
        handles.hmepmean{i}(j) = plot(handles.data.xs{i}{j}(1:handles.data.s1-handles.data.s0),handles.data.mepmean{i}{j},'Visible','off');
        
        if length(handles.data.mepmax{i}{j}) ~= 0
            pos_max = find(handles.data.mepmean{i}{j} == handles.data.mepmax{i}{j});
            
            if length(pos_max) > 0
                handles.hmepmax{i}(j) = plot(handles.data.xs{i}{j}(pos_max(1)),...
                    handles.data.mepmax{i}{j},'+r','Visible','off');
            else
                handles.hmepmax{i}(j) = nan;
            end
        end
        
        if length(handles.data.mepmin{i}{j}) ~= 0
            pos_min = find(handles.data.mepmean{i}{j} == handles.data.mepmin{i}{j});
            
            if length(pos_min) > 0
                handles.hmepmin{i}(j) = plot(handles.data.xs{i}{j}(pos_min(1)),...
                    handles.data.mepmin{i}{j},'+r','Visible','off');
            else
                handles.hmepmin{i}(j) = nan;
            end
        end
        
        update_waitbar(handles,contwb/(5*length(handles.data.Sinal.Mapa)))
        contwb = contwb + 1;
        
    end
end

set(handles.checkbox_ViewSignal,'Value',1)
checkbox_ViewSignal_Callback(hObject, eventdata, handles)


% Update handles structure
guidata(hObject, handles);
