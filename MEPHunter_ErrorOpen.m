function varargout = MEPHunter_ErrorOpen(varargin)
% MEPHUNTER_ERROROPEN M-file for MEPHunter_ErrorOpen.fig
%      MEPHUNTER_ERROROPEN, by itself, creates a new MEPHUNTER_ERROROPEN or raises the existing
%      singleton*.
%
%      H = MEPHUNTER_ERROROPEN returns the handle to a new MEPHUNTER_ERROROPEN or the handle to
%      the existing singleton*.
%
%      MEPHUNTER_ERROROPEN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MEPHUNTER_ERROROPEN.M with the given input arguments.
%
%      MEPHUNTER_ERROROPEN('Property','Value',...) creates a new MEPHUNTER_ERROROPEN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MEPHunter_ErrorOpen_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MEPHunter_ErrorOpen_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MEPHunter_ErrorOpen

% Last Modified by GUIDE v2.5 26-Aug-2011 17:48:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MEPHunter_ErrorOpen_OpeningFcn, ...
                   'gui_OutputFcn',  @MEPHunter_ErrorOpen_OutputFcn, ...
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


% --- Executes just before MEPHunter_ErrorOpen is made visible.
function MEPHunter_ErrorOpen_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MEPHunter_ErrorOpen (see VARARGIN)

% Choose default command line output for MEPHunter_ErrorOpen
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MEPHunter_ErrorOpen wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MEPHunter_ErrorOpen_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in button_ErrorOpen.
function button_ErrorOpen_Callback(hObject, eventdata, handles)
% hObject    handle to button_ErrorOpen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume
close


