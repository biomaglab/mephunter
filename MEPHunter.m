function varargout = MEPHunter(varargin)
% MEPHUNTER M-file for MEPHunter.fig
%      MEPHUNTER, by itself, creates a new MEPHUNTER or raises the existing
%      singleton*.
%
%      H = MEPHUNTER returns the handle to a new MEPHUNTER or the handle to
%      the existing singleton*.
%
%      MEPHUNTER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MEPHUNTER.M with the given input arguments.
%
%      MEPHUNTER('Property','Value',...) creates a new MEPHUNTER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MEPHunter_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MEPHunter_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%
% Created by André Salles Cunha Peres 02/09/2011 

% Edit the above text to modify the response to help MEPHunter

% Last Modified by GUIDE v2.5 31-Aug-2011 13:16:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MEPHunter_OpeningFcn, ...
                   'gui_OutputFcn',  @MEPHunter_OutputFcn, ...
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


% --- Executes just before MEPHunter is made visible.
function MEPHunter_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MEPHunter (see VARARGIN)
% MEPHnter_logo = load('MEPHunter_logo.mat');
MEPHnter_logo = load('mephunter_logo.mat');
axes(handles.axes1)
image(MEPHnter_logo.logo)
pause(2)
MEPHunter_Main
close



% Choose default command line output for MEPHunter
%handles.output = hObject;

% Update handles structure
%guidata(hObject, handles);

% UIWAIT makes MEPHunter wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MEPHunter_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% uiwait(hObject);
% handles = guidata(hObject);
% 
% delete(hObject);

% Get default command line output from handles structure
%varargout{1} = handles.output;
