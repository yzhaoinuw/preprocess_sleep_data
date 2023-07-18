function varargout = Bin2EdfGUI(varargin)
% BIN2EDFGUI MATLAB code for Bin2EdfGUI.fig
%      BIN2EDFGUI, by itself, creates a new BIN2EDFGUI or raises the existing
%      singleton*.
%
%      H = BIN2EDFGUI returns the handle to a new BIN2EDFGUI or the handle to
%      the existing singleton*.
%
%      BIN2EDFGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BIN2EDFGUI.M with the given input arguments.
%
%      BIN2EDFGUI('Property','Value',...) creates a new BIN2EDFGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Bin2EdfGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Bin2EdfGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Bin2EdfGUI

% Last Modified by GUIDE v2.5 24-Apr-2019 18:57:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Bin2EdfGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @Bin2EdfGUI_OutputFcn, ...
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


% --- Executes just before Bin2EdfGUI is made visible.
function Bin2EdfGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Bin2EdfGUI (see VARARGIN)

% Choose default command line output for Bin2EdfGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Bin2EdfGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Bin2EdfGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in OneEdfperBin.
function OneEdfperBin_Callback(hObject, eventdata, handles)
% hObject    handle to OneEdfperBin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of OneEdfperBin


% --- Executes on button press in LoadExpFile.
function LoadExpFile_Callback(hObject, eventdata, handles)
% hObject    handle to LoadExpFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Bin2EDF([],get(handles.OneEdfperBin,'value'));
