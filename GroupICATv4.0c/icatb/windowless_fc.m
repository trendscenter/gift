function varargout = windowless_fc(varargin)
% WINDOWLESS_FC MATLAB code for windowless_fc.fig
%      WINDOWLESS_FC, by itself, creates a new WINDOWLESS_FC or raises the existing
%      singleton*.
%
%      H = WINDOWLESS_FC returns the handle to a new WINDOWLESS_FC or the handle to
%      the existing singleton*.
%
%      WINDOWLESS_FC('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in WINDOWLESS_FC.M with the given input arguments.
%
%      WINDOWLESS_FC('Property','Value',...) creates a new WINDOWLESS_FC or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before windowless_fc_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to windowless_fc_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help windowless_fc

% Last Modified by GUIDE v2.5 16-Jun-2018 12:22:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @windowless_fc_OpeningFcn, ...
    'gui_OutputFcn',  @windowless_fc_OutputFcn, ...
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


% --- Executes just before windowless_fc is made visible.
function windowless_fc_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to windowless_fc (see VARARGIN)

% Choose default command line output for windowless_fc
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

movegui(hObject, 'center');

% UIWAIT makes windowless_fc wait for user response (see UIRESUME)
% uiwait(handles.windowless_fc);


% --- Outputs from this function are returned to the command line.
function varargout = windowless_fc_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in setup_wfc.
function setup_wfc_Callback(hObject, eventdata, handles)
% hObject    handle to setup_wfc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%dynamic_coherence_setup;
icatb_windowless_FC;


% --- Executes on button press in display_wfc.
function display_wfc_Callback(hObject, eventdata, handles)
% hObject    handle to display_wfc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

icatb_display_windowless_FC;

% --- Executes on button press in exit_dyn_coh.
function exit_dyn_coh_Callback(hObject, eventdata, handles)
% hObject    handle to exit_dyn_coh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(get(0, 'children'));
