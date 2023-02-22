function varargout = dfdc_toolbox(varargin)
% DFDC_TOOLBOX MATLAB code for dfdc_toolbox.fig
%      DFDC_TOOLBOX, by itself, creates a new DFDC_TOOLBOX or raises the existing
%      singleton*.
%
%      H = DFDC_TOOLBOX returns the handle to a new DFDC_TOOLBOX or the handle to
%      the existing singleton*.
%
%      DFDC_TOOLBOX('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DFDC_TOOLBOX.M with the given input arguments.
%
%      DFDC_TOOLBOX('Property','Value',...) creates a new DFDC_TOOLBOX or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dfdc_toolbox_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dfdc_toolbox_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help dfdc_toolbox

% Last Modified by GUIDE v2.5 25-Sep-2019 15:44:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @dfdc_toolbox_OpeningFcn, ...
    'gui_OutputFcn',  @dfdc_toolbox_OutputFcn, ...
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


% --- Executes just before dfdc_toolbox is made visible.
function dfdc_toolbox_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dfdc_toolbox (see VARARGIN)

% Choose default command line output for dfdc_toolbox
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

movegui(hObject, 'center');

% UIWAIT makes dfdc_toolbox wait for user response (see UIRESUME)
% uiwait(handles.dfdc_toolbox);


% --- Outputs from this function are returned to the command line.
function varargout = dfdc_toolbox_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in setup_dfdc.
function setup_dfdc_Callback(hObject, eventdata, handles)
% hObject    handle to setup_dfdc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%dynamic_coherence_setup;
icatb_setup_dfdc;


% --- Executes on button press in stats_dfdc.
function stats_dfdc_Callback(hObject, eventdata, handles)
% hObject    handle to stats_dfdc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

icatb_stats_dfdc;

% --- Executes on button press in exit_dyn_coh.
function exit_dyn_coh_Callback(hObject, eventdata, handles)
% hObject    handle to exit_dyn_coh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(get(0, 'children'));
