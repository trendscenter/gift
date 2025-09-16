function varargout = dynamic_functional_conn(varargin)
% DYNAMIC_FUNCTIONAL_CONN MATLAB code for dynamic_functional_conn.fig
%      DYNAMIC_FUNCTIONAL_CONN, by itself, creates a new DYNAMIC_FUNCTIONAL_CONN or raises the existing
%      singleton*.
%
%      H = DYNAMIC_FUNCTIONAL_CONN returns the handle to a new DYNAMIC_FUNCTIONAL_CONN or the handle to
%      the existing singleton*.
%
%      DYNAMIC_FUNCTIONAL_CONN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DYNAMIC_FUNCTIONAL_CONN.M with the given input arguments.
%
%      DYNAMIC_FUNCTIONAL_CONN('Property','Value',...) creates a new DYNAMIC_FUNCTIONAL_CONN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dynamic_functional_conn_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dynamic_functional_conn_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help dynamic_functional_conn

% Last Modified by GUIDE v2.5 23-Jan-2020 22:58:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @dynamic_functional_conn_OpeningFcn, ...
                   'gui_OutputFcn',  @dynamic_functional_conn_OutputFcn, ...
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


% --- Executes just before dynamic_functional_conn is made visible.
function dynamic_functional_conn_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dynamic_functional_conn (see VARARGIN)

% Choose default command line output for dynamic_functional_conn
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

movegui(hObject, 'center');

% UIWAIT makes dynamic_functional_conn wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = dynamic_functional_conn_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in spatial_dfc.
function spatial_dfc_Callback(hObject, eventdata, handles)
% hObject    handle to spatial_dfc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sdfnc_toolbox;

% --- Executes on button press in spatial_chronnectome.
function spatial_chronnectome_Callback(hObject, eventdata, handles)
% hObject    handle to spatial_chronnectome (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

spatial_chronnectome;

% --- Executes on button press in spatial_dynamics_hierarchy.
function spatial_dynamics_hierarchy_Callback(hObject, eventdata, handles)
% hObject    handle to spatial_dynamics_hierarchy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

spatial_dynamics_hierarchy;

% --- Executes on button press in temporal_dfnc.
function temporal_dfnc_Callback(hObject, eventdata, handles)
% hObject    handle to temporal_dfnc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dfnc_toolbox;

% --- Executes on button press in windowless_fc.
function windowless_fc_Callback(hObject, eventdata, handles)
% hObject    handle to windowless_fc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

windowless_fc;

% --- Executes on button press in temporal_dfc.
function temporal_dfc_Callback(hObject, eventdata, handles)
% hObject    handle to temporal_dfc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dyn_fc_toolbox;

% --- Executes on button press in dynamic_coherence.
function dynamic_coherence_Callback(hObject, eventdata, handles)
% hObject    handle to dynamic_coherence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dynamic_coherence;
