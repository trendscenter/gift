function varargout = spatial_dynamics_hierarchy(varargin)
% SPATIAL_DYNAMICS_HIERARCHY MATLAB code for spatial_dynamics_hierarchy.fig
%      SPATIAL_DYNAMICS_HIERARCHY, by itself, creates a new SPATIAL_DYNAMICS_HIERARCHY or raises the existing
%      singleton*.
%
%      H = SPATIAL_DYNAMICS_HIERARCHY returns the handle to a new SPATIAL_DYNAMICS_HIERARCHY or the handle to
%      the existing singleton*.
%
%      SPATIAL_DYNAMICS_HIERARCHY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPATIAL_DYNAMICS_HIERARCHY.M with the given input arguments.
%
%      SPATIAL_DYNAMICS_HIERARCHY('Property','Value',...) creates a new SPATIAL_DYNAMICS_HIERARCHY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before spatial_dynamics_hierarchy_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to spatial_dynamics_hierarchy_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help spatial_dynamics_hierarchy

% Last Modified by GUIDE v2.5 04-Mar-2020 16:02:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @spatial_dynamics_hierarchy_OpeningFcn, ...
    'gui_OutputFcn',  @spatial_dynamics_hierarchy_OutputFcn, ...
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


% --- Executes just before spatial_dynamics_hierarchy is made visible.
function spatial_dynamics_hierarchy_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to spatial_dynamics_hierarchy (see VARARGIN)

% Choose default command line output for spatial_dynamics_hierarchy
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

movegui(hObject, 'center');

% UIWAIT makes spatial_dynamics_hierarchy wait for user response (see UIRESUME)
% uiwait(handles.spatial_dynamics_hierarchy);


% --- Outputs from this function are returned to the command line.
function varargout = spatial_dynamics_hierarchy_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in setup_sdh.
function setup_sdh_Callback(hObject, eventdata, handles)
% hObject    handle to setup_sdh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%dynamic_coherence_setup;
icatb_setup_sdh;


% --- Executes on button press in display_sdh.
function display_sdh_Callback(hObject, eventdata, handles)
% hObject    handle to display_sdh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

icatb_display_sdh;

% --- Executes on button press in exit_dyn_coh.
function exit_dyn_coh_Callback(hObject, eventdata, handles)
% hObject    handle to exit_dyn_coh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(get(0, 'children'));


% --- Executes on button press in post_processing.
function post_processing_Callback(hObject, eventdata, handles)
% hObject    handle to post_processing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

icatb_postprocess_sdh;