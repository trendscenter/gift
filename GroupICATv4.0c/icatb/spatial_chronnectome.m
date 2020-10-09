function varargout = spatial_chronnectome(varargin)
% SPATIAL_CHRONNECTOME MATLAB code for spatial_chronnectome.fig
%      SPATIAL_CHRONNECTOME, by itself, creates a new SPATIAL_CHRONNECTOME or raises the existing
%      singleton*.
%
%      H = SPATIAL_CHRONNECTOME returns the handle to a new SPATIAL_CHRONNECTOME or the handle to
%      the existing singleton*.
%
%      SPATIAL_CHRONNECTOME('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPATIAL_CHRONNECTOME.M with the given input arguments.
%
%      SPATIAL_CHRONNECTOME('Property','Value',...) creates a new SPATIAL_CHRONNECTOME or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before spatial_chronnectome_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to spatial_chronnectome_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help spatial_chronnectome

% Last Modified by GUIDE v2.5 02-Aug-2020 22:44:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @spatial_chronnectome_OpeningFcn, ...
    'gui_OutputFcn',  @spatial_chronnectome_OutputFcn, ...
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


% --- Executes just before spatial_chronnectome is made visible.
function spatial_chronnectome_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to spatial_chronnectome (see VARARGIN)

% Choose default command line output for spatial_chronnectome
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

movegui(hObject, 'center');

% UIWAIT makes spatial_chronnectome wait for user response (see UIRESUME)
% uiwait(handles.spatial_chronnectome);


% --- Outputs from this function are returned to the command line.
function varargout = spatial_chronnectome_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in setup_schron.
function setup_schron_Callback(hObject, eventdata, handles)
% hObject    handle to setup_schron (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%dynamic_coherence_setup;
icatb_setup_spatial_chronnectome;


% --- Executes on button press in postprocess_schron.
function postprocess_schron_Callback(hObject, eventdata, handles)
% hObject    handle to postprocess_schron (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

icatb_postprocess_spatial_chronnectome;

% --- Executes on button press in exit_dyn_coh.
function exit_dyn_coh_Callback(hObject, eventdata, handles)
% hObject    handle to exit_dyn_coh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(get(0, 'children'));


% --- Executes on button press in display_schronn.
function display_schronn_Callback(hObject, eventdata, handles)
% hObject    handle to display_schronn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

icatb_display_spatial_chronnectome;


% --- Executes on button press in stats.
function stats_Callback(hObject, eventdata, handles)
% hObject    handle to stats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

icatb_spatial_chronn_stats_gui;
