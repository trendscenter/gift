function varargout = dynamic_coherence(varargin)
% DYNAMIC_COHERENCE MATLAB code for dynamic_coherence.fig
%      DYNAMIC_COHERENCE, by itself, creates a new DYNAMIC_COHERENCE or raises the existing
%      singleton*.
%
%      H = DYNAMIC_COHERENCE returns the handle to a new DYNAMIC_COHERENCE or the handle to
%      the existing singleton*.
%
%      DYNAMIC_COHERENCE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DYNAMIC_COHERENCE.M with the given input arguments.
%
%      DYNAMIC_COHERENCE('Property','Value',...) creates a new DYNAMIC_COHERENCE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dynamic_coherence_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dynamic_coherence_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help dynamic_coherence

% Last Modified by GUIDE v2.5 14-Oct-2016 10:58:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @dynamic_coherence_OpeningFcn, ...
    'gui_OutputFcn',  @dynamic_coherence_OutputFcn, ...
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


% --- Executes just before dynamic_coherence is made visible.
function dynamic_coherence_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dynamic_coherence (see VARARGIN)

% Choose default command line output for dynamic_coherence
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

movegui(hObject, 'center');

% UIWAIT makes dynamic_coherence wait for user response (see UIRESUME)
% uiwait(handles.dynamic_coherence);


% --- Outputs from this function are returned to the command line.
function varargout = dynamic_coherence_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in setup_dyn_coh.
function setup_dyn_coh_Callback(hObject, eventdata, handles)
% hObject    handle to setup_dyn_coh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%dynamic_coherence_setup;
icatb_setup_dyn_coh;


% --- Executes on button press in display_dyn_coh.
function display_dyn_coh_Callback(hObject, eventdata, handles)
% hObject    handle to display_dyn_coh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

icatb_display_dynamic_coherence;

% --- Executes on button press in exit_dyn_coh.
function exit_dyn_coh_Callback(hObject, eventdata, handles)
% hObject    handle to exit_dyn_coh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(get(0, 'children'));
