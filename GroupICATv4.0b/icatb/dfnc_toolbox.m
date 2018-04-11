function varargout = dfnc_toolbox(varargin)
% DFNC_TOOLBOX M-file for dfnc_toolbox.fig
%      DFNC_TOOLBOX, by itself, creates a new DFNC_TOOLBOX or raises the existing
%      singleton*.
%
%      H = DFNC_TOOLBOX returns the handle to a new DFNC_TOOLBOX or the handle to
%      the existing singleton*.
%
%      DFNC_TOOLBOX('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DFNC_TOOLBOX.M with the given input arguments.
%
%      DFNC_TOOLBOX('Property','Value',...) creates a new DFNC_TOOLBOX or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dfnc_toolbox_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dfnc_toolbox_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help dfnc_toolbox

% Last Modified by GUIDE v2.5 27-Mar-2017 16:38:02


% if (icatb_get_matlab_version < 2006)
%     error('Mancova toolbox works from R2006a');
% end

% v = ver;
% v = lower(cellstr(char(v.Name)));
% vars = lower({'Optimization Toolbox', 'Signal Processing Toolbox', 'Image Processing Toolbox', 'Statistics Toolbox'});
%
% for nV = 1:length(vars)
%     chk = strmatch(vars{nV}, v, 'exact');
%     if (isempty(chk))
%         error('Mancova toolbox requires optimization, image processing, stats and signal processing toolboxes');
%     end
% end

try
    % add this statement to fix the button color
    feature('JavaFigures', 0);
catch
end

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @dfnc_toolbox_OpeningFcn, ...
    'gui_OutputFcn',  @dfnc_toolbox_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before dfnc_toolbox is made visible.
function dfnc_toolbox_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dfnc_toolbox (see VARARGIN)

% Choose default command line output for dfnc_toolbox
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

group_ica_modality = 'fmri';
setappdata(0, 'group_ica_modality', group_ica_modality);

icatb_check_path;

% move the gui at the center of the screen
movegui(hObject, 'center');

%% Turn off finite warning
warning('off', 'MATLAB:FINITE:obsoleteFunction');

%% Turn off nifti class warning when loading SPM mat files
warning('off', 'MATLAB:unknownElementsNowStruc');

warning('off', 'MATLAB:dispatcher:ShadowedMEXExtension');

% UIWAIT makes dfnc_toolbox wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = dfnc_toolbox_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in setup_dfnc.
function setup_dfnc_Callback(hObject, eventdata, handles)
% hObject    handle to setup_dfnc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

icatb_setup_dfnc;


% --- Executes on button press in run_dfnc.
function run_dfnc_Callback(hObject, eventdata, handles)
% hObject    handle to run_dfnc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

icatb_run_dfnc;


% --- Executes on button press in postprocess_dfnc.
function postprocess_dfnc_Callback(hObject, eventdata, handles)
% hObject    handle to postprocess_dfnc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

icatb_post_process_dfnc;


% --- Executes on button press in display.
function display_Callback(hObject, eventdata, handles)
% hObject    handle to display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


icatb_display_dfnc;

% --- Executes on button press in display.
function stats_Callback(hObject, eventdata, handles)
% hObject    handle to display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


icatb_dfnc_stats;

% --- Executes on button press in exit.
function exit_Callback(hObject, eventdata, handles)
% hObject    handle to exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(get(0, 'children'));


% --------------------------------------------------------------------
function est_clusters_menu_Callback(hObject, eventdata, handles)
% hObject    handle to est_clusters_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


icatb_optimal_clusters;

% --------------------------------------------------------------------
function close_menu_Callback(hObject, eventdata, handles)
% hObject    handle to close_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


delete(get(0, 'children'));
