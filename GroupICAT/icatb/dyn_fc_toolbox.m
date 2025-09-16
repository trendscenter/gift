function varargout = dyn_fc_toolbox(varargin)
% DYN_FC_TOOLBOX M-file for dyn_fc_toolbox.fig
%      DYN_FC_TOOLBOX, by itself, creates a new DYN_FC_TOOLBOX or raises the existing
%      singleton*.
%
%      H = DYN_FC_TOOLBOX returns the handle to a new DYN_FC_TOOLBOX or the handle to
%      the existing singleton*.
%
%      DYN_FC_TOOLBOX('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DYN_FC_TOOLBOX.M with the given input arguments.
%
%      DYN_FC_TOOLBOX('Property','Value',...) creates a new DYN_FC_TOOLBOX or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dfnc_toolbox_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dyn_fc_toolbox_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help dyn_fc_toolbox

% Last Modified by GUIDE v2.5 24-Mar-2020 00:53:09


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
    'gui_OpeningFcn', @dyn_fc_toolbox_OpeningFcn, ...
    'gui_OutputFcn',  @dyn_fc_toolbox_OutputFcn, ...
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


% --- Executes just before dyn_fc_toolbox is made visible.
function dyn_fc_toolbox_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dyn_fc_toolbox (see VARARGIN)

% Choose default command line output for dyn_fc_toolbox
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

% UIWAIT makes dyn_fc_toolbox wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = dyn_fc_toolbox_OutputFcn(hObject, eventdata, handles)
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

icatb_setup_dfc_roi;


% --- Executes on button press in run_dfnc.
function run_dfnc_Callback(hObject, eventdata, handles)
% hObject    handle to run_dfnc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%icatb_run_dfc_roi;


% --- Executes on button press in postprocess_dfnc.
function postprocess_dfnc_Callback(hObject, eventdata, handles)
% hObject    handle to postprocess_dfnc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

icatb_postprocess_dfc_roi;


% --- Executes on button press in display.
function display_Callback(hObject, eventdata, handles)
% hObject    handle to display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


icatb_display_dfc_roi;

% --- Executes on button press in display.
function stats_Callback(hObject, eventdata, handles)
% hObject    handle to display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%icatb_dfnc_stats;

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


%icatb_optimal_clusters;

% --------------------------------------------------------------------
function close_menu_Callback(hObject, eventdata, handles)
% hObject    handle to close_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


delete(get(0, 'children'));
