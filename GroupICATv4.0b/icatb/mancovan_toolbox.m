function varargout = mancovan_toolbox(varargin)
% MANCOVAN_TOOLBOX M-file for mancovan_toolbox.fig
%      MANCOVAN_TOOLBOX, by itself, creates a new MANCOVAN_TOOLBOX or raises the existing
%      singleton*.
%
%      H = MANCOVAN_TOOLBOX returns the handle to a new MANCOVAN_TOOLBOX or the handle to
%      the existing singleton*.
%
%      MANCOVAN_TOOLBOX('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MANCOVAN_TOOLBOX.M with the given input arguments.
%
%      MANCOVAN_TOOLBOX('Property','Value',...) creates a new MANCOVAN_TOOLBOX or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mancovan_toolbox_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mancovan_toolbox_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help mancovan_toolbox

% Last Modified by GUIDE v2.5 13-Nov-2017 12:36:24


if (icatb_get_matlab_version < 2006)
    error('Mancova toolbox works from R2006a');
end

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
    'gui_OpeningFcn', @mancovan_toolbox_OpeningFcn, ...
    'gui_OutputFcn',  @mancovan_toolbox_OutputFcn, ...
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


% --- Executes just before mancovan_toolbox is made visible.
function mancovan_toolbox_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mancovan_toolbox (see VARARGIN)

% Choose default command line output for mancovan_toolbox
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

% UIWAIT makes mancovan_toolbox wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = mancovan_toolbox_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in create_mancovan_design.
function create_mancovan_design_Callback(hObject, eventdata, handles)
% hObject    handle to create_mancovan_design (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

icatb_setup_mancovan_design;


% --- Executes on button press in setup_features.
function setup_features_Callback(hObject, eventdata, handles)
% hObject    handle to setup_features (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

icatb_setup_mancovan;


% --- Executes on button press in run_mancovan.
function run_mancovan_Callback(hObject, eventdata, handles)
% hObject    handle to run_mancovan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%% Select nuisance covariates?
check = icatb_questionDialog('title', 'Remove nuisance covariates', 'textbody', 'Do you want to remove variance associated with nuisance covariates? For example, you could remove sites variance from the data before computing mancovan.');
nuisance_cov_file = '';
drawnow;
if (check)
    nuisance_cov_file = icatb_selectEntry('title', 'Select nuisance covariates file', 'typeEntity', 'file', 'typeSelection', 'single', 'filter', '*.asc;*.dat');
end
drawnow;
icatb_run_mancovan([], 3, nuisance_cov_file);


% --- Executes on button press in display.
function display_Callback(hObject, eventdata, handles)
% hObject    handle to display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


icatb_display_mancovan;

% --- Executes on button press in exit.
function exit_Callback(hObject, eventdata, handles)
% hObject    handle to exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(get(0, 'children'));


% --------------------------------------------------------------------
function fnc_domain_avg_Callback(hObject, eventdata, handles)
% hObject    handle to fnc_domain_avg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function open_menu_Callback(hObject, eventdata, handles)
% hObject    handle to open_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%% Select Mancovan file
if (~exist('param_file', 'var'))
    param_file = icatb_selectEntry('title', 'Select ICA/Mancovan Parameter File', 'typeEntity', 'file', 'typeSelection', 'single', 'filter', '*mancovan.mat');
    drawnow;
    if (isempty(param_file))
        error('Mancovan parameter file is not selected');
    end
end

icatb_setup_mancovan(param_file);

% --------------------------------------------------------------------
function close_menu_Callback(hObject, eventdata, handles)
% hObject    handle to close_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


delete(get(0, 'children'));
