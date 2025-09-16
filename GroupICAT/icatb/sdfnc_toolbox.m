function varargout = sdfnc_toolbox(varargin)
% SDFNC_TOOLBOX M-file for sdfnc_toolbox.fig
%      SDFNC_TOOLBOX, by itself, creates a new SDFNC_TOOLBOX or raises the existing
%      singleton*.
%
%      H = SDFNC_TOOLBOX returns the handle to a new SDFNC_TOOLBOX or the handle to
%      the existing singleton*.
%
%      SDFNC_TOOLBOX('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SDFNC_TOOLBOX.M with the given input arguments.
%
%      SDFNC_TOOLBOX('Property','Value',...) creates a new SDFNC_TOOLBOX or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before sdfnc_toolbox_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to sdfnc_toolbox_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help sdfnc_toolbox

% Last Modified by GUIDE v2.5 31-Mar-2014 16:12:43


try
    % add this statement to fix the button color
    feature('JavaFigures', 0);
catch
end

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @sdfnc_toolbox_OpeningFcn, ...
    'gui_OutputFcn',  @sdfnc_toolbox_OutputFcn, ...
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


% --- Executes just before sdfnc_toolbox is made visible.
function sdfnc_toolbox_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to sdfnc_toolbox (see VARARGIN)

% Choose default command line output for sdfnc_toolbox
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

% UIWAIT makes sdfnc_toolbox wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = sdfnc_toolbox_OutputFcn(hObject, eventdata, handles)
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

icatb_setup_spatial_dfnc;


function open_sdfnc_Callback(hObject, eventdata, handles)

icatb_setup_spatial_dfnc;

% [filename, pathname] = uigetfile('*_sdfnc.mat', 'Select spatial dFNC Info File');
%
% if (~filename)
%     error('Spatial dfnc file is not selected');
% end
%
% param_file = fullfile(pathname, filename);
%
% icatb_setup_spatial_dfnc(param_file);

% --- Executes on button press in run_dfnc.
function run_dfnc_Callback(hObject, eventdata, handles)
% hObject    handle to run_dfnc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

icatb_run_spatial_dfnc;


% --- Executes on button press in postprocess_dfnc.
function postprocess_dfnc_Callback(hObject, eventdata, handles)
% hObject    handle to postprocess_dfnc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

icatb_post_process_spatial_dfnc;


% --- Executes on button press in display.
function display_Callback(hObject, eventdata, handles)
% hObject    handle to display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


icatb_spatial_dfnc_results;

% param_file = icatb_selectEntry('title', 'Select spatial dFNC Parameter File', 'typeEntity', 'file', 'typeSelection', 'single', 'filter', '*_sdfnc.mat');
% drawnow;
% if (isempty(param_file))
%     error('ICA/sdFNC parameter file is not selected');
% end
% 
% global SDFNC_PARAM_FILE;
% SDFNC_PARAM_FILE = param_file;
% 
% load(param_file);
% prefix = sdfncInfo.userInput.prefix;
% outDir = fullfile(fileparts(param_file), [prefix, '_sdfnc_html']);
% 
% publish('icatb_spatial_dfnc_results', 'outputDir', outDir, 'showCode', false, 'useNewFigure', false);
% 
% icatb_openHTMLHelpFile(fullfile(outDir, 'icatb_spatial_dfnc_results.html'));

% --- Executes on button press in exit.
function exit_Callback(hObject, eventdata, handles)
% hObject    handle to exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(get(0, 'children'));
