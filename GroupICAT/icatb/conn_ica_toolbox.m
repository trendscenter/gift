function varargout = conn_ica_toolbox(varargin)
% CONN_ICA_TOOLBOX MATLAB code for conn_ica_toolbox.fig
%      CONN_ICA_TOOLBOX, by itself, creates a new CONN_ICA_TOOLBOX or raises the existing
%      singleton*.
%
%      H = CONN_ICA_TOOLBOX returns the handle to a new CONN_ICA_TOOLBOX or the handle to
%      the existing singleton*.
%
%      CONN_ICA_TOOLBOX('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CONN_ICA_TOOLBOX.M with the given input arguments.
%
%      CONN_ICA_TOOLBOX('Property','Value',...) creates a new CONN_ICA_TOOLBOX or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before conn_ica_toolbox_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to conn_ica_toolbox_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help conn_ica_toolbox

% Last Modified by GUIDE v2.5 19-May-2025 23:30:47

% Begin initialization code - DO NOT EDIT
icatb_delete_gui({'groupica', 'gift', 'eegift', 'sbm', 'fnc'});

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @conn_ica_toolbox_OpeningFcn, ...
    'gui_OutputFcn',  @conn_ica_toolbox_OutputFcn, ...
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


% --- Executes just before conn_ica_toolbox is made visible.
function conn_ica_toolbox_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to conn_ica_toolbox (see VARARGIN)


group_ica_modality = 'CONN';

setappdata(0, 'group_ica_modality', group_ica_modality);

% Choose default command line output for conn_ica_toolbox
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

movegui(hObject, 'center');

% UIWAIT makes conn_ica_toolbox wait for user response (see UIRESUME)
% uiwait(handles.conn_ica_toolbox);


% --- Outputs from this function are returned to the command line.
function varargout = conn_ica_toolbox_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in setup_conn_ica.
function setup_conn_ica_Callback(hObject, eventdata, handles)
% hObject    handle to setup_conn_ica (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

icatb_setup_conn_ica;


% --- Executes on button press in display_fnc_ica.
function display_fnc_ica_Callback(hObject, eventdata, handles)
% hObject    handle to display_fnc_ica (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



icatb_defaults;
global PARAMETER_INFO_MAT_FILE;
global GICA_PARAM_FILE;

filterP = ['*', PARAMETER_INFO_MAT_FILE, '*.mat'];
param_file = icatb_selectEntry('typeEntity', 'file', 'title', 'Select Parameter File', 'filter', filterP);

if (isempty(param_file))
    error('Parameter file is not selected for analysis');
end

load(param_file);

if (~exist('sesInfo', 'var'))
    error('Selected file is not a valid parameter file');
end


temporal_stats_betas = [];

results = results_summary_gui('num_subjects', sesInfo.numOfSub);

formatName = results.format;
results.formatName = formatName;

drawnow;

giftPath = fileparts(which('gift.m'));

resultsFile = fullfile(fileparts(param_file), [sesInfo.userInput.prefix, '_tmp_results_struct.mat']);
save(resultsFile, 'results');

% Run second matlab instance (matlab is installed on system)
disp('Generating summary with new matlab process ...');
commandStr = ['matlab -nodesktop -nosplash -r "addpath(genpath(''', giftPath, ''')); icatb_report_generator(''', param_file, ...
    ''', ''', resultsFile, ''');exit;', '"'];
[status, message] = system(commandStr);

if (status ~= 0)
    error(message);
end

disp('Done');


% param_file = icatb_selectEntry('typeEntity', 'file', 'title', 'Select Parameter File', 'filter', '*ica*parameter*.mat');
% 
% load(param_file);
% 
% drawnow;
% 
% if (isempty(param_file))
%     error('Parameter file is not selected for display');
% end
% 
% formatName = questdlg('Select results format', 'Results format', 'HTML', 'PDF', 'HTML');
% 
% results.format = formatName;
% results.formatName = formatName;
% 
% drawnow;
% 
% giftPath = fileparts(which('gift.m'));
% 
% resultsFile = fullfile(fileparts(param_file), [sesInfo.userInput.prefix, '_tmp_results_struct.mat']);
% save(resultsFile, 'results');
% 
% % Run second matlab instance (matlab is installed on system)
% disp('Generating summary with new matlab process ...');
% commandStr = ['matlab -nodesktop -nosplash -r "addpath(genpath(''', giftPath, ''')); icatb_report_generator(''', param_file, ...
%     ''', ''', resultsFile, ''');exit;', '"'];
% [status, message] = system(commandStr);
% 
% if (status ~= 0)
%     error(message);
% end
% 
% disp('Done');

%icatb_report_generator(param_file);

% --- Executes on button press in exit.
function exit_Callback(hObject, eventdata, handles)
% hObject    handle to exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(get(0, 'children'));


% --- Executes on button press in run_fnc_ica.
function run_fnc_ica_Callback(hObject, eventdata, handles)
% hObject    handle to run_fnc_ica (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%icatb_run_conn_ica;
icatb_runAnalysis;


% --- Executes on button press in stats.
function stats_Callback(hObject, eventdata, handles)
% hObject    handle to stats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%icatb_stats_loadings;


% --- Executes on button press in import_data.
function import_data_Callback(hObject, eventdata, handles)
% hObject    handle to import_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

icatb_import_data_conn_ica;


% --------------------------------------------------------------------
function display_tools_Callback(hObject, eventdata, handles)
% hObject    handle to display_tools (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function image_viewer_Callback(hObject, eventdata, handles)
% hObject    handle to image_viewer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

icatb_image_viewer;
