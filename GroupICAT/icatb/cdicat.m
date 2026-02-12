function varargout = cdicat(varargin)
% cdicat MATLAB code for cdicat.fig
%      cdicat, by itself, creates a new cdicat or raises the existing
%      singleton*.
%
%      H = cdicat returns the handle to a new cdicat or the handle to
%      the existing singleton*.
%
%      cdicat('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in cdicat.M with the given input arguments.
%
%      cdicat('Property','Value',...) creates a new cdicat or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before cdicat_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to cdicat_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help cdicat

% Last Modified by GUIDE v2.5 19-May-2025 23:30:47

% Begin initialization code - DO NOT EDIT
icatb_delete_gui({'groupica', 'gift', 'eegift', 'sbm', 'fnc'});

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @cdicat_OpeningFcn, ...
    'gui_OutputFcn',  @cdicat_OutputFcn, ...
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


% --- Executes just before cdicat is made visible.
function cdicat_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to cdicat (see VARARGIN)


group_ica_modality = 'CONN';

setappdata(0, 'group_ica_modality', group_ica_modality);

% Choose default command line output for cdicat
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

movegui(hObject, 'center');

% UIWAIT makes cdicat wait for user response (see UIRESUME)
% uiwait(handles.cdicat);


% --- Outputs from this function are returned to the command line.
function varargout = cdicat_OutputFcn(hObject, eventdata, handles)
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
disp('Done running Connectivity Domain Analysis');


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
