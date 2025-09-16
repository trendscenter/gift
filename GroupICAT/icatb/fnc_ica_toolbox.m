function varargout = fnc_ica_toolbox(varargin)
% FNC_ICA_TOOLBOX MATLAB code for fnc_ica_toolbox.fig
%      FNC_ICA_TOOLBOX, by itself, creates a new FNC_ICA_TOOLBOX or raises the existing
%      singleton*.
%
%      H = FNC_ICA_TOOLBOX returns the handle to a new FNC_ICA_TOOLBOX or the handle to
%      the existing singleton*.
%
%      FNC_ICA_TOOLBOX('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FNC_ICA_TOOLBOX.M with the given input arguments.
%
%      FNC_ICA_TOOLBOX('Property','Value',...) creates a new FNC_ICA_TOOLBOX or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fnc_ica_toolbox_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fnc_ica_toolbox_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help fnc_ica_toolbox

% Last Modified by GUIDE v2.5 19-Oct-2024 10:36:57

% Begin initialization code - DO NOT EDIT
icatb_delete_gui({'groupica', 'gift', 'eegift', 'sbm'});

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @fnc_ica_toolbox_OpeningFcn, ...
    'gui_OutputFcn',  @fnc_ica_toolbox_OutputFcn, ...
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


% --- Executes just before fnc_ica_toolbox is made visible.
function fnc_ica_toolbox_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to fnc_ica_toolbox (see VARARGIN)


group_ica_modality = 'FNC';

setappdata(0, 'group_ica_modality', group_ica_modality);

% Choose default command line output for fnc_ica_toolbox
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

movegui(hObject, 'center');

% UIWAIT makes fnc_ica_toolbox wait for user response (see UIRESUME)
% uiwait(handles.fnc_ica_toolbox);


% --- Outputs from this function are returned to the command line.
function varargout = fnc_ica_toolbox_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in setup_fnc_ica.
function setup_fnc_ica_Callback(hObject, eventdata, handles)
% hObject    handle to setup_fnc_ica (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

icatb_setup_fnc_ica;


% --- Executes on button press in display_fnc_ica.
function display_fnc_ica_Callback(hObject, eventdata, handles)
% hObject    handle to display_fnc_ica (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

param_file = icatb_selectEntry('typeEntity', 'file', 'title', 'Select Parameter File', 'filter', '*ica*parameter*.mat');

load(param_file);

drawnow;

if (isempty(param_file))
    error('Parameter file is not selected for display');
end

formatName = questdlg('Select results format', 'Results format', 'HTML', 'PDF', 'HTML');

results.format = formatName;
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

icatb_run_fnc_ica;


% --- Executes on button press in stats.
function stats_Callback(hObject, eventdata, handles)
% hObject    handle to stats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

icatb_stats_loadings;