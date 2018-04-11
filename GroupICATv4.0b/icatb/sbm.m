function varargout = sbm(varargin)
% SBM M-file for sbm.fig
%      SBM, by itself, creates a new SBM or raises the existing
%      singleton*.
%
%      H = SBM returns the handle to a new SBM or the handle to
%      the existing singleton*.
%
%      SBM('Property','Value',...) creates a new SBM using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to sbm_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      SBM('CALLBACK') and SBM('CALLBACK',hObject,...) call the
%      local function named CALLBACK in SBM.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help sbm

% Last Modified by GUIDE v2.5 01-Mar-2016 14:55:11

try
    % add this statement to fix the button color
    feature('JavaFigures', 0);
catch
end

icatb_delete_gui({'groupica', 'gift', 'eegift'});

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @sbm_OpeningFcn, ...
    'gui_OutputFcn',  @sbm_OutputFcn, ...
    'gui_LayoutFcn',  [], ...
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


% --- Executes just before sbm is made visible.
function sbm_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

group_ica_modality = 'smri';

setappdata(0, 'group_ica_modality', group_ica_modality);

% Choose default command line output for sbm
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

icatb_check_path;

% move the gui at the center of the screen
movegui(hObject, 'center');

%% Turn off finite warning
warning('off', 'MATLAB:FINITE:obsoleteFunction');

%% Turn off nifti class warning when loading SPM mat files
warning('off', 'MATLAB:unknownElementsNowStruc');

warning('off', 'MATLAB:dispatcher:ShadowedMEXExtension');

warning('off', 'MATLAB:pfileOlderThanMfile');

% UIWAIT makes sbm wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = sbm_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function groupAnalysis_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Setup ICA Analysis Callback
icatb_enterParametersGUI;

function runAnalysis_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Run Analysis Callback
icatb_runAnalysis;

function analysisInfo_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Analysis Info Callback
icatb_displaySesInfo;

function dispGUI_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Display GUI Callback
icatb_displayGUI;

function utilities_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% all the strings
allStrings = get(hObject, 'string');

% get the value
getValue = get(hObject, 'value');

if ~iscell(allStrings)
    % get the selected string
    selectedString = deblank(allStrings(getValue, :));
else
    % if the selected string is cell
    selectedString = allStrings{getValue};
end

% if the selected string is other than the utilities
if getValue > 1
    % call the function
    icatb_utilities(lower(selectedString));
end

function componentExplorer_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Component Explorer callback
icatb_componentExplore;

function compositeViewer_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Composite Viewer Callback
icatb_compositeViewer;

function orthogonalViewer_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Orthogonal Viewer Callback
icatb_display_orthoviews;

function about_Callback(hObject, eventdata, handles)
% hObject    handle to about (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%icatb_directions('sbm-help');

% About Callback
icatb_titleDialog;

function exit_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Exit Callback

icatb_exit;

function instructions_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Help Callback
icatb_openHelp;

% --- Executes on button press in results_summary.
function results_summary_Callback(hObject, eventdata, handles)
% hObject    handle to results_summary (see GCBO)
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

results = results_summary_gui('num_subjects', sesInfo.numOfScans);

formatName = results.format;

drawnow;

if (isempty(formatName))
    error('Format is not selected');
end

GICA_PARAM_FILE = param_file;
outDir = fullfile(fileparts(param_file), [sesInfo.userInput.prefix, '_sbm_results']);
opts.outputDir = outDir;
opts.showCode = false;
opts.useNewFigure = false;
opts.format = lower(formatName);
opts.createThumbnail = true;
if (strcmpi(opts.format, 'pdf'))
    opt.useNewFigure = false;
end
assignin('base', 'param_file', param_file);
assignin('base', 'results', results);
opts.codeToEvaluate = 'icatb_sbm_html_report(param_file, results);';
disp('Generating reults summary. Please wait ....');
drawnow;
publish('icatb_sbm_html_report', opts);
clear global GICA_PARAM_FILE;

close all;

if (strcmpi(opts.format, 'html'))
    icatb_openHTMLHelpFile(fullfile(outDir, 'icatb_sbm_html_report.html'));
else
    open(fullfile(outDir, 'icatb_sbm_html_report.pdf'));
end

disp('Done');
