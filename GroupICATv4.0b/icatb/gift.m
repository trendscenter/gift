function varargout = gift(varargin)
%%%%%%%%%%%%%%% Group ICA of fMRI Toolbox (GIFT) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIFT uses Independent Component Analysis to make group inferences from fMRI data
% For more information use the help button in the toolbox

try
    % add this statement to fix the button color
    feature('JavaFigures', 0);
catch
end

icatb_delete_gui({'groupica', 'eegift', 'sbm'});

% Begin initialization code - DO NOT EDIT

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @gift_OpeningFcn, ...
    'gui_OutputFcn',  @gift_OutputFcn, ...
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


% --- Executes just before gift is made visible.
function gift_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gift (see VARARGIN)


group_ica_modality = 'fmri';

setappdata(0, 'group_ica_modality', group_ica_modality);

% Choose default command line output for gift
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

%%%%%%% Object Callbacks %%%%%%%%%%%%%%%

% --- Outputs from this function are returned to the command line.
function varargout = gift_OutputFcn(hObject, eventdata, handles)
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

function about_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%icatb_directions('gift-help');

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


function html_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
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

chkRegressors = icatb_questionDialog('title', 'Select Regressors?', 'textbody', 'Do You Want To Select Regressors For Temporal Sort?');
selectedRegressors = {};
temporal_stats_betas = [];

if (chkRegressors)
    selectedRegressors = selectRegressors(param_file);
end

temporal_stats_betas = [];
if (~isempty(selectedRegressors))
    sesInfo = icatb_temporal_regress(param_file, selectedRegressors, []);
    temporal_stats_betas = sesInfo.temporal_stats_betas;
end

results = results_summary_gui('num_subjects', sesInfo.numOfSub);
results.temporal_stats_betas = temporal_stats_betas;

formatName = results.format;

%
% formatName = questdlg('Which format do you want to use for results summary?', ...
%     'File format', ...
%     'PDF', 'HTML', 'HTML');

drawnow;

if (isempty(formatName))
    error('Format is not selected');
end

GICA_PARAM_FILE = param_file;
outDir = fullfile(fileparts(param_file), [sesInfo.userInput.prefix, '_gica_results']);
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
opts.codeToEvaluate = 'icatb_gica_html_report(param_file, results);';
%publish('icatb_gica_html_report', 'outputDir', outDir, 'showCode', false, 'useNewFigure', false);
disp('Generating reults summary. Please wait ....');
drawnow;
publish('icatb_gica_html_report', opts);
clear global GICA_PARAM_FILE;

close all;

if (strcmpi(opts.format, 'html'))
    icatb_openHTMLHelpFile(fullfile(outDir, 'icatb_gica_html_report.html'));
else
    open(fullfile(outDir, 'icatb_gica_html_report.pdf'));
end

disp('Done');

function selectedRegressors = selectRegressors(param_file)
% Select regressors

load(param_file);

if (isempty(sesInfo.userInput.designMatrix.name))
    tmpSPMFiles = icatb_selectEntry('typeEntity', 'file', 'title', 'Select SPM design matrix/matrices', 'filter', 'SPM.mat', 'typeSelection', 'multiple');
    if (isempty(tmpSPMFiles))
        error('Design matrix/matrices is/are not selected');
    end
    sesInfo.userInput.designMatrix.name = tmpSPMFiles;
    drawnow;
    save(param_file, 'sesInfo');
end

spm_files = cellstr(char(sesInfo.userInput.designMatrix.name));

if (length(spm_files) > sesInfo.numOfSub)
    error('Number of design matrix/matrices should not exceed the number of subjects');
end

regressor_names = cell(1, length(spm_files));
for nF = 1:length(spm_files)
    
    load(spm_files{nF});
    names = strtrim(regexprep(cellstr(char(SPM.xX.name)),'Sn\(\d+\)',''));
    regressor_names{nF} = char(lower(names));
    
end


regressor_names = cellstr(char(regressor_names));
[dd, inds] = unique(regressor_names);
regressor_names = regressor_names(sort(inds));

dd = icatb_listdlg('promptstring', 'Select regressors of interest', 'liststring', regressor_names, 'selectionmode', 'multiple');
selectedRegressors = regressor_names(dd);