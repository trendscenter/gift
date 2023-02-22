function varargout = setup_reference_ica(varargin)
%SETUP_REFERENCE_ICA MATLAB code file for setup_reference_ica.fig
%      SETUP_REFERENCE_ICA, by itself, creates a new SETUP_REFERENCE_ICA or raises the existing
%      singleton*.
%
%      H = SETUP_REFERENCE_ICA returns the handle to a new SETUP_REFERENCE_ICA or the handle to
%      the existing singleton*.
%
%      SETUP_REFERENCE_ICA('Property','Value',...) creates a new SETUP_REFERENCE_ICA using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to setup_reference_ica_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      SETUP_REFERENCE_ICA('CALLBACK') and SETUP_REFERENCE_ICA('CALLBACK',hObject,...) call the
%      local function named CALLBACK in SETUP_REFERENCE_ICA.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help setup_reference_ica

% Last Modified by GUIDE v2.5 25-Oct-2021 17:46:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @setup_reference_ica_OpeningFcn, ...
    'gui_OutputFcn',  @setup_reference_ica_OutputFcn, ...
    'gui_LayoutFcn',  [], ...
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


% --- Executes just before setup_reference_ica is made visible.
function setup_reference_ica_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for setup_reference_ica
handles.output = hObject;

set(handles.output, 'visible', 'off');

outputDir = icatb_selectEntry('typeEntity', 'directory', 'title', 'Select Analysis Output Directory');

pause(1);

set(handles.output, 'visible', 'on');

cm = uicontextmenu;
TRmenu = uimenu(cm, 'Label', 'TR', 'callback', {@openTRWindow, handles});
set(handles.TR, 'UIContextMenu', cm);

handles.outputDir = outputDir;

set(handles.template_mask_file, 'string', fullfile(fileparts(which('gift.m')), 'icatb_templates', 'RSN_28.nii'));

set(handles.template_labels_file, 'string', fullfile(fileparts(which('gift.m')), 'icatb_templates',  'RSN_28.txt'));

% Update handles structure
guidata(hObject, handles);

movegui(hObject, 'center');

% UIWAIT makes setup_reference_ica wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = setup_reference_ica_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
%varargout{1} = handles.output;
waitfor(hObject);
sesInfo = [];
appName = 'setRefAppData';
if (isappdata(0, appName))
    sesInfo = getappdata(0, appName);
    rmappdata(0, appName);
end

varargout{1} = sesInfo;

function template_mask_file_Callback(hObject, eventdata, handles)
% hObject    handle to template_mask_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of template_mask_file as text
%        str2double(get(hObject,'String')) returns contents of template_mask_file as a double


% --- Executes during object creation, after setting all properties.
function template_mask_file_CreateFcn(hObject, eventdata, handles)
% hObject    handle to template_mask_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in template_selection.
function template_selection_Callback(hObject, eventdata, handles)
% hObject    handle to template_selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

results = guidata(handles.output);
P = icatb_selectEntry('filter', '*nii;*.img', 'title', 'Select template file ...');
if (~isempty(P))
    set(results.template_mask_file, 'string', deblank(P));
end


function template_labels_file_Callback(hObject, eventdata, handles)
% hObject    handle to template_labels_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of template_labels_file as text
%        str2double(get(hObject,'String')) returns contents of template_labels_file as a double


% --- Executes during object creation, after setting all properties.
function template_labels_file_CreateFcn(hObject, eventdata, handles)
% hObject    handle to template_labels_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in template_labels.
function template_labels_Callback(hObject, eventdata, handles)
% hObject    handle to template_labels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

results = guidata(handles.output);
P = icatb_selectEntry('filter', '*.txt', 'title', 'Select atlas labels file ...', 'typeSelection', 'file', 'typeSelection', 'single');
if (~isempty(P))
    set(results.template_labels_file, 'string', deblank(P));
end

% --- Executes on button press in select_data.
function select_data_Callback(hObject, eventdata, handles)
% hObject    handle to select_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

filesInfo = [];
try
    filesInfo = handles.filesInfo;
catch
end

use_bids = icatb_questionDialog('title', 'BIDS format?', 'textbody', 'Is the data in BIDS format?');
drawnow;

if (use_bids)
    
    
    
    data_setDir = icatb_selectEntry('typeEntity', 'directory', 'title', ...
        'Select root folder for subjects and sessions');
    drawnow;
    
    
    bids_info = icatb_select_bids_params(struct('root_dir', data_setDir, 'modality_dir', 'func', 'listFiles', 0));
    input_data_file_patterns = {};
    filesInfo.filesList = input_data_file_patterns;
    filesInfo.bids_info = bids_info;
    
else
    
    filesInfo = icatb_fileSelector(filesInfo, '*.nii');
    
end
if (~isempty(filesInfo.filesList) || ~isempty(filesInfo.bids_info))
    set(handles.select_data, 'ForegroundColor', [0, 1, 0]);
end

handles.filesInfo = filesInfo;
handles.numOfDataSets = size(filesInfo.filesList, 1);
guidata(handles.output, handles);


% --- Executes on selection change in mask.
function mask_Callback(hObject, eventdata, handles)
% hObject    handle to mask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns mask contents as cell array
%        contents{get(hObject,'Value')} returns selected item from mask


maskFilter = '*.img;*.nii';
maskFigTitle = 'Select a mask file in analyze or nifti format';
fileType = 'image';
maskFile = [];

getMask = get(hObject, 'value');
getStr = cellstr(get(hObject, 'string'));
selMask = deblank(getStr{getMask});

try
    
    if (strcmpi(selMask, 'custom mask'))
        [P] = icatb_selectEntry('filter', maskFilter, 'title', maskFigTitle, 'fileType', fileType, ...
            'fileNumbers', 1);
        if ~isempty(P)
            maskFile = P;
        else
            disp('No mask file selected. Setting to default mask ...');
            set(hObject, 'value', 1);
        end
    end
    
catch
    
end

getMask = get(hObject, 'value');
getStr = cellstr(get(hObject, 'string'));
selMask = deblank(getStr{getMask});

if (isempty(maskFile))
    maskFile = selMask;
end

handles.maskFile = maskFile;

guidata(handles.output, handles);


function TR_Callback(hObject, eventdata, handles)
% hObject    handle to TR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TR as text
%        str2double(get(hObject,'String')) returns contents of TR as a double


% --- Executes during object creation, after setting all properties.
function TR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in done.
function done_Callback(hObject, eventdata, handles)
% hObject    handle to done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

results = guidata(handles.output);

outputDir = pwd;

try
    outputDir = results.outputDir;
catch
end

handles.outputDir = outputDir;

% Specify data
filesList = '';
try
    filesList = results.filesInfo.filesList;
catch
end

bids_info = [];
try
    bids_info = results.filesInfo.bids_info;
catch
end

if (isempty(bids_info))
    if (isempty(filesList))
        error('data is not specified');
    end
end

filesList = cellstr(filesList);
numOfDataSets = length(filesList);
numOfSess = 1;
try
    numOfSess = results.filesInfo.numOfSess;
catch
end
numOfSub = numOfDataSets/numOfSess;

handles.dataSelectionMethod = 4;

if (isempty(bids_info))
    input_data_file_patterns = cell(numOfSub, numOfSess);
    endTp = 0;
    for nSub = 1:numOfSub
        startTp = endTp + 1;
        endTp = endTp + numOfSess;
        input_data_file_patterns(nSub, :) = filesList(startTp:endTp);
    end
else
    input_data_file_patterns = {};
end

handles.input_data_file_patterns = input_data_file_patterns;

% Mask
maskFile = [];
try
    maskFile = results.maskFile;
catch
end
handles.maskFile = maskFile;


% TR
TR = str2num(get(results.TR, 'string'));
if (isempty(TR))
    error('Specify TR of the experiment');
end
handles.TR = TR;


prefix = get(results.prefix, 'string');
if (isempty(prefix))
    error('Specify an output prefix');
end

handles.prefix = prefix;

algoList = get(handles.algorithm, 'string');
algoVal = get(handles.algorithm, 'value');
algorithm = lower(deblank(algoList{algoVal}));

handles.algorithm = algorithm;

% Templates info
atlas_mask_file = strtrim(get(results.template_mask_file, 'string'));
if (isempty(atlas_mask_file))
    error('Specify template nifti file');
end

handles.refFiles = atlas_mask_file;

atlas_labels_file = strtrim(get(results.template_labels_file, 'string'));
if (isempty(atlas_labels_file))
    error('Specify template labels file');
end

compLabels = getCompLabels(atlas_labels_file);
networkOpts = cell(length(compLabels), 2);
for n = 1:length(compLabels)
    networkOpts{n, 1} = compLabels(n).name;
    networkOpts{n, 2} = compLabels(n).value;
end

network_summary_opts.threshold = 1.5;
% convert spatial maps to z-scores
network_summary_opts.convert_to_z = 'yes';

try
    network_summary_opts = handles.network_summary_opts;
catch
end

% format: options are html and pdf
network_summary_opts.format = 'html';


% structural file
network_summary_opts.structFile = fullfile(fileparts(which('gift.m')), 'icatb_templates', 'ch2bet_3x3x3.nii');


network_summary_opts.comp_network_names = networkOpts;


network_summary_opts.save_info = 1;

% connectivity threshold (fnc in r-values)
%network_summary_opts.conn_threshold = 0.25;

handles.network_summary_opts = network_summary_opts;

if (~isempty(bids_info))
    handles.bids_info = bids_info;
end

generateBatch(handles);

%labels = icatb_textscan(atlas_labels_file);

%setappdata(0, 'setRefAppData', sesInfo);

delete(gcbf);

drawnow;


function prefix_Callback(hObject, eventdata, handles)
% hObject    handle to prefix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of prefix as text
%        str2double(get(hObject,'String')) returns contents of prefix as a double


% --- Executes during object creation, after setting all properties.
function prefix_CreateFcn(hObject, eventdata, handles)
% hObject    handle to prefix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in analysis_type_grp.
function analysis_type_grp_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in analysis_type_grp
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in algorithm.
function algorithm_Callback(hObject, eventdata, handles)
% hObject    handle to algorithm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns algorithm contents as cell array
%        contents{get(hObject,'Value')} returns selected item from algorithm


% --- Executes during object creation, after setting all properties.
function algorithm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to algorithm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function defaults_Callback(hObject, eventdata, handles)
% hObject    handle to defaults (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

icatb_defaults;
global PREPROC_DEFAULT;
global SCALE_DEFAULT;

preproc_default = matchString(PREPROC_DEFAULT, icatb_preproc_data, 1);

scaleOptions = icatb_scaleICA;

sd = SCALE_DEFAULT;

if (isempty(sd))
    sd = 0;
end

if ~ischar(sd)
    sd = sd + 1;
end

scale_default = matchString(sd, scaleOptions, 4);

dlg_title = 'Select options for reference ica';

numParameters = 1;

inputText(numParameters).promptString =  'Select Type Of Data Pre-processing';
inputText(numParameters).answerString = char(icatb_preproc_data);
inputText(numParameters).uiType = 'popup';
inputText(numParameters).dataType = 'string';
inputText(numParameters).tag = 'preproc_type';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = preproc_default;
inputText(numParameters).flag = 'scalar';
inputText(numParameters).help = struct('title', 'Pre-processing', 'string', char('Data is pre-processed prior to the first data reduction. Options are discussed below:', '1) Remove Mean Per Timepoint - At each time point, image mean is removed.', ...
    '2) Remove Mean Per Voxel - Time-series mean is removed at each voxel', '3) Intensity Normalization - At each voxel, time-series is scaled to have a mean of 100. Since the data is already scaled to percent signal change, there is no need to scale the components.', ...
    '4) Variance Normalization - At each voxel, time-series is linearly detrended and converted to Z-scores.'));

numParameters = numParameters + 1;

inputText(numParameters).promptString = 'Do You Want To Scale The Results?';
inputText(numParameters).answerString = scaleOptions;
inputText(numParameters).uiType = 'popup';
inputText(numParameters).dataType = 'string';
inputText(numParameters).tag = 'scaleType';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = scale_default;
inputText(numParameters).flag = 'scalar';
inputText(numParameters).help = struct('title', 'Scaling', 'string', char('1) Scale to Original Data(%) - Components are scaled to original data units to represent percent signal change.', ...
    '2) Z-scores - Components are scaled to z-scores.', ...
    '3) Scaling in Timecourses - Normalize spatial maps using maximum value (not absolute value) and apply it to timecourses.', ...
    '4) Scaling in Maps and Timecourses - Spatial maps are scaled using the standard deviation of timecourses and timecourses are scaled using the maximum spatial intensity value.'));

numParameters = numParameters + 1;

inputText(numParameters).promptString = 'Enter number of matlab workers. Use value of 1 for serial processing.';
inputText(numParameters).answerString = '4';
inputText(numParameters).uiType = 'edit';
inputText(numParameters).dataType = 'numeric';
inputText(numParameters).tag = 'numWorkers';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = 1;
inputText(numParameters).flag = 'scalar';
inputText(numParameters).help = struct('title', 'Run Parallel?', 'string', char('1. Serial - Group ICA will be run in serial mode. This option is preferred for analyzing smaller data-sets.', ...
    '2. Parallel - Group ICA will be run in parallel mode. Very useful option for analyzing large data. If you don''t have parallel computing toolbox parts of code will be run in sessions. Option is provided to enter number of sessions.'));

answers = icatb_inputDialog('inputtext', inputText, 'Title', dlg_title, 'handle_visibility', 'on');

if (~isempty(answers))
    
    handles.preproc_type = answers{1};
    handles.scaleType = answers{2};
    handles.numWorkers = answers{3};
    
end


guidata(handles.output, handles);


function generateBatch(handles)


% Enter no. of dummy scans to exclude from the group ICA analysis. If you have no dummy scans leave it as 0.
try
    dummy_scans = handles.filesInfo.file_numbers;
catch
    dummy_scans = 0;
end

if (isempty(dummy_scans))
    dummy_scans = 0;
end

preproc_type = 'Remove Mean Per Timepoint';
try
    preproc_type = handles.preproc_type;
catch
end

scaleType = 'z-scores';
try
    scaleType = handles.scaleType;
catch
end

dataFileName = fullfile( handles.outputDir, [handles.prefix, '_files.mat']);

if (~isfield(handles, 'bids_info'))
    files = handles.input_data_file_patterns;
else
    bids_info = handles.bids_info;
    files = {};
end

display_results.formatName = 'html';

network_summary_opts = [];
try
    network_summary_opts = handles.network_summary_opts;
catch
end

display_results.network_summary_opts = network_summary_opts;

numWorkers = 4;
try
    numWorkers = handles.numWorkers;
catch
end

if (length(files) == 1)
    numWorkers = 1;
end

parallel_info.mode = 'serial';

if (numWorkers > 1)
    parallel_info.mode = 'parallel';
end

parallel_info.num_workers = numWorkers;
if (~exist('bids_info', 'var') || isempty(bids_info))
    save(dataFileName, 'files', 'display_results', 'parallel_info');
else
    save(dataFileName, 'bids_info', 'display_results', 'parallel_info');
end


clear display_results files parallel_info;

%% Write batch file with necessary info
batchInfo.modalityType = 'fmri';
comments{1} = '%% Modality Type';
batchInfo.outputDir = handles.outputDir;
comments{end + 1} = '%% Output Directory';
batchInfo.prefix = handles.prefix;
comments{end + 1} = '%% All the output files will be preprended with the specified prefix.';
%batchInfo.parallel_info = 'parallel_info;';
%comments{end + 1} = '%% Enter number of parallel workers.';
batchInfo.perfType = 1;
comments{end + 1} = '%% Group PCA performance settings. Best setting for each option will be selected based on variable MAX_AVAILABLE_RAM in icatb_defaults.m.';
batchInfo.dataSelectionMethod = 4;

if (exist('bids_info', 'var') && ~isempty(bids_info))
    
    comments{end + 1} = '%% Bids Info';
    batchInfo.bids_info = 'bids_info;';
    comments{end + 1} = '%% Specify BIDS info. Data file patterns is read from the bids structure and input_data_file_patterns variable is skipped';
    
else
    
    comments{end + 1} = '%% Data selection option. If option 4 is specified, file names must be entered in input_data_file_patterns';
    batchInfo.input_data_file_patterns = 'files;';
    comments{end + 1} = '%% Input data file pattern for data-sets must be in a cell array. The no. of rows of cell array correspond to no. of data-sets.';
    
end

batchInfo.input_design_matrices = {};
comments{end + 1} = '%% Design matrix/matrices.';
batchInfo.dummy_scans = dummy_scans;
comments{end + 1} = '%% Number of dummy scans.';

batchInfo.maskFile = handles.maskFile;
comments{end + 1} = '%% Full file path of the mask file.';

batchInfo.TR = handles.TR;
comments{end + 1} = '%% TR in seconds';

batchInfo.preproc_type = preproc_type;
comments{end + 1} = '%% Data pre-processing option. By default, Remove Mean Per Timepoint is used';

batchInfo.scaleType = scaleType;
comments{end + 1} = '%% Scale components. By default, components are converted to z-scores.';
batchInfo.algoType = handles.algorithm;
comments{end + 1} = '%% ICA algorithm to use.';
batchInfo.refFiles = handles.refFiles;
comments{end + 1} = '%% Reference templates';
%batchInfo.display_results = 'display_results;';
%comments{end + 1} = '%% Display results';

batchFileName = fullfile(handles.outputDir, [handles.prefix, '_gift_batch_file.m']);

fnames = fieldnames(batchInfo);

%defs = repmat({''}, length(fnames), 1);
defs{1} = sprintf('load (''%s'');', dataFileName);
for nF = 1:length(fnames)
    tmp = batchInfo.(fnames{nF});
    if (isempty(tmp))
        tmp = [fnames{nF}, ' = {};'];
    else
        if (~strcmpi(fnames{nF}, 'input_data_file_patterns') && ~strcmpi(fnames{nF}, 'parallel_info') && ~strcmpi(fnames{nF}, 'display_results') ...
                && ~strcmpi(fnames{nF}, 'bids_info'))
            if (isnumeric(tmp))
                tmp = sprintf('%s = %s;', fnames{nF}, mat2str(tmp, 4));
            else
                tmp = sprintf('%s = ''%s'';', fnames{nF}, tmp);
            end
        else
            
            tmp = [fnames{nF}, ' = ', tmp];
        end
    end
    
    defs{end + 1} = comments{nF};
    defs{end + 1} = tmp;
    defs{end + 1} = '';
    
end


dlmwrite(batchFileName, char(defs), '');
disp(char(['Input information is saved in batch file ', batchFileName, '.'], 'Edit the batch file according to your needs and ', ...
    ['use icatb_batch_file_run(''', batchFileName, '''); to run the batch analysis']));
fprintf('\n\n');

icatb_batch_file_run(batchFileName);


function matchInd = matchString(selectedStr, options, def)

if (ischar(selectedStr))
    matchInd = strmatch(lower(selectedStr), lower(cellstr(options)), 'exact');
else
    matchInd = selectedStr;
end

if (isempty(matchInd))
    matchInd = def;
end


function compLabels = getCompLabels(txtFile)

fid = fopen(txtFile, 'r');
if (fid == -1)
    error(['File ', txtFile, ' cannot be opened for reading']);
end
try
    dd = textscan(fid, '%s', 'delimiter', '\t\n,', 'multipleDelimsAsOne', 1, 'whitespace', ' ');
    val = dd{1};
catch
    val = [];
end
fclose(fid);
val = val(icatb_good_cells(val));
chk = cellfun('isempty', regexp(val, '^\d+$'));

inds = find(chk == 1);

compLabels = repmat(struct('name', '', 'value', []), 1, length(inds));
for nI = 1:length(inds)
    compLabels(nI).name = val{inds(nI)};
    if (nI == length(inds))
        endT = length(val);
    else
        endT = inds(nI + 1) - 1;
    end
    dd = str2num(char(val{inds(nI) + 1:endT}));
    compLabels(nI).value = dd(:)';
end


% --- Executes on key press with focus on template_mask_file and none of its controls.
function template_mask_file_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to template_mask_file (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


function openTRWindow(hObject, event_data, InputHandle)
%% Open TR WIndows


TR = get(InputHandle.TR, 'string');

prompt = {'Enter TR in seconds (One per subject):'};
name = 'TR';
numlines = 20;
defaultanswer = {TR};
answer = inputdlg(prompt,name,numlines,defaultanswer);

try
    set(InputHandle.TR, 'string', answer{1});
catch
end


% --------------------------------------------------------------------
function disp_defaults_Callback(hObject, eventdata, handles)
% hObject    handle to disp_defaults (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


conn_threshold = -Inf;
try
    conn_threshold = handles.network_summary_opts.conn_threshold;
catch
end

display_type = 'slices';
try
    display_type = handles.network_summary_opts.display_type;
catch
end

slice_plane = 'axial';
try
    slice_plane = handles.network_summary_opts.slice_plane;
catch
end

display_threshold = 1.0;
try
    display_threshold = handles.network_summary_opts.threshold;
catch
end

ZOptions = {'No', 'Yes'};
convert_to_z = 'yes';
try
    convert_to_z = handles.network_summary_opts.convert_to_z;
catch
end

displayTypeOpts = {'Slices', 'Render'};

display_type_def = matchString(display_type, displayTypeOpts, 1);

planeOptions = {'Axial', 'Sagittal', 'Coronal'};

plane_default = matchString(slice_plane, planeOptions, 1);

z_default = matchString(convert_to_z, ZOptions, 1);

dlg_title = 'Select options for network summary';

numParameters = 1;

inputText(numParameters).promptString =  'Enter FNC connectivity threshold (connectogram)';
inputText(numParameters).answerString = num2str(conn_threshold);
inputText(numParameters).uiType = 'edit';
inputText(numParameters).dataType = 'numeric';
inputText(numParameters).tag = 'conn_threshold';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = 1;
inputText(numParameters).flag = 'scalar';

numParameters = numParameters + 1;

inputText(numParameters).promptString = 'Select display type';
inputText(numParameters).answerString = displayTypeOpts;
inputText(numParameters).uiType = 'popup';
inputText(numParameters).dataType = 'string';
inputText(numParameters).tag = 'display_type';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = display_type_def;
inputText(numParameters).flag = 'scalar';


numParameters = numParameters + 1;

inputText(numParameters).promptString = 'Select slice plane';
inputText(numParameters).answerString = planeOptions;
inputText(numParameters).uiType = 'popup';
inputText(numParameters).dataType = 'string';
inputText(numParameters).tag = 'slice_plane';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = plane_default;
inputText(numParameters).flag = 'scalar';


numParameters = numParameters + 1;

inputText(numParameters).promptString = 'Convert to z-scores';
inputText(numParameters).answerString = ZOptions;
inputText(numParameters).uiType = 'popup';
inputText(numParameters).dataType = 'string';
inputText(numParameters).tag = 'z_options';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = z_default;
inputText(numParameters).flag = 'scalar';


numParameters = numParameters + 1;

inputText(numParameters).promptString = 'Select display threshold (spatial maps)';
inputText(numParameters).answerString = num2str(display_threshold);
inputText(numParameters).uiType = 'edit';
inputText(numParameters).dataType = 'numeric';
inputText(numParameters).tag = 'threshold';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = 1;
inputText(numParameters).flag = 'scalar';


% display_results.network_summary_opts.threshold = display_results.threshold;
% display_results.network_summary_opts.convert_to_z = display_results.convert_to_zscores;


answers = icatb_inputDialog('inputtext', inputText, 'Title', dlg_title, 'handle_visibility', 'on');

if (~isempty(answers))
    
    handles.network_summary_opts.conn_threshold  = answers{1};
    handles.network_summary_opts.display_type = answers{2};
    handles.network_summary_opts.slice_plane = answers{3};
    handles.network_summary_opts.convert_to_z = answers{4};
    handles.network_summary_opts.threshold = answers{5};
end

guidata(handles.output, handles);
