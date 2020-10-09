function varargout = setup_dyn_fc(varargin)
%SETUP_DYN_FC MATLAB code file for setup_dyn_fc.fig
%      SETUP_DYN_FC, by itself, creates a new SETUP_DYN_FC or raises the existing
%      singleton*.
%
%      H = SETUP_DYN_FC returns the handle to a new SETUP_DYN_FC or the handle to
%      the existing singleton*.
%
%      SETUP_DYN_FC('Property','Value',...) creates a new SETUP_DYN_FC using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to setup_dyn_fc_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      SETUP_DYN_FC('CALLBACK') and SETUP_DYN_FC('CALLBACK',hObject,...) call the
%      local function named CALLBACK in SETUP_DYN_FC.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help setup_dyn_fc

% Last Modified by GUIDE v2.5 18-Sep-2020 18:48:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @setup_dyn_fc_OpeningFcn, ...
    'gui_OutputFcn',  @setup_dyn_fc_OutputFcn, ...
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


% --- Executes just before setup_dyn_fc is made visible.
function setup_dyn_fc_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for setup_dyn_fc
handles.output = hObject;

set(handles.output, 'visible', 'off');

outputDir = icatb_selectEntry('typeEntity', 'directory', 'title', 'Select dFC ROI Analysis Output Directory');

pause(1);

set(handles.output, 'visible', 'on');

handles.dfcRoiInfo.userInput.outputDir = outputDir;

set(handles.atlas_mask_file, 'string', fullfile(fileparts(which('gift.m')), 'toolbox', 'noisecloud', 'mr', 'raw', 'aal2mni152.nii'));

set(handles.atlas_labels_file, 'string', fullfile(fileparts(which('gift.m')), 'toolbox',  'noisecloud', 'mr', 'raw', 'AALLabels.txt'));

% Update handles structure
guidata(hObject, handles);

movegui(hObject, 'center');

% UIWAIT makes setup_dyn_fc wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = setup_dyn_fc_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
%varargout{1} = handles.output;
waitfor(hObject);
dfcRoiInfo = [];
appName = 'roiAppData';
if (isappdata(0, appName))
    dfcRoiInfo = getappdata(0, appName);
    rmappdata(0, appName);
end

varargout{1} = dfcRoiInfo;

function atlas_mask_file_Callback(hObject, eventdata, handles)
% hObject    handle to atlas_mask_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of atlas_mask_file as text
%        str2double(get(hObject,'String')) returns contents of atlas_mask_file as a double


% --- Executes during object creation, after setting all properties.
function atlas_mask_file_CreateFcn(hObject, eventdata, handles)
% hObject    handle to atlas_mask_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in atlas_mask.
function atlas_mask_Callback(hObject, eventdata, handles)
% hObject    handle to atlas_mask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

results = guidata(handles.output);
P = icatb_selectEntry('filter', '*nii;*.img', 'title', 'Select atlas file ...', 'fileType', 'image', ...
    'fileNumbers', 1);
if (~isempty(P))
    set(results.atlas_mask_file, 'string', deblank(P));
end


function atlas_labels_file_Callback(hObject, eventdata, handles)
% hObject    handle to atlas_labels_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of atlas_labels_file as text
%        str2double(get(hObject,'String')) returns contents of atlas_labels_file as a double


% --- Executes during object creation, after setting all properties.
function atlas_labels_file_CreateFcn(hObject, eventdata, handles)
% hObject    handle to atlas_labels_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in atlas_labels.
function atlas_labels_Callback(hObject, eventdata, handles)
% hObject    handle to atlas_labels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

results = guidata(handles.output);
P = icatb_selectEntry('filter', '*.txt', 'title', 'Select atlas labels file ...', 'typeSelection', 'file', 'typeSelection', 'single');
if (~isempty(P))
    set(results.atlas_labels_file, 'string', deblank(P));
end

% --- Executes on button press in select_data.
function select_data_Callback(hObject, eventdata, handles)
% hObject    handle to select_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

filesInfo = [];
try
    filesInfo = handles.dfcRoiInfo.userInput.filesInfo;
catch
end

filesInfo = icatb_fileSelector(filesInfo, '*.nii');
if (~isempty(filesInfo.filesList))
    set(handles.select_data, 'ForegroundColor', [0, 1, 0]);
end

handles.dfcRoiInfo.userInput.filesInfo = filesInfo;
handles.dfcRoiInfo.userInput.numOfDataSets = size(filesInfo.filesList, 1);
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

try
    
    if (getMask == 2)
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

handles.dfcRoiInfo.userInput.maskFile = maskFile;

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

dfcRoiInfo = results.dfcRoiInfo;

outputDir = pwd;

try
    outputDir = dfcRoiInfo.userInput.outputDir;
catch
end

dfcRoiInfo.userInput.outputDir = outputDir;

% Specify data
filesList = '';
try
    filesList = dfcRoiInfo.userInput.filesInfo.filesList;
catch
end

if (isempty(filesList))
    error('data is not specified');
end

dfcRoiInfo.userInput.numOfDataSets = size(filesList, 1);


% Mask
maskFile = [];
try
    maskFile = dfcRoiInfo.userInput.maskFile;
catch
end

dfcRoiInfo.userInput.maskFile = maskFile;

% TR
TR = str2num(get(results.TR, 'string'));
if (isempty(TR))
    error('Specify TR of the experiment');
end

dfcRoiInfo.userInput.TR = TR;

prefix = get(results.prefix, 'string');
if (isempty(prefix))
    error('Specify an output prefix');
end

dfcRoiInfo.userInput.prefix = prefix;

% Analysis type
analysisType = 'roi-roi';
if (~isempty(icatb_findstr(lower(get(get(results.analysis_type_grp, 'SelectedObject'), 'string')), 'roi-voxel')))
    analysisType = 'roi-voxel';
end

dfcRoiInfo.userInput.analysisType = analysisType;


% Atlas info
atlas_mask_file = strtrim(get(results.atlas_mask_file, 'string'));
if (isempty(atlas_mask_file))
    error('Specify atlas nifti file');
end

atlas_labels_file = strtrim(get(results.atlas_labels_file, 'string'));
if (isempty(atlas_labels_file))
    error('Specify atlas labels file');
end

dfcRoiInfo.userInput.atlas_mask_file = atlas_mask_file;
dfcRoiInfo.userInput.atlas_labels_file = atlas_labels_file;

labels = icatb_textscan(atlas_labels_file);
index = icatb_listdlg('PromptString', 'Select labels of interest', 'SelectionMode', 'multiple',...
    'ListString', labels, 'movegui', 'center', 'windowStyle', 'modal');
if (isempty(index))
    error('Labels are not selected');
end

index = sort(index);

if ((length(labels) == 1) && strcmpi(dfcRoiInfo.userInput.analysisType , 'roi-roi'))
    error('You need to select sufficient number of labels to do roi-roi analysis');
end

labels = cellstr(labels);

dfcRoiInfo.userInput.roi_info = {labels, index};

setappdata(0, 'roiAppData', dfcRoiInfo);

delete(gcbf);

drawnow;


% --------------------------------------------------------------------
function dfnc_defaults_Callback(hObject, eventdata, handles)
% hObject    handle to dfnc_defaults (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

defaultsCallback(hObject, [], handles);


function defaultsCallback(hObject, ed, handles, figVisibility)
%% Defaults callback

if (~exist('figVisibility', 'var'))
    figVisibility = 'on';
end

res = guidata(handles.output);
dfncInfo = res.dfcRoiInfo;

try
    covInfo = dfncInfo.userInput.feature_params.final.tc_covariates;
catch
end

try
    covInfo.numOfDataSets = dfncInfo.userInput.numOfDataSets;
catch
    error('data needs to be selected before accessing defaults');
end

default_params = icatb_dfnc_options('covInfo', covInfo);

if (~isfield(dfncInfo.userInput, 'feature_params'))
    feature_params = default_params;
else
    feature_params = dfncInfo.userInput.feature_params;
end

tags = cellstr(char(feature_params(1).inputParameters(1).options.tag));
motion_inds = strmatch('tc_covariates', tags, 'exact');
if (isempty(motion_inds))
    feature_params(1).inputParameters(1).options(end+1) = feature_params(1).inputParameters(1).options(end);
    feature_params(1).inputParameters(1).options(end).callback = [];
    feature_params(1).inputParameters(1).options(end) = default_params(1).inputParameters(1).options(end);
end


for nF = 1:length(feature_params.inputParameters)
    for nO = 1:length(feature_params.inputParameters(nF).options)
        if (~isfield(feature_params.inputParameters(nF).options(nO), 'enable') || isempty(feature_params.inputParameters(nF).options(nO).enable))
            feature_params.inputParameters(nF).options(nO).enable = 'on';
            feature_params.defaults(nF).options(nO).enable = 'on';
        end
    end
end

tags = cellstr(char(feature_params(1).inputParameters(2).options.tag));
check_window_type = strmatch('window_type', tags, 'exact');
check_window_alpha = strmatch('window_alpha', tags, 'exact');
if (~isempty(check_window_type))
    if (feature_params(1).inputParameters(2).options(check_window_type).value == 1)
        feature_params(1).inputParameters(2).options(check_window_alpha).answerString = '3';
    end
    feature_params(1).inputParameters(2).options(check_window_type) = [];
end

out = icatb_OptionsWindow(feature_params.inputParameters, feature_params.defaults, figVisibility, ...
    'title', 'dFC ROI Options');

drawnow;

feature_params.inputParameters = out.inputParameters;
feature_params.final = out.results;
dfncInfo.userInput.feature_params = feature_params;


res.dfcRoiInfo = dfncInfo;
guidata(handles.output, res);



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
