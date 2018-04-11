function input_parameters = icatb_mancovan_feature_options(varargin)
%% Feature selection options
%

icatb_defaults;
global DETRENDNUMBER;
global MANCOVA_DEFAULTS;


fnc_lag = 3;
fnc_shift_resolution = 25;

try
    fnc_lag = MANCOVA_DEFAULTS.fnc.lag;
catch
end

try
    fnc_shift_resolution = MANCOVA_DEFAULTS.fnc.shift_resolution;
catch
end

TR = 1;
dims = [];
covInfo = [];

for i = 1:2:length(varargin)
    if (strcmpi(varargin{i}, 'tr'))
        TR = varargin{i + 1};
    elseif (strcmpi(varargin{i}, 'mask_dims'))
        dims = varargin{i + 1};
    elseif (strcmpi(varargin{i}, 'covInfo'))
        covInfo = varargin{i + 1};
    end
end

TR = min(TR);

% List includes Spatial Maps, Timecourses spectra,

%%%%%%%%% Input Parameters Structure
numParameters = 1;
inputParameters(numParameters).listString = 'Spatial Maps';
optionNumber = 1;
% Option 1 of parameter 1
options(optionNumber).promptString = 'Type Of Mask?';
options(optionNumber).answerString = char('Default', 'User specified');
options(optionNumber).uiType = 'popup'; options(optionNumber).value = 1;
options(optionNumber).tag = 'sm_mask'; options(optionNumber).answerType = 'string';
options(optionNumber).flag = 'delete';
options(optionNumber).callback = {@selectMaskCallback, struct('type', 'mask', 'dims', dims)};
options(optionNumber).uiPos = [];
options(optionNumber).enable = [];

optionNumber = optionNumber + 1;
% Option 2 of parameter 1
options(optionNumber).promptString = 'Center spatial maps?';
options(optionNumber).answerString = char('Yes', 'No');
options(optionNumber).uiType = 'popup'; options(optionNumber).value = 1;
options(optionNumber).tag = 'sm_center'; options(optionNumber).answerType = 'string';
options(optionNumber).flag = 'delete';
options(optionNumber).uiPos = [0.12 0.045];
inputParameters(numParameters).options = options; % will be used in plotting the controls in a frame

optionNumber = optionNumber + 1;
% Option 3 of parameter 1
options(optionNumber).promptString = 'Statistic for thresholding';
options(optionNumber).answerString = char('T', 'Z');
options(optionNumber).uiType = 'popup'; options(optionNumber).value = 1;
options(optionNumber).tag = 'stat_threshold_maps'; options(optionNumber).answerType = 'string';
options(optionNumber).flag = 'delete';
options(optionNumber).uiPos = [0.12 0.045];
options(optionNumber).callback = @selectStatCallback;
inputParameters(numParameters).options = options; % will be used in plotting the controls in a frame

optionNumber = optionNumber + 1;
% Option 4 of parameter 1
options(optionNumber).promptString = 'Select Z Threshold';
options(optionNumber).answerString = '1.0';
options(optionNumber).uiType = 'edit'; options(optionNumber).value = 1;
options(optionNumber).tag = 'z_threshold_maps'; options(optionNumber).answerType = 'numeric';
options(optionNumber).flag = 'delete';
options(optionNumber).uiPos = [0.12 0.045];
options(optionNumber).enable = 'inactive';
inputParameters(numParameters).options = options; % will be used in plotting the controls in a frame
clear options;

numParameters = numParameters + 1; % second parameter
inputParameters(numParameters).listString = 'Timecourses Spectra';
optionNumber = 1;
% Option 1 of parameter 2
options(optionNumber).promptString = 'TC Detrend number';
options(optionNumber).answerString = str2mat('0', '1', '2', '3');
matchedIndex = strmatch(num2str(DETRENDNUMBER), options(optionNumber).answerString, 'exact');
options(optionNumber).uiType = 'popup'; options(optionNumber).value = matchedIndex;
options(optionNumber).tag = 'spectra_detrend'; options(optionNumber).answerType = 'numeric';
options(optionNumber).flag = 'delete';
options(optionNumber).uiPos = [0.12 0.045];

optionNumber = optionNumber + 1;
% Option 2 of parameter 2
options(optionNumber).promptString = 'Tapers from dpss';
options(optionNumber).answerString = '3, 5';
options(optionNumber).uiType = 'edit'; options(optionNumber).value = 0;
options(optionNumber).tag = 'spectra_tapers'; options(optionNumber).answerType = 'numeric';
options(optionNumber).flag = 'delete';
options(optionNumber).uiPos = [0.16 0.045];

optionNumber = optionNumber + 1;
% Option 2 of parameter 2
options(optionNumber).promptString = 'Sampling Frequency';
options(optionNumber).answerString = num2str(1/TR);
options(optionNumber).uiType = 'edit'; options(optionNumber).value = 1;
options(optionNumber).tag = 'spectra_sampling_freq'; options(optionNumber).answerType = 'numeric';
options(optionNumber).flag = 'delete';
options(optionNumber).uiPos = [0.16 0.045];

optionNumber = optionNumber + 1;
% Option 3 of parameter 2
options(optionNumber).promptString = 'Frequency band';
options(optionNumber).answerString = num2str([0.0, 1/(2*TR)]);
options(optionNumber).uiType = 'edit'; options(optionNumber).value = 0;
options(optionNumber).tag = 'spectra_freq_band'; options(optionNumber).answerType = 'numeric';
options(optionNumber).flag = 'delete';
options(optionNumber).uiPos = [0.16 0.045];

optionNumber = optionNumber + 1;

options(optionNumber).promptString = 'Use Fractional Amplitude?';
options(optionNumber).answerString = char('Yes', 'No');
options(optionNumber).uiType = 'popup'; options(optionNumber).value = 1;
options(optionNumber).tag = 'spectra_normalize_subs'; options(optionNumber).answerType = 'string';
options(optionNumber).flag = 'delete';
options(optionNumber).uiPos = [0.12 0.045];

optionNumber = optionNumber + 1;

options(optionNumber).promptString = 'Log transform spectra?';
options(optionNumber).answerString = char('Yes', 'No');
options(optionNumber).uiType = 'popup'; options(optionNumber).value = 1;
options(optionNumber).tag = 'spectra_transform'; options(optionNumber).answerType = 'string';
options(optionNumber).flag = 'delete';
options(optionNumber).uiPos = [0.12 0.045];


inputParameters(numParameters).options = options; % will be used in plotting the controls in a frame
clear options;

numParameters = numParameters + 1; % Third parameter
inputParameters(numParameters).listString = 'FNC Correlations';
optionNumber = 1;
% Option 1 of parameter 3
options(optionNumber).promptString = 'Detrend number';
options(optionNumber).answerString = str2mat('0', '1', '2', '3');
matchedIndex = strmatch(num2str(DETRENDNUMBER), options(optionNumber).answerString, 'exact');
options(optionNumber).uiType = 'popup'; options(optionNumber).value = matchedIndex;
options(optionNumber).tag = 'fnc_tc_detrend'; options(optionNumber).answerType = 'numeric';
options(optionNumber).flag = 'delete';
options(optionNumber).uiPos = [0.12 0.045];

optionNumber = optionNumber + 1;
% Option 2 of parameter 3
options(optionNumber).promptString = 'Despike Timecourses?';
options(optionNumber).answerString = char('Yes', 'No');
options(optionNumber).uiType = 'popup'; options(optionNumber).value = 1;
options(optionNumber).tag = 'fnc_tc_despike'; options(optionNumber).answerType = 'string';
options(optionNumber).flag = 'delete';
options(optionNumber).uiPos = [0.12 0.045];

optionNumber = optionNumber + 1;
% Option 2 of parameter 3
options(optionNumber).promptString = 'Filter cutoff (Hz)';
options(optionNumber).answerString = '0.15';
options(optionNumber).uiType = 'edit'; options(optionNumber).value = 1;
options(optionNumber).tag = 'fnc_tc_filter'; options(optionNumber).answerType = 'numeric';
options(optionNumber).flag = 'delete';
options(optionNumber).uiPos = [0.12 0.045];

optionNumber = optionNumber + 1;
% Option 3 of parameter 3
options(optionNumber).promptString = 'Regress covariates?';
options(optionNumber).answerString =  char('None', 'Select');
options(optionNumber).uiType = 'popup'; options(optionNumber).value = 1;
options(optionNumber).tag = 'fnc_tc_covariates'; options(optionNumber).answerType = 'string';
options(optionNumber).flag = 'delete';
options(optionNumber).uiPos = [0.16 0.045];
options(optionNumber).enable = [];
options(optionNumber).callback = {@selectParams, covInfo};

optionNumber = optionNumber + 1;
% Option 4 of parameter 3
options(optionNumber).promptString = 'Min/Max Lag in seconds?';
options(optionNumber).answerString =  num2str(fnc_lag);
options(optionNumber).uiType = 'edit'; options(optionNumber).value = 1;
options(optionNumber).tag = 'fnc_tc_lag'; options(optionNumber).answerType = 'numeric';
options(optionNumber).flag = 'delete';
options(optionNumber).uiPos = [0.16 0.045];
options(optionNumber).enable = [];

optionNumber = optionNumber + 1;
% Option 5 of parameter 3
options(optionNumber).promptString = 'Shift resolution (lag) in ms';
options(optionNumber).answerString =  num2str(fnc_shift_resolution);
options(optionNumber).uiType = 'edit'; options(optionNumber).value = 1;
options(optionNumber).tag = 'fnc_tc_shift_resolution'; options(optionNumber).answerType = 'numeric';
options(optionNumber).flag = 'delete';
options(optionNumber).uiPos = [0.16 0.045];
options(optionNumber).enable = [];



inputParameters(numParameters).options = options; % will be used in plotting the controls in a frame
clear options;

% store the input parameters in two fields
input_parameters.inputParameters = inputParameters;
input_parameters.defaults = inputParameters;

function selectMaskCallback(hObject, event_data, info)
% Select mask callback

gVal = get(hObject, 'value');
gStr = cellstr(get(hObject, 'string'));

hh = gcbf;
ho = findobj(hh, 'tag', ['answer', 'stat_threshold_maps']);
threshH = findobj(hh, 'tag', ['answer', 'z_threshold_maps']);

selectedval = get(ho, 'value');
strs = cellstr(get(ho, 'string'));
selectedStr = deblank(strs{selectedval});

if (strcmpi(gStr{gVal}, 'default'))
    set(threshH, 'enable', 'on');
    if (strcmpi(selectedStr, 't'))
        set(threshH, 'enable', 'inactive');
    end
    set(ho, 'enable', 'on');
    return;
end

set(threshH, 'enable', 'inactive');
set(ho, 'enable', 'off');

dims = info.dims;
mask = icatb_selectEntry('typeEntity', 'file', 'typeSelection', 'multiple', 'filter', '*.img;*.nii', 'filetype', 'image', 'title', ...
    ['Select mask of dimensions (', num2str(dims), ')']);
drawnow;

if (isempty(mask))
    set(hObject, 'value', 1);
    set(threshH, 'enable', 'on');
    if (strcmpi(selectedStr, 't'))
        set(threshH, 'enable', 'inactive');
    end
    set(ho, 'enable', 'on');
    return;
end

data = icatb_loadData(mask); % Mask data
size_data = size(data);
if length(size_data) == 2
    size_data(3) = 1;
end

%% Check the dimensions of the mask w.r.t data
if length(find(size_data == dims)) ~= length(dims)
    msg = sprintf('Mask dimensions ([%s]) doesn''t match that of data dimensions ([%s])', num2str(size_data), num2str(dims));
    disp(msg);
    disp('');
    icatb_errorDialog(msg, 'Mask Error', 'Modal');
    set(hObject, 'value', 1);
    set(threshH, 'enable', 'on');
    set(ho, 'enable', 'on');
    return;
end

set(hObject, 'userdata', deblank(mask));

function selectStatCallback(hObject, event_data, handles)
% Select Stat
%

hh = gcbf;
ho = findobj(hh, 'tag', ['answer', 'stat_threshold_maps']);
val = get(ho, 'value');
strs = cellstr(get(ho, 'string'));
selectedStr = deblank(strs{val});
threshH = findobj(hh, 'tag', ['answer', 'z_threshold_maps']);
set(threshH, 'enable', 'inactive');
if (strcmpi(selectedStr, 'z'))
    set(threshH, 'enable', 'on');
end

function selectParams(hObject, event_data, covInfo)
%% Select covariates of interest
%


if (get(hObject, 'value') == 2)
    
    if (~isempty(get(hObject, 'userdata')))
        covInfo = get(hObject, 'userdata');
    end
    
    icatb_defaults;
    global UI_FS;
    
    filesList = '';
    file_numbers = '';
    numOfDataSets = covInfo.numOfDataSets;
    
    try
        filesList = covInfo.filesList;
    catch
    end
    
    try
        file_numbers = covInfo.file_numbers;
    catch
    end
    
    if (isnumeric(file_numbers))
        file_numbers = num2str(file_numbers);
    end
    
    figTag = 'select_motion_params';
    
    % Delete figures that have tag
    if ~isempty(findobj(0, 'tag', figTag))
        delete(findobj(0, 'tag', figTag));
    end
    
    graphicsHandle = icatb_getGraphics('Select motion parameters', 'normal', figTag, 'on');
    set(graphicsHandle, 'menubar', 'none');
    set(graphicsHandle, 'CloseRequestFcn', {@figCloseCallback, hObject});
    
    % if ispc
    %     set(graphicsHandle, 'windowstyle', 'modal');
    % end
    
    % Offsets
    xOffset = 0.05; yOffset = 0.05;
    buttonHeight = 0.052; promptHeight = 0.052;
    editTextHeight = 0.05; editTextWidth = 0.7; yPos = 0.92;
    
    editTextPos = [xOffset, yPos - 0.5*yOffset - 0.5*promptHeight, editTextWidth, promptHeight];
    
    % Plot Text
    promptH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'text', ...
        'position', editTextPos, 'String', 'Select covariates for all subjects and sessions ...', 'fontsize', UI_FS - 1, ...
        'horizontalalignment', 'center');
    
    promptH = icatb_wrapStaticText(promptH);
    
    editTextPos = get(promptH, 'position');
    
    editTextHeight = 0.52;
    editTextPos(2) = editTextPos(2)  - yOffset - editTextHeight;
    editTextPos(4) = editTextHeight;
    
    % Plot edit
    editTextH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'edit', ...
        'position', editTextPos, 'String', filesList, 'fontsize', UI_FS - 1, ...
        'horizontalalignment', 'left', 'min', 0, 'max', 2, 'tag', 'regress_files');
    
    % Browse
    buttonPos = [editTextPos(1) + editTextPos(3) + xOffset, editTextPos(2) + 0.5*editTextPos(4) - 0.5*yOffset, 0.12, buttonHeight];
    buttonH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'pushbutton', ...
        'position', buttonPos, 'String', 'browse', 'fontsize', UI_FS - 1, 'horizontalalignment', 'center', 'tooltipstring', ...
        'Select files or paste files in the editbox ...', 'callback', {@selectFiles, editTextH});
    
    % Prompt
    editTextPos(2) = editTextPos(2) - 1.5*yOffset - promptHeight;
    editTextPos(4) = promptHeight;
    
    % Plot Text
    promptH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'text', ...
        'position', editTextPos, 'String', 'Enter scans to include. Leave empty to select all scans.', 'fontsize', UI_FS - 1, ...
        'horizontalalignment', 'center');
    
    promptH = icatb_wrapStaticText(promptH);
    
    editTextPos = get(promptH, 'position');
    
    editTextPos(1) = editTextPos(1) + editTextPos(3) + xOffset;
    editTextPos(3) = 0.15;
    editTextPos(2) = editTextPos(2) + (0.5*editTextPos(4) - 0.5*promptHeight);
    editTextPos(4) = promptHeight;
    
    % Plot edit
    editTextH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'edit', ...
        'position', editTextPos, 'String', file_numbers, 'fontsize', UI_FS - 1, ...
        'horizontalalignment', 'left', 'tag', 'file_numbers');
    
    % Plot done
    buttonWidth = 0.12;
    editTextPos(4) = buttonHeight;
    editTextPos(1) = 0.5 - 0.5*buttonWidth;
    editTextPos(2) = yOffset + 0.01;
    editTextPos(3) = buttonWidth;
    
    % Plot done
    icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'pushbutton', ...
        'position', editTextPos, 'String', 'Done', 'fontsize', UI_FS - 1, 'horizontalalignment', 'center', 'callback', ...
        {@doneCallback, graphicsHandle, covInfo, hObject});
    
end


function doneCallback(hObject, event_data, handles, covInfo, dropdownH)
%% Done callback
%

editH = findobj(handles, 'tag', 'regress_files');

files = strtrim(deblank(get(editH, 'string')));

if (isempty(files))
    error('Files are not selected');
end

if (size(files, 1) ~= covInfo.numOfDataSets)
    error(['Number of files doesn''t match the no. of data-sets (', num2str(covInfo.numOfDataSets), ')']);
end

fileNumH = findobj(handles, 'tag', 'file_numbers');
file_numbers = str2num(deblank(get(fileNumH, 'string')));
%file_numbers(file_numbers > covInfo.numOfDataSets) = [];

covInfo.file_numbers = file_numbers;
covInfo.filesList = files;

set(dropdownH, 'userdata', covInfo);

delete(handles);

function selectFiles(hObject, event_data, handles)
%% Select files
%

parentH = get(hObject, 'parent');

set(parentH, 'pointer', 'watch');

files = icatb_selectEntry('title', 'Select files ...', 'typeEntity', 'file', 'typeSelection', 'multiple', 'filter', '*.txt;*.dat');

if (~isempty(files))
    set(handles, 'string', files);
end

set(parentH, 'pointer', 'arrow');

function figCloseCallback(hObject, event_data, dropdownH)
%% Fig close callback

set(dropdownH, 'value', 1);
delete(hObject);