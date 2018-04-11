function [parameters, inputText] = icatb_displayGUI_defaults(cmd_str, inputText, parameters, handle_visibility)
%
% Purpose: show display GUI defaults
%
% Input:
% 1. cmd_str - Command String
% 2. inputText - inputText structure for plotting
% 3. parameters - display parameters
%
% Output:
% parameters - Ouput variable
%

if ~exist('cmd_str', 'var')
    cmd_str = 'init';
end

if ~exist('handle_visibility', 'var')
    handle_visibility = 'on';
end

structuralFile = parameters.structFile;

% run defaults file to get global variables
icatb_defaults;
global PARAMETER_INFO_MAT_FILE; %Holds information for group session parameters
global FUNCTIONAL_DATA_FILTER;

% Screen Color Defaults
global BG_COLOR;
global BG2_COLOR;
global BUTTON_COLOR;
global BUTTON_FONT_COLOR;
global FONT_COLOR;
global AXES_COLOR;
global FONT_COLOR;

% Fonts
global UI_FONTNAME;
global UI_FONTUNITS;
global UI_FS;

% display_defaults
global SORT_COMPONENTS;
global IMAGE_VALUES;
global THRESHOLD_VALUE;
global CONVERT_Z;
global IMAGES_PER_FIGURE;
global ANATOMICAL_PLANE;
global COMPLEX_IMAGES_PER_FIGURE;

% get the data type
if isfield(parameters, 'dataType')
    dataType = parameters.dataType;
else
    dataType = 'real';
end

if strcmpi(cmd_str, 'init')
    % get the input text using the icatb_defaults.m file
    %--Image Values
    options(1).str = 'Positive and Negative';
    options(2).str = 'Positive';
    options(3).str = 'Absolute Value';
    options(4).str = 'Negative'; % Added negative image values

    % Apply defaults from icatb_defaults
    imageOptions = {'Positive and Negative', 'Positive', 'Absolute Value', 'Negative'};
    options(1).str =  IMAGE_VALUES;

    % Apply defaults from icatb_defaults
    imageOptions = checkDefaults_gui(options(1).str, imageOptions, 'exact');

    for jj = 1:length(imageOptions)
        options(jj).str = imageOptions{jj};
    end

    % needed for all the visualization methods
    numParameters = 1;
    inputText(numParameters).promptString = 'Image Values';
    inputText(numParameters).answerString = str2mat(options.str);
    inputText(numParameters).uiType = 'popup';
    inputText(numParameters).dataType = 'string';
    inputText(numParameters).tag = 'image_values';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).uiPos = [0.4 0.05];
    inputText(numParameters).value = 1;
    clear options;


    %--Convert To Z Scores
    % Apply defaults from icatb_defaults
    zOptions = {'No', 'Yes'};
    options(1).str =  CONVERT_Z;
    % Apply defaults from icatb_defaults
    zOptions = checkDefaults_gui(options(1).str, zOptions, 'exact');
    for jj = 1:length(zOptions)
        options(jj).str = zOptions{jj};
    end

    numParameters = numParameters + 1;
    % used for all the visualization methods
    inputText(numParameters).promptString = 'Convert To Z Scores';
    inputText(numParameters).answerString = str2mat(options.str);
    inputText(numParameters).uiType = 'popup';
    inputText(numParameters).dataType = 'string';
    inputText(numParameters).tag = 'Convert_To_ZScores';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).uiPos = [0.25 0.05];
    inputText(numParameters).value = 1;
    clear options;

    numParameters = numParameters + 1;
    inputText(numParameters).promptString = 'Threshold Value';
    inputText(numParameters).answerString = THRESHOLD_VALUE;
    inputText(numParameters).uiType = 'edit';
    inputText(numParameters).dataType = 'numeric';
    inputText(numParameters).tag = 'Threshold_Value';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).uiPos = [0.25 0.05];
    inputText(numParameters).value = 1;

    if strcmp(dataType, 'real')
        % options for number of images per figure
        otherOptions = {'1', '4', '9', '16', '25'};
        options(1).str =  IMAGES_PER_FIGURE;
    else
        % handle complex data
        % options for number of images per figure
        otherOptions = {'1', '2', '8'};
        options(1).str =  COMPLEX_IMAGES_PER_FIGURE;
    end
    % Apply defaults from icatb_defaults
    otherOptions = checkDefaults_gui(options(1).str, otherOptions, 'exact');

    for jj = 1:length(otherOptions)
        options(jj).str = otherOptions{jj};
    end
    numParameters = numParameters + 1;
    inputText(numParameters).promptString = 'Images Per Figure';
    inputText(numParameters).answerString = str2mat(options.str);
    inputText(numParameters).uiType = 'popup';
    inputText(numParameters).dataType = 'numeric';
    inputText(numParameters).tag = 'images_per_figure';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).uiPos = [0.25 0.05];
    inputText(numParameters).value = 1;
    clear options;

    planeOptions = {'Axial', 'Sagittal', 'Coronal'};
    planeOptions = checkDefaults_gui(lower(ANATOMICAL_PLANE), lower(planeOptions), 'optional');

    for jj = 1:length(planeOptions)
        options(jj).str = planeOptions{jj};
    end

    numParameters = numParameters + 1;
    inputText(numParameters).promptString = 'Anatomical Plane';
    inputText(numParameters).answerString = str2mat(options.str);
    inputText(numParameters).uiType = 'popup';
    inputText(numParameters).dataType = 'string';
    inputText(numParameters).tag = 'anatomical_plane';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).uiPos = [0.4 0.05];
    inputText(numParameters).value = 1;

    clear options;

    slicePlane = lower(ANATOMICAL_PLANE);
    imagVol = icatb_get_vol_nifti(structuralFile);
    % get the slices in mm for the corresponding plane
    [sliceParameters] = icatb_get_slice_def(imagVol, slicePlane);
    % get the slices in mm
    slices_in_mm = sliceParameters.slices;
    clear sliceParameters;
    % construct string
    slices_in_mm = icatb_constructString(slices_in_mm);
    options(1).str = slices_in_mm;

    numParameters = numParameters + 1;
    inputText(numParameters).promptString = 'Slice Range';
    inputText(numParameters).answerString = options.str;
    inputText(numParameters).uiType = 'edit';
    inputText(numParameters).dataType = 'numeric';
    inputText(numParameters).tag = 'slice_range';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).uiPos = [0.4 0.05];
    inputText(numParameters).value = 1;
    clear options;
    %%%%%%%% done with defining the parameter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
% end for getting the defaults

figureData.inputText = inputText;
figureData.structFile = structuralFile;
figureData.parameters = parameters;

% Setup figure for GUI
[InputHandle] = icatb_getGraphics('Display Options', 'normal', 'figure');

set(InputHandle, 'windowstyle', 'modal');

% Set no menu bar for the figure
set(InputHandle, 'Menubar', 'None', 'tag', 'disp_parameters_gui', 'userdata', figureData);

[InputHandle] = icatb_plot_controls_fig(inputText, get(InputHandle, 'tag'), handle_visibility, 'OK', 'Cancel', ...
    cellstr(str2mat(inputText.tag)), InputHandle);

cancelHandle = findobj(InputHandle, 'tag', 'Cancel');
okHandle = findobj(InputHandle, 'tag', 'OK');

%%%%% set function callbacks %%%%%

% Cancel callback
set(cancelHandle, 'callback', {@closeCallback, InputHandle});

% Done callback
set(okHandle, 'callback', {@applyCallback, InputHandle});

% anatomical plane Callback
set(findobj(InputHandle, 'tag', 'anatomical_plane'), 'callback', ...
    {@anatomicalCallback, InputHandle});

if strcmpi(handle_visibility, 'on')

    try
        set(InputHandle, 'visible', 'on');
        waitfor(InputHandle);
    catch
        if ishandle(InputHandle)
            delete(InputHandle);
        end
    end

else
    % call apply callback
    applyCallback(okHandle, [], InputHandle);
end

% get the parameters from the figure
if isappdata(0, 'inputParaDisplay')
    temp = getappdata(0, 'inputParaDisplay');
    inputText = temp.inputText;
    parameters = temp.parameters;
    % remove the application data
    rmappdata(0, 'inputParaDisplay');
end


function optionsString = checkDefaults_gui(choiceString, optionsString, stringChar)
% put the defaults at the Top

temp = optionsString; % Assign options string to a temporary variable

if strcmp(stringChar, 'exact')
    matchIndex = strmatch(lower(choiceString), lower(temp), stringChar); % find the index of the choice string
else
    matchIndex = strmatch(lower(choiceString), lower(temp)); % find the index of the choice string
end

jj = 1; optionsString{1} = temp{matchIndex}; % choice string

for numOptions = 1:length(optionsString)
    if numOptions ~= matchIndex
        jj = jj + 1;
        optionsString{jj} = temp{numOptions}; % collect other options below the choice string
    end
end

function anatomicalCallback(hObject, evd, handles)
% anatomical plane callback

% purpose: updates the editbox slice plane
getString = lower(get(hObject, 'string'));
getValue = get(hObject, 'value');
slicePlane = deblank(getString(getValue, :));

figureData = get(handles, 'userdata');
% get the name of the structural file
structFile = figureData.structFile;
imagVol = icatb_get_vol_nifti(structFile);
% get the slices in mm for the corresponding plane
[sliceParameters] = icatb_get_slice_def(imagVol, slicePlane);
% get the slices in mm
slices_in_mm = sliceParameters.slices;
clear sliceParameters;
% construct string
slices_in_mm = icatb_constructString(slices_in_mm);

% set the editbox to the display the default size of slice plane
set(findobj(handles, 'tag', 'slice_range'), 'string', slices_in_mm);

function applyCallback(handleObj, event_data, handles)
% get the required parameters

figureData = get(handles, 'userdata');
% get the input parameters
inputText = figureData.inputText;

parameters = figureData.parameters; %struct;
% loop over number of input args
for ii = 1:length(inputText)
    dataType = inputText(ii).dataType;
    answerTag = [inputText(ii).tag];
    objHandle = findobj(handles, 'tag', answerTag);
    if strcmp(lower(inputText(ii).uiType), 'edit')
        % check the answer type
        if strcmp(lower(dataType), 'numeric')
            answerVal = str2num(get(objHandle, 'string'));
            inputText(ii).answerString = get(objHandle, 'string');
        else
            answerVal = get(objHandle, 'string');
            inputText(ii).answerString = answerVal;
        end
    elseif strcmp(lower(inputText(ii).uiType), 'popup')
        % check the answer type
        getString = get(objHandle, 'string');
        getVal = get(objHandle, 'value');
        inputText(ii).value = getVal;
        answerVal = deblank(getString(getVal, :));
        % convert to numeric
        if strcmp(lower(dataType), 'numeric')
            answerVal = str2num(answerVal);
        end
    end
    % set field to the parameters
    parameters = setfield(parameters, strrep(deblank(lower(inputText(ii).tag)), '_', ''), answerVal);
end

% get the return Value
if strcmp(lower(parameters.imagevalues), 'positive')
    returnValue = 2;
elseif strcmp(lower(parameters.imagevalues), 'absolute value')
    returnValue = 3;
elseif strcmp(lower(parameters.imagevalues), 'negative')
    returnValue = 4;
else
    returnValue = 1;
end

if strcmp(lower(parameters.converttozscores), 'yes')
    convertToZ = 1;
else
    convertToZ = 0;
end

% pass the return values and convert to z options
parameters.returnValue = returnValue;
parameters.convertToZ = convertToZ;

figureData.parameters = parameters; figureData.inputText = inputText;
% set application data
setappdata(0, 'inputParaDisplay', figureData);

delete(handles);

% close callback
function closeCallback(hObject, event_data, handles)

delete(handles);