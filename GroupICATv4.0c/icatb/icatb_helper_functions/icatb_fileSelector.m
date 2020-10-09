function fileInfo = icatb_fileSelector(fileInfo, filePattern)
%% File selector
%
% Outputs:
%
% fileInfo -
%
icatb_defaults;
global UI_FS;

filesList = '';
file_numbers = '';
numOfSess = 1;

try
    filesList = fileInfo.filesList;
catch
end

try
    file_numbers = fileInfo.file_numbers;
catch
end

try
    numOfSess = fileInfo.numOfSess;
catch
end

if (isnumeric(file_numbers))
    file_numbers = num2str(file_numbers);
end

if (~exist('filePattern', 'var'))
    filePattern = '*nii';
end

fileInfo.filesList = filesList;
fileInfo.file_numbers = file_numbers;
fileInfo.filePattern = filePattern;

figTag = 'file_selector_gui';

% Delete figures that have tag
if ~isempty(findobj(0, 'tag', figTag))
    delete(findobj(0, 'tag', figTag));
end

graphicsHandle = icatb_getGraphics('File Selector', 'displaygui', figTag, 'on');
set(graphicsHandle, 'menubar', 'none');
set(graphicsHandle, 'CloseRequestFcn', @figCloseCallback);
set(graphicsHandle, 'userdata', fileInfo);


% Offsets
xOffset = 0.05; yOffset = 0.05;
buttonHeight = 0.052; promptHeight = 0.052;
editTextHeight = 0.05; editTextWidth = 0.7; yPos = 0.92;

editTextPos = [xOffset, yPos - 0.5*yOffset - 0.5*promptHeight, editTextWidth, promptHeight];

% Plot Text
promptH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'text', ...
    'position', editTextPos, 'String', 'Select files or file patterns for all subjects and sessions ...', 'fontsize', UI_FS - 1, ...
    'horizontalalignment', 'center');

promptH = icatb_wrapStaticText(promptH);

editTextPos = get(promptH, 'position');

editTextHeight = 0.38;
editTextPos(2) = editTextPos(2)  - yOffset - editTextHeight;
editTextPos(4) = editTextHeight;

% Plot edit
editTextH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'edit', ...
    'position', editTextPos, 'String', filesList, 'fontsize', UI_FS - 1, ...
    'horizontalalignment', 'left', 'min', 0, 'max', 2, 'tag', 'filesList');

% Browse
buttonPos = [editTextPos(1) + editTextPos(3) + xOffset, editTextPos(2) + 0.5*editTextPos(4) - 0.5*yOffset, 0.12, buttonHeight];
buttonH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'pushbutton', ...
    'position', buttonPos, 'String', 'browse', 'fontsize', UI_FS - 1, 'horizontalalignment', 'center', 'tooltipstring', ...
    'Select files or paste files in the editbox ...', 'callback', {@selectFiles, editTextH});



% Prompt
editTextPos(2) = editTextPos(2) - 1.5*yOffset - promptHeight;
editTextPos(4) = promptHeight;

promptTextPos = editTextPos;

% Plot Text
promptH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'text', ...
    'position', promptTextPos, 'String', 'Enter number of sessions', 'fontsize', UI_FS - 1, ...
    'horizontalalignment', 'center');

promptH = icatb_wrapStaticText(promptH);

editTextPos = get(promptH, 'position');

editTextPos(1) = editTextPos(1) + editTextPos(3) + xOffset;
editTextPos(3) = 0.15;
editTextPos(2) = editTextPos(2) + (0.5*editTextPos(4) - 0.5*promptHeight);
editTextPos(4) = promptHeight;

% Plot edit
editTextH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'edit', ...
    'position', editTextPos, 'String', num2str(numOfSess), 'fontsize', UI_FS - 1, ...
    'horizontalalignment', 'left', 'tag', 'numOfSess');


% Prompt
promptTextPos(2) = editTextPos(2) - 1.5*yOffset - promptHeight;
%editTextPos(4) = promptHeight;

%promptTextPos = editTextPos;

% Plot Text
promptH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'text', ...
    'position', promptTextPos, 'String', 'Enter scans to include. Leave empty to select all scans.', 'fontsize', UI_FS - 1, ...
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
    {@doneCallback, graphicsHandle});

try
    waitfor(graphicsHandle);
catch
end

appName = 'fileAppData';
if (isappdata(0, appName))
    fileInfo = getappdata(0, appName);
    rmappdata(0, appName);
end




function selectFiles(hObject, event_data, handles)
%% Select files
%

parentH = get(hObject, 'parent');
fileInfo = get(parentH, 'userdata');

set(parentH, 'pointer', 'watch');

files = icatb_selectEntry('title', 'Select files ...', 'typeEntity', 'file', 'typeSelection', 'multiple', 'filter', fileInfo.filePattern);

if (~isempty(files))
    set(handles, 'string', files);
end

set(parentH, 'pointer', 'arrow');

function doneCallback(hObject, event_data, handles)
%% Done callback
%

fileInfo = get(handles, 'userdata');

editH = findobj(handles, 'tag', 'filesList');

files = strtrim(deblank(get(editH, 'string')));

if (isempty(files))
    error('Files are not selected');
end

fileNumH = findobj(handles, 'tag', 'file_numbers');
file_numbers = str2num(deblank(get(fileNumH, 'string')));
%file_numbers(file_numbers > fileInfo.numOfDataSets) = [];


sessH = findobj(handles, 'tag', 'numOfSess');
numOfSess = str2num(deblank(get(sessH, 'string')));

fileInfo.file_numbers = file_numbers;
fileInfo.filesList = files;
fileInfo.numOfSess = numOfSess;

set(handles, 'userdata', fileInfo);
setappdata(0, 'fileAppData', fileInfo);
delete(handles);


function figCloseCallback(hObject, event_data, handles)
%% Fig close callback

delete(hObject);
