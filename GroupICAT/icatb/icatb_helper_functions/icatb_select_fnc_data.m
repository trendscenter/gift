function sesInfo = icatb_select_fnc_data(sesInfo)
%% Select FNC data
%

icatb_defaults;
global UI_FS;

sessionNames = '';

% if ~exist('sesInfo', 'var')
%     sesInfo.userInput.dataInfo = [];
% end

if (~isfield(sesInfo.userInput, 'dataInfo'))
    sesInfo.userInput.dataInfo = [];
end

try
    sessionNames = (cellstr(char(sesInfo.userInput.dataInfo.name)));
catch
end

fnc_variable_name = 'fnc_corrs_all';
try
    fnc_variable_name = sesInfo.userInput.fnc_variable_mat_file;
catch
end

contrast_vector = '';
try
    contrast_vector = sesInfo.userInput.contrast_vector;
catch
end

if isnumeric(contrast_vector)
    contrast_vector = num2str(contrast_vector);
end

InputHandle = icatb_getGraphics('Select FNC data', 'normal', 'select_fnc_data', 'on');
set(InputHandle, 'menubar', 'none');
set(InputHandle, 'userdata', sesInfo);

controlWidth = 0.52;
promptHeight = 0.05;
promptWidth = controlWidth;
listboxHeight = 0.4; listboxWidth = controlWidth;
xOffset = 0.02; yOffset = promptHeight; yPos = 0.95;
okWidth = 0.12; okHeight = promptHeight;

%% Features text and listbox
promptPos = [0.5 - 0.5*controlWidth, yPos - 0.5*yOffset, promptWidth, promptHeight];

textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Enter FNC data ...', 'tag', ...
    'prompt_fnc', 'fontsize', UI_FS - 1);
icatb_wrapStaticText(textH);
listboxXOrigin = promptPos(1);
listboxYOrigin = promptPos(2) - yOffset - listboxHeight;
listboxPos = [listboxXOrigin, listboxYOrigin, listboxWidth, listboxHeight];
icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'listbox', 'position', listboxPos, 'string', sessionNames, 'tag', ...
    'sessions', 'fontsize', UI_FS - 1, 'min', 0, 'max', 1, 'value', 1, 'callback', {@addSessions, InputHandle});

addButtonPos = [listboxPos(1) + listboxPos(3) + xOffset, listboxPos(2) + 0.5*listboxPos(4) + 0.5*promptHeight, promptHeight + 0.01, promptHeight - 0.01];
removeButtonPos = [listboxPos(1) + listboxPos(3) + xOffset, listboxPos(2) + 0.5*listboxPos(4) - 0.5*promptHeight, promptHeight + 0.01, promptHeight - 0.01];
icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', addButtonPos, 'string', '+', 'tag', 'add_sessions', 'fontsize',...
    UI_FS - 1, 'callback', {@addSessions, InputHandle});
icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', removeButtonPos, 'string', '-', 'tag', 'remove_sessions', 'fontsize',...
    UI_FS - 1, 'callback', {@removeSessions, InputHandle});

promptPos = listboxPos;
promptPos(1) = 0.05;
promptPos(2)= promptPos(2) - 2*yOffset - promptHeight;
promptPos(4) = promptHeight;

textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Enter FNC variable name (Leave empty if parameter file is selected)', 'tag', ...
    'prompt_fnc_var', 'fontsize', UI_FS - 1);
icatb_wrapStaticText(textH);

editPos = promptPos;
editPos(1) = editPos(1) + editPos(3) + xOffset;
editPos(3) = 0.25;
icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'edit', 'position', editPos, 'string', fnc_variable_name, 'tag', ...
    'fnc_variable_name', 'fontsize', UI_FS - 1);

promptPos(2) = promptPos(2) - 2*yOffset - promptHeight;
textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', ...
    'Enter contrast vector for sessions (Leave empty for average)', 'tag', ...
    'prompt_fnc_sess', 'fontsize', UI_FS - 1);
icatb_wrapStaticText(textH);

editPos = promptPos;
editPos(1) = editPos(1) + editPos(3) + xOffset;
editPos(3) = 0.25;
icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'edit', 'position', editPos, 'string', contrast_vector, 'tag', ...
    'contrast_vector', 'fontsize', UI_FS - 1);

%% add ok
okPos = [0.5 - 0.5*okWidth, promptPos(2) - 1.25*yOffset - 0.5*okHeight, okWidth, okHeight];
okPos(2) = okPos(2) - 0.5*okPos(4);
icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', okPos, 'string', 'Done', 'tag', 'done_button', 'fontsize',...
    UI_FS - 1, 'callback', {@doneCallback, InputHandle});

try
    waitfor(InputHandle);
catch
end


param_file = fullfile(sesInfo.userInput.pwd, [sesInfo.userInput.prefix, '_ica_parameter_info.mat']);
load(param_file);


if (isempty(sesInfo.userInput.dataInfo))
    error('Data is not selected');
end



function addSessions(hObject, event_data, figH)
%% Add Covariates
%

icatb_defaults;
global UI_FS;

figureTag = 'add_sess_fig';
sesFigHandle = findobj('tag', figureTag);
if (~isempty(sesFigHandle))
    delete(sesFigHandle);
end

sesInfo = get(figH, 'userdata');

listH = findobj(figH, 'tag', 'sessions');


sessName = '';
files = '';

if (listH == hObject)
    if (~strcmpi(get(figH, 'selectionType'), 'open'))
        return;
    end
    val = get(listH, 'value');
    try
        sessName = sesInfo.userInput.dataInfo(val).name;
        files = sesInfo.userInput.dataInfo(val).files;
    catch
    end
end


sesFigHandle = icatb_getGraphics('Select Files', 'normal', figureTag);
set(sesFigHandle, 'menubar', 'none');

promptWidth = 0.52; controlWidth = 0.32; promptHeight = 0.05;
xOffset = 0.02; yOffset = promptHeight; yPos = 0.9;
okWidth = 0.12; okHeight = promptHeight;

%% Covariate name and value
promptPos = [xOffset, yPos - 0.5*yOffset, promptWidth, promptHeight];
textH = icatb_uicontrol('parent', sesFigHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Enter Session Name', 'fontsize', UI_FS - 1);
icatb_wrapStaticText(textH);

promptPos = get(textH, 'position');

editPos = promptPos;
editPos(1) = editPos(1) + editPos(3) + xOffset;
editPos(3) = controlWidth;
editH = icatb_uicontrol('parent', sesFigHandle, 'units', 'normalized', 'style', 'edit', 'position', editPos, 'string', sessName, 'tag', 'session_name', 'fontsize', UI_FS - 1);

editTextPos = get(textH, 'position');

editTextHeight = 0.38;
editTextPos(1) = 0.5 - 0.5*editTextPos(3);
editTextPos(2) = editTextPos(2)  - yOffset - editTextHeight;
editTextPos(4) = editTextHeight;

% Plot edit
editTextH = icatb_uicontrol('parent', sesFigHandle, 'units', 'normalized', 'style', 'edit', ...
    'position', editTextPos, 'String', files, 'fontsize', UI_FS - 1, ...
    'horizontalalignment', 'left', 'min', 0, 'max', 2, 'tag', 'filesList');

buttonHeight = 0.052; promptHeight = 0.052;

% Browse
buttonPos = [editTextPos(1) + editTextPos(3) + xOffset, editTextPos(2) + 0.5*editTextPos(4) - 0.5*yOffset, 0.12, buttonHeight];
buttonH = icatb_uicontrol('parent', sesFigHandle, 'units', 'normalized', 'style', 'pushbutton', ...
    'position', buttonPos, 'String', 'browse', 'fontsize', UI_FS - 1, 'horizontalalignment', 'center', 'tooltipstring', ...
    'Select files or paste files in the editbox ...', 'callback', {@selectFiles, editTextH});

buttonWidth = 0.12;
editTextPos(4) = buttonHeight;
editTextPos(1) = 0.5 - 0.5*buttonWidth;
editTextPos(2) = yOffset + 0.01;
editTextPos(3) = buttonWidth;

% Plot done
icatb_uicontrol('parent', sesFigHandle, 'units', 'normalized', 'style', 'pushbutton', ...
    'position', editTextPos, 'String', 'Done', 'fontsize', UI_FS - 1, 'horizontalalignment', 'center', 'callback', ...
    {@fileDoneCallback, sesFigHandle});

try
    waitfor(sesFigHandle);
catch
    delete(sesFigHandle);
end

appName = 'fileAppData';
if (isappdata(0, appName))
    fileInfo = getappdata(0, appName);
    rmappdata(0, appName);
    files = fileInfo.filesList;
    sessName = fileInfo.sessName;
end

if (length(sesInfo.userInput.dataInfo) > 0)
    chk = strmatch(lower(sessName), lower(cellstr(char(sesInfo.userInput.dataInfo.name))), 'exact');
    if (~isempty(chk))
        ind = chk;
    end
end

if (~exist('ind', 'var'))
    ind = length(sesInfo.userInput.dataInfo) + 1;
end

if (~isempty(files))
    sesInfo.userInput.dataInfo(ind).name = sessName;
    sesInfo.userInput.dataInfo(ind).files = files;
    set(figH, 'userdata', sesInfo);
end


set(listH, 'string', char(sesInfo.userInput.dataInfo.name));
set(listH, 'value', 1);

function removeSessions(hObject, event_data, figH)
%% Remove Sessions
%

sesInfo = get(figH, 'userdata');
listH = findobj(figH, 'tag', 'sessions');
val = get(listH, 'value');

if (~isempty(val))
    check = icatb_questionDialog('title', 'Remove Session?', 'textbody', 'Do you want to remove the session from the list?');
    if (~check)
        return;
    end
end

try
    strs = cellstr(char(sesInfo.userInput.dataInfo.name));
    sesInfo.userInput.dataInfo(val) = [];
    strs(val) = [];
    set(listH, 'value', 1);
    set(listH, 'string', strs);
    set(figH, 'userdata', sesInfo);
catch
end


function selectFiles(hObject, event_data, handles)
%% Select files
%

parentH = get(hObject, 'parent');
fileInfo = get(parentH, 'userdata');

set(parentH, 'pointer', 'watch');

files = icatb_selectEntry('title', 'Select files ...', 'typeEntity', 'file', 'typeSelection', 'multiple', 'filter', '*ica*param*mat;*.mat');

if (~isempty(files))
    set(handles, 'string', files);
end

set(parentH, 'pointer', 'arrow');

function fileDoneCallback(hObject, event_data, handles)
%% Done callback
%

fileInfo = get(handles, 'userdata');

editH = findobj(handles, 'tag', 'filesList');

files = strtrim(deblank(get(editH, 'string')));

if (isempty(files))
    error('Files are not selected');
end

sessNameH = findobj(handles, 'tag', 'session_name');
sessName = deblank(get(sessNameH, 'string'));

if (isempty(sessName))
    error('Session name is not selected');
end

fileInfo.filesList = files;
fileInfo.sessName = sessName;

set(handles, 'userdata', fileInfo);
setappdata(0, 'fileAppData', fileInfo);
delete(handles);


function doneCallback(hObject, event_data, handles)
%% Done callback

sesInfo = get(handles, 'userdata');

set(handles, 'pointer', 'watch');

if (~isfield(sesInfo.userInput, 'dataInfo'))
    error('Data is not selected');
end

cVecH = findobj(handles, 'tag', 'contrast_vector');
contrast_vector = strtrim(get(cVecH, 'string'));
sesInfo.userInput.contrast_vector = contrast_vector;


fncVarH = findobj(handles, 'tag', 'fnc_variable_name');
fnc_variable_name = strtrim(get(fncVarH, 'string'));
sesInfo.userInput.fnc_variable_mat_file = fnc_variable_name;

param_file = fullfile(sesInfo.userInput.pwd, [sesInfo.userInput.prefix, '_ica_parameter_info.mat']);
icatb_save(param_file, 'sesInfo');

set(handles, 'pointer', 'arrow');

try
    set(handles, 'pointer', 'arrow');
    delete(handles);
catch
end
