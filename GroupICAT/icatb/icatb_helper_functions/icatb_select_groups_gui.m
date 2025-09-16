function [groupName, groupVal] = icatb_select_groups_gui(subjectString, titleString, figTag, groupName, groupVal)
%% GUI for selecting data sets

if (~exist('groupName', 'var'))
    groupName = '';
end

if (~exist('groupVal', 'var'))
    groupVal = [];
end

icatb_defaults;
global UI_FS;

if ~exist(figTag, 'var')
    figTag = 'select_subjects';
end

% Delete figures that have tag select_subjects_spm_stats
if ~isempty(findobj(0, 'tag', figTag))
    delete(findobj(0, 'tag', figTag));
end

figTitle = ['Select data sets for ', titleString];

% Plot Figure
[graphicsHandle] = icatb_getGraphics(figTitle, 'normal', figTag, 'off');
set(graphicsHandle, 'menubar', 'none');

if ispc
    set(graphicsHandle, 'windowstyle', 'modal');
end

% Offsets
xOffset = 0.05; yOffset = 0.05;
editTextHeight = 0.05; editTextWidth = 0.45; yPos = 0.95;

editTextPos = [xOffset, yPos - yOffset - editTextHeight, editTextWidth, editTextHeight];

% Plot Text
editTextH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'text', ...
    'position', editTextPos, 'String', ['Enter name for ', titleString], 'fontsize', UI_FS - 1, ...
    'horizontalalignment', 'center');

[editTextH] = icatb_wrapStaticText(editTextH);

editTextPos = get(editTextH, 'position');

%%% Edit Box
editWidth = 0.4;
editPos = editTextPos;
editPos(1) = editPos(1) + editPos(3) + xOffset;
editPos(3) = editWidth;

% Plot edit control
editH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'edit', ...
    'position', editPos, 'String', groupName, 'fontsize', UI_FS - 1, 'tag', 'group_name');

% List text width and height
listTextWidth = 0.6; listTextHeight = 0.05;
listTextPos = [0.5 - 0.5*listTextWidth, editTextPos(2) - 2*yOffset - 0.5*listTextHeight, listTextWidth, ...
    listTextHeight];

% Plot listbox
listTextH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'text', ...
    'position', listTextPos, 'String', 'Select data sets', 'fontsize', UI_FS - 1);


[listTextH] = icatb_wrapStaticText(listTextH);

listTextPos = get(listTextH, 'position');

% Listbox position
listWidth = listTextWidth; listHeight = 0.4;
listPos = listTextPos;
listPos(2) = listPos(2) - yOffset - listHeight;
listPos(3) = listWidth;
listPos(4) = listHeight;

% Plot listbox
listH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'listbox', ...
    'position', listPos, 'String', subjectString, 'fontsize', UI_FS - 1, 'tag', 'subject_listbox', 'min', 0, ...
    'max', 2, 'value', groupVal);

editWidth = 0.3;
editHeight = 0.05;
editPos = [listPos(1) + 0.5*listPos(3) - 0.5*editWidth, listPos(2) - yOffset - 0.5*editHeight, ...
    editWidth, editHeight];

% Plot edit control
editH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'edit', ...
    'position', editPos, 'String', num2str(get(listH, 'value')), 'fontsize', UI_FS - 1, 'tag', ...
    'subject_editbox');


set(listH, 'callback', {@listSubjectsCallback, editH});
set(editH, 'callback', {@editSubjectsCallback, listH});

% Ok button
okWidth = 0.2; okHeight = 0.05;
okPos = [0.5 - 0.5*okWidth, yOffset + 0.5*okHeight, okWidth, okHeight];
okH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'pushbutton', ...
    'position', okPos, 'String', 'OK', 'fontsize', UI_FS - 1, 'tag', 'ok_select_subjects', 'callback', ...
    {@okSubjectsCallback, graphicsHandle,listH,  titleString});

try
    % Set graphics handle visibility on
    set(graphicsHandle, 'visible', 'on');
    waitfor(graphicsHandle);
catch
end

if isappdata(0, 'groupValData')
    groupData = getappdata(0, 'groupValData');
    groupName = groupData.groupName;
    groupVal = groupData.groupVal;
    rmappdata(0, 'groupValData');
end


function okSubjectsCallback(hObject, event_data, handles, listH, titleString)
% Ok button for selecting data sets

try
    
    groupNameH = findobj(handles, 'tag', 'group_name');
    
    groupName = get(groupNameH, 'string');
    
    groupVal = get(listH, 'value');
    
    if isempty(groupVal)
        error(['Data sets are not selected for ', titleString]);
    end
    
    groupData.groupName = groupName;
    groupData.groupVal = groupVal;
    
    setappdata(0, 'groupValData', groupData);
    
    delete(handles);
    
catch
    disp(lasterr);
    icatb_errorDialog(lasterr, 'Error Found');
end

function editSubjectsCallback(hObject, event_data, subjectListH)
% Set the listbox value

editString = deblank(get(hObject, 'string'));

editVal = str2num(editString);

numItems = size(get(subjectListH, 'string'), 1);

%%%% Do Error checking %%%%
if isempty(editVal)
    error('Error:EditBox', 'Check the edit box string (%s) \nas it must generate a valid numeric value', editString);
end

if ~isempty(find(editVal == 0))
    error('Error:EditBox', 'Edit box string (%s) cannot accept a value of 0', editString);
end

if numItems < length(editVal)
    error('Error:EditBox', 'Number of items in edit box (%d) is larger than the \nnumber of items (%d) in listbox ', ...
        length(editVal), numItems);
end

if numItems < max(editVal)
    error('Error:EditBox', 'Maximum value in editbox string (%s) is larger than the \nnumber of items (%d) in listbox ', ...
        editString, numItems);
end
%%%% End for doing error checking %%%%

checkInteger = [];
try
    checkInteger = strread(num2str(editVal), '%d');
catch
end
if isempty(checkInteger)
    error('Error:EditBox', 'Editbox string (%s) does not contain integer items', editString);
end

set(subjectListH, 'value', editVal);


function listSubjectsCallback(hObject, event_data, editH)
% Set the listbox value

editString = num2str(get(hObject, 'value'));

set(editH, 'string', editString);