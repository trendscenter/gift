function refAppdata = icatb_selectRef_EventAvg(varargin)
% function used to select reference functions. Shows the selected reference
% functions at the right hand side.
% Input:
% Output:

% defaults:
typeSelection = 'single';
helpDetails.str = sprintf(['Click on the left top listbox and the corresponding reference functions will be displayed in', ...
    ' the left bottom listbox. The selected reference functions from the left bottom listbox will be displayed at the right top listbox.', ...
    ' After the selection process is done press the Ok button.']);
titleFig = 'Reference Function Window';
helpDetails.title = titleFig;

handle_visibility = 'on';

% loop over the number of args
for ii = 1:2:nargin
    if strcmp(lower(varargin{ii}), 'substring')
        sub_string = varargin{ii + 1}; % subject listbox string
    elseif strcmp(lower(varargin{ii}), 'numsubjects')
        numSubjects = varargin{ii + 1}; % number of subjects
    elseif strcmp(lower(varargin{ii}), 'numsessions')
        numSessions = varargin{ii + 1}; % number of sessions
    elseif strcmp(lower(varargin{ii}), 'help_details')
        helpDetails = varargin{ii + 1}; % help string for the help button
    elseif strcmp(lower(varargin{ii}), 'title_figure')
        titleFig = varargin{ii + 1}; % Figure tag and figure title
    elseif strcmp(lower(varargin{ii}), 'typeselection')
        typeSelection = varargin{ii + 1};
    elseif strcmp(lower(varargin{ii}), 'regressor_data')
        refData = varargin{ii + 1};
    elseif strcmp(lower(varargin{ii}), 'handle_visibility')
        handle_visibility = lower(varargin{ii + 1});
    end
end

% Check the existence of the number of subjects and sessions
if ~exist('numSubjects', 'var')
    numSubjects = 1;
end

if ~exist('numSessions', 'var')
    numSessions = 1;
end
% caculate data sets
numDataSets = numSubjects*numSessions;

% check the subject string
if ~exist('sub_string', 'var')
    count = 0;
    for ii = 1:numSubjects
        for jj = 1:numSessions
            count = count + 1;
            sub_string(count).name = ['Subject ', num2str(ii), ' Session ', num2str(jj)]; % default subject string
        end
    end
end

if strcmp(lower(typeSelection), 'single')
    maxEntries = 1;
else
    maxEntries = 2;
end

% Procedure: draw left top and bottom listboxes and right top listbox and
% beneath right top listbox draw textbox and disable the textbox but update
% the textbox with the number of entries selected.

% Load defaults
icatb_defaults;

% Screen Color Defaults
global BG_COLOR;
global BG2_COLOR;
global BUTTON_COLOR;
global BUTTON_FONT_COLOR;
global FONT_COLOR;
global AXES_COLOR;
global FONT_COLOR;
global HELP_FONT_COLOR;
% font defaults
global UI_FONTNAME;
global UI_FONTUNITS;
global UI_FS;

% create three listboxes one for subjects, one for reference functions,
% one for selected reference functions
subjectString = str2mat(sub_string.name);

clear sub_string;

% Figure
[InputHandle] = icatb_getGraphics(titleFig, 'normal', titleFig, handle_visibility);

% Set no menu bar for the figure
set(InputHandle, 'Menubar', 'None');

% tolerance
ytol = 0.03; xtol = 0.05;

%%%%%%%%%%% ok positions %%%%%%%%%%%%%%%
% define push button positions here
okWidth = 0.2; okHeight = 0.05;

okPos = [0.75 - 0.5*okWidth okHeight + 0.5*ytol  okWidth okHeight];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%% Help Position %%%%%%%%%%%%%%%%%%%
helpPos = [okPos(1) okPos(2) + okPos(4) + ytol okWidth okHeight];

%%%% selected Text Position %%%%%%%%%%%%%%%%%%%%%
selectedTextPos = helpPos;
selectedTextPos(2) = helpPos(2) + helpPos(4) + ytol;
selectedTextPos(4) = 0.05;

listHeight = 0.4;

listWidth = 0.42;

%%%%%%%%%%%%% Text height width %%%%%%%%%%%%%%%%%
textWidth = 0.3;
textHeight = 0.05;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

textWidth = listWidth;

%%%%%%%%%%%% Subject Listbox Positions %%%%%%%%%%%%%%%%%%%%%%%%%%%

% subject text position
subjectTextPos = [xtol 0.95 - 0.5*textHeight textWidth textHeight];

listboxPos = selectedTextPos(4) + selectedTextPos(2) + ytol; %okPos(4) + okPos(2) + ytol;

subjectPos = [xtol  subjectTextPos(2) - listHeight - ytol listWidth listHeight];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% Reference Listbox Positions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reference Text Position
refTextPos = [subjectPos(1) subjectPos(2) - ytol - textHeight textWidth textHeight];

referencePos = [xtol  ytol listWidth (refTextPos(2) - 2*ytol)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% Selected Reference Listbox Positions
% Selected Reference Text Position
selRefTextPos = [xtol + subjectPos(1) + subjectPos(3) subjectTextPos(2) textWidth textHeight];

% selected Text Position
selectedTextPos(1) = selRefTextPos(1);
selectedTextPos(3) = listWidth;

listHeight = (subjectTextPos(2) - listboxPos - ytol);

selReferencePos = [xtol + subjectPos(1) + subjectPos(3)  listboxPos listWidth listHeight];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%% Plot UI Controls %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot subject Text and listbox
subjectTextH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', ...
    subjectTextPos, 'horizontalalignment', 'center', 'string', 'Subjects', 'tag', 'subject_text');

subjectH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'listbox', 'position', subjectPos, ...
    'horizontalalignment', 'left', 'min', 0,  'max', 1, 'string', subjectString, 'callback', ...
    {@subjectCallback, InputHandle}, 'tag', 'subjectListbox');

% get the extent of the text under static text
getExtent = get(subjectTextH, 'extent');

if subjectTextPos(3) < getExtent(3)
    subjectTextPos(3) = getExtent(3);
end

set(subjectTextH, 'position', subjectTextPos);

% Plot Reference text and listbox
refTextH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', refTextPos, ...
    'horizontalalignment', 'center', 'string', 'Ref. Functions');
refListH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'listbox', 'position', referencePos, ...
    'horizontalalignment', 'left', 'min', 0, 'max', maxEntries, 'callback', {@refListCallback, InputHandle}, ...
    'tag', 'referenceListbox');

% get the extent of the text under static text
getExtent = get(refTextH, 'extent');

if refTextPos(3) < getExtent(3)
    refTextPos(3) = getExtent(3);
end

set(refTextH, 'position', refTextPos);

% Plot Selected text and listbox
selrefTextH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', selRefTextPos, ...
    'horizontalalignment', 'center', 'string', 'Selected Ref Function');
selrefListH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'listbox', 'position', selReferencePos, ...
    'horizontalalignment', 'left', 'tag', 'selectListbox');
% selected text box
selectedTextH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', selectedTextPos, ...
    'horizontalalignment', 'center', 'tag', 'selectedTextBox');

% get the extent of the text under static text
getExtent = get(selrefTextH, 'extent');

if selRefTextPos(3) < getExtent(3)
    selRefTextPos(3) = getExtent(3);
end

set(selrefTextH, 'position', selRefTextPos);

% Draw push button
okH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', ...
    okPos, 'string', 'OK', 'tooltipstring', 'Select Entries...', 'callback', {@okCallback, InputHandle}, 'tag', 'ok');

% Help push button
helpH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', helpPos, ...
    'string', '?', 'BackgroundColor', BUTTON_COLOR, 'ForegroundColor', HELP_FONT_COLOR, 'fontweight', 'bold', ...
    'fontunits', UI_FONTUNITS, 'fontname', UI_FONTNAME, 'FontSize', UI_FS, 'tooltipstring', 'help', 'tag', ...
    'help', 'callback', {@helpCallback, InputHandle}, 'userdata', helpDetails);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set the userdata to the figure
set(InputHandle, 'userdata', refData);

% set the application data here
for ii = 1:numSubjects*numSessions
    refAppData(ii).name = [];
end

set(okH, 'userdata', refAppData);

% % perform the callback action here
% subjectCallback(subjectH, [], InputHandle);
%
% refListCallback(refListH, [], InputHandle);

for ii = 1:numSubjects*numSessions
    set(subjectH, 'value', ii);
    % perform the callback action here
    subjectCallback(subjectH, [], InputHandle);
    refListCallback(refListH, [], InputHandle);
end

set(subjectH, 'value', 1);
% call the reference listbox again
subjectCallback(subjectH, [], InputHandle);
refListCallback(refListH, [], InputHandle);

clear refData;

clear refAppData;

%set(InputHandle, 'visible', 'on');


if strcmpi(handle_visibility, 'on')
    try
        uiwait(InputHandle);
    catch
        if ishandle(InputHandle)
            delete(InputHandle);
        end
    end
else
    % Execute the ok callback for invisible figure windows
    okCallback(okH, [], InputHandle);
end


% application data
if isappdata(0, 'refAppData')
    refAppdata = getappdata(0, 'refAppData');
    rmappdata(0, 'refAppData');
else
    error('Figure window was quit');
end

%%%%%%%%%%%%% Callbacks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define callbacks here

% Subject Listbox
function subjectCallback(handleObj, evd, handles)

% Purpose: get the string under selection and update the reference function
% listbox and the selected reference listbox

% Reference List box
refListH = findobj(handles, 'tag', 'referenceListbox');

% selection Listbox
selectH = findobj(handles, 'tag', 'selectListbox');

% Initialise reference Listbox
set(refListH, 'value', 1);
set(refListH, 'string', []);

% Initialise select String Listbox
set(selectH, 'value', 1);
set(selectH, 'string', []);

% userdata
refData = get(handles, 'userdata');

% get the values
getValue = get(handleObj, 'value');

% set the string of reference listbox to the selected subject string
set(refListH, 'string', str2mat(refData(getValue).refNames));

clear refData;

refAppData = get(findobj(handles, 'tag', 'ok'), 'userdata');
selListH = findobj(handles, 'tag', 'selectListbox');
set(selListH, 'string', refAppData(getValue).name);
set(findobj(handles, 'tag', 'selectedTextBox'), 'string', ...
    ['Selected Ref. Func: ', num2str(size(refAppData(getValue).name, 1))]);

clear refAppData;

% Reference Listbox
function refListCallback(handleObj, evd, handles)

% Purpose: get the string under selection and update the selected reference
% function

% figure data
figureData = get(handles, 'userdata');

subListbox = findobj(handles, 'tag', 'subjectListbox');

% get the string
getString = get(handleObj, 'string');

% get the values
getselValue = get(handleObj, 'value');

selectStrings = getString(getselValue, :);

% get the userdata
refAppData = get(findobj(handles, 'tag', 'ok'), 'userdata');

if ~isempty(refAppData)
    % get the subject under consideration
    getValue = get(subListbox, 'value');
    refAppData(getValue).name = selectStrings;
    refAppData(getValue).value = getselValue;

    % set the string
    set(findobj(handles, 'tag', 'selectListbox'), 'string', refAppData(getValue).name);
    % show the number of selected entries in the text box
    set(findobj(handles, 'tag', 'selectedTextBox'), 'string', ['Selected Ref. Func: ', ...
        num2str(length(getselValue))]);
    set(findobj(handles, 'tag', 'ok'), 'userdata', refAppData);
end

% Ok callback
function okCallback(handleObj, evd, handles)

% get the selected reference functions and delete the current figure

try
    % get the userdata
    refAppData = get(handleObj, 'userdata');

    subListH = findobj(handles, 'tag', 'subjectListbox'); % subject listbox

    maxEntries = get(findobj(handles, 'tag', 'referenceListbox'), 'max');

    subListStrings = get(subListH, 'string'); % get the strings of the subject listbox

    subjectTextH = findobj(handles, 'tag', 'subject_text'); % get the subject text

    % check if any of the selected model time courses are not selected
    for ii = 1:size(refAppData, 2)
        if isempty(refAppData(ii).name)
            error(['Reference Functions are not selected for data set ', subListStrings(ii, :)]);
        end
        % check if it is multiple regression criteria
        if ii > 1 & maxEntries > 1
            previousLength = size(refAppData(ii - 1).name, 1);
            currentLength = size(refAppData(ii).name, 1);
            if previousLength ~= currentLength
                error(['Select same number of regressors for data set ', subListStrings(ii, :)]);
            end
        end
    end

    % set the application data
    setappdata(0, 'refAppData', refAppData);
    set(handleObj, 'userdata', []);
    delete(handles);
catch
    icatb_errorDialog(lasterr, 'Reference Function Err');
end

function helpCallback(hObject, evd, handles)

% help button callback
helpDetails = get(hObject, 'userdata');

str = helpDetails.str;
titleString = helpDetails.title;

figHandle = icatb_dialogBox('title', titleString, 'textBody', str, 'textType', 'large');

waitfor(figHandle);

function loadCall(hObject, evd, handles)
% load callback:
% loads the reference functions from a file:
% if inexact macth is en
