function Answer = icatb_inputDialog(varargin)
% input dialog box - Dialog box for getting the user input.
% The uicontrols are pop up and edit box.
%
% Input:
% 1. inputText - structure containg the type of uicontrol,
%    prompt string, optional arguments (default answers).
% 2. Title - title of the input dialog
%
% Output:
%
% Answer - cell array containing the answers selected by the user,
% returns {} when the cancel button is pressed
%
% Note:
% units for the uicontrols are normalized
% returns the values for the variables can be string or number (By default
% answer is string)

handle_visibility = 'on';
windowStyle = 'normal';

if nargin < 2
    error('Atleast two arguments are required');
end

if mod(nargin, 2) ~= 0
    error('Arguments must be in pairs');
end

for i = 1:2:nargin
    if strcmpi(varargin{i}, 'inputtext')
        inputText = varargin{i + 1};
    elseif strcmpi(varargin{i}, 'title')
        titleFig = varargin{i + 1};
    elseif strcmpi(varargin{i}, 'handle_visibility')
        handle_visibility = varargin{i + 1};
    elseif strcmpi(varargin{i}, 'windowstyle')
        windowStyle = varargin{i + 1};
    end
end

handle_visibility = lower(handle_visibility);

if ~exist('titleFig', 'var')
    titleFig = 'Input Dialog';
end

% number of UIcontrols excluding the push buttons
numUIcontrols = length(inputText);

if ~isfield(inputText, 'enable')
    % loop over input texts
    for ii = 1:length(inputText)
        inputText(ii).enable = 'on';
    end
end
% end for checking enable property

okTag = 'OK';
cancelTag = 'Cancel';

[figHandle] = icatb_plot_controls_fig(inputText, titleFig, handle_visibility, okTag, cancelTag);

set(figHandle, 'windowstyle', windowStyle);

okHandle = findobj(figHandle, 'tag', okTag);
cancelHandle = findobj(figHandle, 'tag', cancelTag);

% done callback
set(okHandle, 'callback', {@okCallback, figHandle});

% cancel callback
set(cancelHandle, 'callback', {@cancelCallback, figHandle});


% set the figure data
set(figHandle, 'UserData', inputText, 'keypressfcn', {@keypressCallback, figHandle, okHandle});

if strcmpi(handle_visibility, 'on')
    
    figure(figHandle);
    
    % wait for figure handle
    try
        waitfor(figHandle);
    catch
        delete(figHandle);
    end
    
else
    % execute the done callback
    okCallback(okHandle, [], figHandle);
end

Answer = {};
% check application data
if isappdata(0, 'inputAppData')
    % get the application data
    Answer = getappdata(0, 'inputAppData');
    rmappdata(0, 'inputAppData'); % remove the application data
end

drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Define function callbacks %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ok callback
function okCallback(handleObj, event_data, handles)

% OK push button callback for the input dialog
% get the figure data
inputParameters = get(handles, 'userdata');
answer = {};
% loop over number of parameters
for ii = 1:length(inputParameters)
    tag = [inputParameters(ii).tag];
    answerH = findobj(handles, 'tag', tag); % find tag
    if isfield(inputParameters(ii), 'dataType')
        dataType = inputParameters(ii).dataType; % answer type
    else
        dataType = 'string';
    end
    getStyle = inputParameters(ii).uiType;
    getString = get(answerH, 'string');  % get the answer string
    if iscell(getString)
        getString = str2mat(getString);
    end
    getValue = get(answerH, 'value'); % get the answer value
    % check the uicontrol type
    if strcmp(lower(getStyle), 'edit')
        answerString = deblank(getString);
    else
        answerString = deblank(getString(getValue, :));
    end
    % check the answer type
    if strcmp(lower(dataType), 'numeric')
        answer{ii} = str2num(answerString);
    else
        answer{ii} = answerString;
    end
end
% set the application data
setappdata(0, 'inputAppData', answer);
delete(handles);

% cancel callback
function cancelCallback(handleObj, evd, handles)

% close figure window
delete(handles);

function keypressCallback(hObject, event_data, handles, okHandle);

keyVal = get(handles, 'currentcharacter');

% if enter key is detected
if keyVal == 13
    okCallback(okHandle, [], handles);
end


function enableCallback(hObject, event_data, handles)
% enable on callback for writing complex images

getString = get(hObject, 'string'); % popup string
getValue = get(hObject, 'value'); % popup value
% get the selected string
if iscell(getString)
    selectedStr = getString{getValue};
else
    selectedStr = deblank(getString(getValue, :));
end
% end for selecting string

% reading complex data
complexReadHandle = findobj(handles, 'tag', ['answer', 'read_complex_images']);
% writing complex data
complexImHandle = findobj(handles, 'tag', ['answer', 'write_complex_images']);

% for writing complex images
if strcmpi(selectedStr, 'complex')
    
    % reading complex data
    if ~isempty(complexReadHandle)
        set(complexReadHandle, 'enable', 'on');
    else
        set(complexReadHandle, 'enable', 'off');
    end
    
    % writing complex data
    if ~isempty(complexImHandle)
        set(complexImHandle, 'enable', 'on');
    else
        set(complexImHandle, 'enable', 'off');
    end
    
else
    
    % reading complex data
    if ~isempty(complexReadHandle)
        set(complexReadHandle, 'enable', 'off');
    end
    
    % writing complex data
    if ~isempty(complexImHandle)
        set(complexImHandle, 'enable', 'off');
    end
    
end
% end for checking the complex images
