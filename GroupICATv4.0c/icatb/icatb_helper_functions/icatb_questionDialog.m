function [value] = icatb_questionDialog(varargin)

% Question dialog with yes and no as the buttons. Input arguments must be in pairs
% Yes means 1 and No means 0.
% Inputs:
% 1. 'title' - title of the figure
% 2. 'text body' - Body of the text

textBody = '';
okString = 'Yes';
cancelString = 'No';

% check the variables
for i = 1:2:nargin
    if strcmp(lower(varargin{i}), 'title')
        titleVar = varargin{i + 1};
    elseif strcmp(lower(varargin{i}), 'textbody')
        textBody = varargin{i + 1};
    elseif strcmp(lower(varargin{i}), 'okstring')
        okString = varargin{i + 1};
    elseif strcmp(lower(varargin{i}), 'cancelstring')
        cancelString = varargin{i + 1};
    end
end

% Title doesn't exist
if ~exist('titleVar', 'var')
    titleVar = '';
end

% text wrap detects only cell strings
if ~iscell(textBody)
    textBody = {textBody};
end


% set up the defaults

icatb_defaults;
global BG_COLOR;
global BG2_COLOR;
global BUTTON_COLOR;
global BUTTON_FONT_COLOR;
global FG_COLOR;
global AXES_COLOR;
global FONT_COLOR;
global UI_FONTNAME;
global UI_FONTUNITS;
global UI_FS;

% set up fonts
titleFont = 14;
textFont = 12;

% Position of the controls
axisPos = [0 0 1 1];
titlePos = [0.5 0.95];
textPos = [0.1 0.65];
xOffSet = 0.01; yOffSet = 0.1;

% ok position
okPos(3) = 0.12; okPos(4) = 0.08;
okPos(1) = 0.4 - 0.5*okPos(3); okPos(2) = yOffSet;

% cancel button position
cancelPos = okPos;

cancelPos(1) = 0.6 - 0.5*okPos(3);

% title color
titleColor = [0 0.9 0.9];

% set the defaults for the figure window
figHandle = figure('Resize','off', ...
    'menubar', 'none', ...
    'DefaultTextColor', FONT_COLOR,...
    'DefaultTextInterpreter', 'none',...
    'DefaultAxesColor', AXES_COLOR,...
    'DefaultAxesXColor', 'k',...
    'DefaultAxesYColor', 'k',...
    'DefaultAxesZColor', 'k',...
    'DefaultPatchFaceColor', 'k',...
    'DefaultPatchEdgeColor', 'k',...
    'DefaultSurfaceEdgeColor', 'k',...
    'DefaultLineColor', 'k',...
    'DefaultUicontrolInterruptible', 'on',...
    'PaperType', 'usletter',...
    'PaperUnits', 'normalized',...
    'PaperPositionmode', 'auto',...
    'InvertHardcopy', 'off',...
    'Renderer', 'zbuffer',...
    'color', BG_COLOR, 'resize', 'off', 'name', ['About ' titleVar], 'windowstyle', 'modal');

set(figHandle, 'keypressfcn', {@keypressCallback, figHandle});

% change the size of the figure window here
pos = get(figHandle, 'position');
screenSize = get(0, 'screensize');
figurePos(1) = 0.5*(screenSize(1) + screenSize(3)) - 0.5*pos(3);
figurePos(2) = 0.5*(screenSize(2) + screenSize(4)) - 0.5*pos(4);
figurePos(3) = 0.95*pos(3);
figurePos(4) =  0.7*pos(4);
set(figHandle, 'position', figurePos);

%set(figHandle, 'position', pos);

% set axis handle to off
axisHandle = axes('Parent', figHandle, 'Position', axisPos, 'Visible', 'off');

% Name of the toolbox
text('units', 'normalized', 'string', titleVar, 'position', titlePos, 'fontsize', titleFont, 'HorizontalAlignment', 'center', ...
    'fontweight', 'bold', 'FontName', 'times', 'color', titleColor, 'FontAngle', 'italic');


textPos(2) = titlePos(2) - 0.3;

position = get(gcf, 'position');

% define the position of the listbox
position(1) = xOffSet; position(2) = okPos(2) + okPos(4) + yOffSet; position(3) = 1 - 2*xOffSet; position(4) = titlePos(2) - position(2) - yOffSet;


handle_scroll = icatb_uicontrol('parent', figHandle, 'units', 'normalized', 'style','listbox', ...
    'position', position, 'string', textBody, 'foregroundcolor', FONT_COLOR, 'horizontalalignment', 'left', ...
    'backgroundcolor', BG2_COLOR, 'fontunits', UI_FONTUNITS, 'fontname', UI_FONTNAME, 'FontSize', UI_FS - 1);

% Apply conditions for dialog box differently for different platforms
if ispc
    set(handle_scroll, 'enable', 'inactive');
else
    set(handle_scroll, 'enable', 'on');
end


% Number of characters uicontrol can include
maxChars = 75;

% Adjust string in a listbox
[newString, newPos] = textwrap(handle_scroll, textBody, maxChars);

% Update new string
set(handle_scroll, 'String', newString);

set(handle_scroll, 'min', 0, 'max', 2);

% make no selection
set(handle_scroll, 'value', []);

okH = icatb_uicontrol('parent', figHandle, 'units', 'normalized', 'style', 'pushbutton', 'string', okString, 'callback', ...
    @okCallback , 'position', okPos, 'tag', 'yes');

cancelH = icatb_uicontrol('parent', figHandle, 'units', 'normalized', 'style', 'pushbutton', 'string', cancelString, ...
    'callback', @cancelCallback , 'position', cancelPos, 'tag', 'no');

try
    set(figHandle, 'visible','on');
    uiwait(figHandle);
catch
    if ishandle(figHandle)
        delete(figHandle);
    end
end



% set the application data here
if isappdata(0, 'ListDialogAppData')
    value = getappdata(0, 'ListDialogAppData');
    rmappdata(0, 'ListDialogAppData')
else
    value = 0;
end

drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define callbacks here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function okCallback(handleObj, evd, handles)

value = 1;
setappdata(0, 'ListDialogAppData', value)
delete(gcbf);

function cancelCallback(handleObj, evd, handles)

value = 0;
setappdata(0, 'ListDialogAppData', value)
delete(gcbf);

function keypressCallback(hObject, event_data, handles);

keyVal = get(handles, 'currentcharacter');

% if enter key is detected
if keyVal == 13
    delete(handles);
end
