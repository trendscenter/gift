function figHandle = icatb_dialogBox(varargin)
% Purpose: Dialog box for handling large text in Group ICA Toolbox. Inputs must be in pairs.
%
% Input:
% 1. 'title' - title of the dialog box
% 2. 'textBody' - Body of the text
% 3. 'textType' - large for larger text and small for smaller text
% 4. 'plotbutton' - specify structure containing fields x, y, title.
% 5. 'windowStyle' - modal or normal.

% Check the number of args
if nargin < 1
    error('No arguments passed');
elseif mod(nargin, 2) ~= 0
    error('Arguments passed must be even in number');
end

textBody = '';
textType = 'small';
windowStyle = 'modal';

% check the variables
for i = 1:2:nargin
    if strcmp(lower(varargin{i}), 'title')
        titleVar = varargin{i + 1};
    elseif strcmp(lower(varargin{i}), 'textbody')
        textBody = varargin{i + 1};
    elseif strcmp(lower(varargin{i}), 'texttype')
        textType = varargin{i + 1};
    elseif strcmp(lower(varargin{i}), 'plotbutton')
        plot_figure = varargin{i + 1};
    elseif strcmp(lower(varargin{i}), 'windowstyle')
        windowStyle = varargin{i + 1};
    end
end

% Title doesn't exist
if ~exist('titleVar', 'var')
    titleVar = '';
end

if ~exist('plot_figure')
    plot_figure = [];
end

% Title doesn't exist
if ~exist('buttonType', 'var')
    buttonType = 'none';
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
global FG_COLOR;
global AXES_COLOR;
global BUTTON_FONT_COLOR;

global UI_FONTNAME;
global UI_FONTUNITS;
global FONT_COLOR;
global UI_FS;


% set up fonts
titleFont = 14;
textFont = 12;

% Position of the controls
axisPos = [0 0 1 1];
titlePos = [0.5 0.95];
textPos = [0.1 0.65];
xOffSet = 0.01; yOffSet = 0.04;

% ok position
okPos(3) = 0.12; okPos(4) = 0.08;
okPos(1) = 0.5 - 0.5*okPos(3); okPos(2) = yOffSet;

% title color
titleColor = [0 0.9 0.9];

keystr1 = sprintf(['keyVal = get(gcbf, ''currentcharacter''); \n']);

keystr2 = sprintf(['if keyVal == 13 \n delete(gcbf); clear keyVal; \n end\n']);

keypressCallback = [keystr1 keystr2];

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
    'PaperPositionMode', 'auto', ...
    'InvertHardcopy', 'off',...
    'Renderer', 'zbuffer',...
    'color', BG_COLOR, 'resize', 'off', 'name', ['About ' titleVar], ...
    'KeyPressFcn', keypressCallback, 'windowstyle', windowStyle, 'visible', 'off');

% change the size of the figure window here
pos = get(figHandle, 'position');
screenSize = get(0, 'screensize');
figurePos(1) = 0.5*(screenSize(1) + screenSize(3)) - 0.5*pos(3);
figurePos(2) = 0.5*(screenSize(2) + screenSize(4)) - 0.5*pos(4);
figurePos(3) = 0.95*pos(3);
figurePos(4) =  0.7*pos(4);

set(figHandle, 'position', figurePos);

% set axis handle to off
axisHandle = axes('Parent', figHandle, 'Position', axisPos, 'Visible', 'off');

% Name of the toolbox
text('units', 'normalized', 'string', titleVar, 'position', titlePos, 'fontsize', titleFont, 'HorizontalAlignment', 'center', ...
    'fontweight', 'bold', 'FontName', UI_FONTNAME, 'color', titleColor, 'FontAngle', 'italic');


textPos(2) = titlePos(2) - 0.3;

if strcmp(textType, 'small')
    text('units', 'normalized', 'string', textBody, 'position', textPos, 'fontsize', textFont, 'HorizontalAlignment', 'left', ...
        'fontweight', 'normal', 'FontName', UI_FONTNAME);
    
    okCallback = 'delete(gcbf)';
    icatb_uicontrol('parent', figHandle, 'units', 'normalized', 'style', 'pushbutton', 'string', 'OK', 'callback', ...
        okCallback, 'fontunits', UI_FONTUNITS, 'fontname', UI_FONTNAME, 'FontSize', UI_FS - 1, 'position', okPos, ...
        'ForegroundColor', BUTTON_FONT_COLOR, 'backgroundcolor', BUTTON_COLOR);
    
elseif strcmp(textType, 'large')
    position = get(gcbf, 'position');
    
    % define the position of the listbox
    position(1) = xOffSet; position(2) = okPos(2) + okPos(4) + yOffSet; position(3) = 1 - 2*xOffSet; position(4) = titlePos(2) - position(2) - yOffSet;
    
    
    handle_scroll = icatb_uicontrol('parent', figHandle, 'units', 'normalized', 'style','listbox', ...
        'position', position, 'string', textBody, 'foregroundcolor', FONT_COLOR, 'fontunits', UI_FONTUNITS, ...
        'fontname', UI_FONTNAME, 'horizontalalignment', 'left', 'backgroundcolor', BG2_COLOR, 'fontsize', UI_FS - 1);
    
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
    
    position = get(gcbf, 'position');
    
    cancelCallback = 'delete(gcbf)';
    
    if isempty(plot_figure)
        icatb_uicontrol('parent', figHandle, 'units', 'normalized', 'style', 'pushbutton', 'string', 'OK', ...
            'callback', cancelCallback, 'position', okPos, 'fontunits', UI_FONTUNITS, 'fontname', UI_FONTNAME, ...
            'FontSize', UI_FS, 'ForegroundColor', BUTTON_FONT_COLOR, 'backgroundcolor', BUTTON_COLOR);
    else
        % change the positions of the push buttons
        okPos(1) = 0.25 - 0.5*okPos(3);
        % Initialise plot position
        plotPos = okPos;
        plotPos(1) = 0.75 - 0.5*okPos(3);
        % change Ok button position
        icatb_uicontrol('parent', figHandle, 'units', 'normalized', 'style', 'pushbutton', 'string', 'OK', 'callback', cancelCallback, 'position', ...
            okPos, 'fontunits', UI_FONTUNITS, 'fontname', UI_FONTNAME, 'FontSize', UI_FS, 'ForegroundColor',BUTTON_FONT_COLOR, ...
            'backgroundcolor', BUTTON_COLOR);
        % plot button
        plotH = icatb_uicontrol('parent', figHandle, 'units', 'normalized', 'style', 'pushbutton', 'string', 'Plot', 'callback', {@plotCallback, figHandle} , 'position', ...
            plotPos, 'fontunits', UI_FONTUNITS, 'fontname', UI_FONTNAME, 'FontSize', UI_FS, 'ForegroundColor',BUTTON_FONT_COLOR, ...
            'backgroundcolor', BUTTON_COLOR, 'tag', 'plot', 'userdata', plot_figure);
    end
    
end


set(figHandle, 'visible', 'on');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% function callbacks %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotCallback(handleObj, evd, handles)

icatb_defaults;

global UI_FONTNAME;
global UI_FONTUNITS;
global UI_FS;
global FONT_COLOR;

% Plot callback
plotX = get(findobj(handles, 'tag', 'plot'), 'userdata');

x = plotX.x;

y = plotX.y;

if isfield(plotX, 'title')
    titleFig = plotX.title;
else
    titleFig = '';
end

%plotH = figure('windowstyle', 'modal', 'resize', 'off');
plotH = icatb_getGraphics('Parameter Options', 'other', 'figure');

plot(x, y, '-c');

%title(titleFig);
if ~iscell(titleFig)
    titleFig = {titleFig};
end

title(titleFig, 'color', FONT_COLOR, 'HorizontalAlignment', 'center', 'fontunits', UI_FONTUNITS, 'fontname', UI_FONTNAME, 'fontSize', UI_FS - 1);

set(gca, 'YColor', FONT_COLOR, 'XColor', FONT_COLOR, 'position', [0.1 0.1 0.8 0.8]);

axis('tight');

set(plotH, 'windowstyle', 'modal');