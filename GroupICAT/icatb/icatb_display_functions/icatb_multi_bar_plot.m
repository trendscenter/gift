function icatb_multi_bar_plot(yAxis, varargin)
%% Show multi bar plots using a slider. Variables must be passed in pairs.
%

%% Number of parameters
numParameters = length(varargin);
if (mod(numParameters, 2) ~= 0)
    error('Input parameters after yAxis must be passed in pairs');
end

icatb_defaults;

% FONT DEFAULTS
global FONT_COLOR;
global UI_FONTNAME; % font name
global UI_FONTUNITS; % font units
global UI_FS; % font size

plotTitle = 'Multi bar plot';
xAxisLabel = 'XAxis';
yAxisLabel = 'YAxis';
axesTitle = 'Multi bar plot';


%% Loop over parameters
for i = 1:2:length(varargin)
    if strcmpi(varargin{i}, 'xticklabel')
        XTickLabel = varargin{i + 1};
    elseif strcmpi(varargin{i}, 'xlabel')
        xAxisLabel = varargin{i + 1};
    elseif strcmpi(varargin{i}, 'ylabel')
        yAxisLabel = varargin{i + 1};
    elseif strcmpi(varargin{i}, 'title')
        plotTitle = varargin{i + 1};
    elseif strcmpi(varargin{i}, 'legend')
        legendString = varargin{i + 1};
    elseif strcmpi(varargin{i}, 'axestitle')
        axesTitle = varargin{i + 1};
    elseif strcmpi(varargin{i}, 'elemperaxes')
        elemsPerAxes = varargin{i + 1};
    end
end

% Convert yAxis to a column vector 
if (size(yAxis, 1) == 1)
    yAxis = yAxis';
end

% Number of X ticks
numXTicks = size(yAxis, 1);

if (size(yAxis, 2) > 1)
    numElem = size(yAxis, 2);
else
    numElem = size(yAxis, 1);
end

figData.legendExists = 0;
if exist('legendString', 'var')
    legendString = cellstr(legendString);
    if (length(legendString) ~= numElem)
        error('Error:Legend', 'legendString length must be %d\n', numElem);
    end
    figData.legendExists = 1;
end

%% Number of plots
numPlots = ceil(numElem/elemsPerAxes);

%% Check XTickLabel
if exist('XTickLabel', 'var')
    XTickLabel = cellstr(XTickLabel);
    if (length(XTickLabel) ~= numXTicks)
        error('Error:XTick', 'XTickLabel length must be %d\n', numXTicks);
    end
else
    XTickLabel = cellstr(num2str((1:numXTicks)'));
end

%% Replicate axis titles
axesTitle = cellstr(axesTitle);
if (length(axesTitle) == 1)
    axesTitle = repmat(axesTitle, numPlots, 1);
end

%% Replicate x-axis labels
xAxisLabel = cellstr(xAxisLabel);
if (length(xAxisLabel) == 1)
    xAxisLabel = repmat(xAxisLabel, numPlots, 1);
end

%% Replicate y-axis labels
yAxisLabel = cellstr(yAxisLabel);
if (length(yAxisLabel) == 1)
    yAxisLabel = repmat(yAxisLabel, numPlots, 1);
end

% Check lengths
if (length(xAxisLabel) ~= numPlots)
    error('Error:Xlabel', 'Xlabel length must be %d\n', numPlots);
end

if (length(yAxisLabel) ~= numPlots)
    error('Error:Ylabel', 'Ylabel length must be %d\n', numPlots);
end

if (length(axesTitle) ~= numPlots)
    error('Error:AxesTitle', 'Axes title length must be %d\n', numPlots);
end


drawnow;

disp('Plotting bar graphs. Please wait ...');
fprintf('\n');
% Set the figure properties here
graphicsHandle = icatb_getGraphics(plotTitle, 'Graphics', 'multi_bar_plot', 'off');
set(graphicsHandle, 'userdata', figData);

%% Initialise vars
startT = 1;
endT = 0;
initialPosition = zeros(numPlots, 4);
xOffset = 0.075;
yOffset = 0.075;
axesWidth = (1 - 3*xOffset)/2;
axesHeight = (1 - 3*yOffset)/2;

axesPos = [xOffset, 1 - yOffset - axesHeight, axesWidth, axesHeight];
axesHandles = zeros(numPlots, 1);
initialYPos = 1 - yOffset - axesHeight;


%% Loop over plots
for nPlot = 1:numPlots
    
    endT = endT + elemsPerAxes;
    
    if (endT > numElem)
        endT = numElem;
    end
    
    if size(yAxis, 2) > 1
        currentYAxis = yAxis(:, startT:endT);
        currentXAxis = XTickLabel;
    else
        currentYAxis = yAxis(startT:endT);
        currentXAxis = XTickLabel(startT:endT);
    end
    
    
    % Row number
    row_number = ceil(nPlot/2);
    
    % Change Y axis position
    axesPos(2) = initialYPos - ((row_number - 1)*(yOffset + axesHeight));
    
    if (mod(nPlot, 2) ~= 0)
        axesPos(1) = xOffset;
    else
        axesPos(1) = 2*xOffset + axesWidth;
    end
    
    %% Bar plot
    axes_h = axes('parent', graphicsHandle, 'units', 'normalized', 'position', axesPos);
    bar_h = bar(currentYAxis);
    ylabel(yAxisLabel{nPlot}, 'parent', axes_h);
    xlabel(xAxisLabel{nPlot}, 'parent', axes_h);
    title(axesTitle{nPlot}, 'parent', axes_h);
    set(axes_h, 'YColor', FONT_COLOR, 'XColor', FONT_COLOR);
    set(axes_h, 'fontname', UI_FONTNAME, 'fontunits', UI_FONTUNITS, 'fontsize', UI_FS);
    set(axes_h, 'XTickLabel', currentXAxis);
    if exist('legendString', 'var')
        tempL = str2mat(legendString(startT:endT));
        icatb_legend(tempL);
        axisUserData.legend = tempL;
        set(axes_h, 'userdata', axisUserData);
        clear axisUserData;
    end
    
    % Store initial positions
    initialPosition(nPlot, :) = get(axes_h, 'position');
    axesHandles(nPlot) = axes_h;
    
    startT = endT + 1;
    
end
%% End loop over plots

drawnow;

%% Add Vertical Slider if necessary (axes position is less than 0.01
if (initialPosition(end, 2) < 0.01)
    sliderStep = [0.05 0.2];
    % Vertical Slider Origin & dimensions
    verSliderXOrigin = 0.97; verSliderYOrigin = 0;
    verSliderWidth = 1 - verSliderXOrigin; verSliderHeight = 1 - verSliderYOrigin;
    vertSliderPos = [verSliderXOrigin, verSliderYOrigin, verSliderWidth, verSliderHeight];
    % Store this information in userdata for both the sliders
    sliderData.initialPositions = initialPosition;
    sliderData.axesHandles = axesHandles;
    % Plot slider
    maxVal = 0;
    minVal = -abs(initialPosition(end, 2) - yOffset);
    vertSliderH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'slider', 'position', ...
        vertSliderPos, 'tag', 'verticalSlider', 'max', maxVal, 'min', minVal, 'value', maxVal, 'userdata', ...
        sliderData, 'callback', {@verSliderCallback, graphicsHandle}, 'sliderstep', sliderStep);
end

set(graphicsHandle, 'WindowButtonDownFcn', @buttonDownFcn);

set(graphicsHandle, 'visible', 'on');

fprintf('Done plotting bar graphs');
fprintf('\n');

%% Function callbacks

function verSliderCallback(handleObj, evd, handles)
%% Vertical Slider Callback
%

set(handles, 'pointer', 'watch');
figData = get(handles, 'userdata');

% Get the user data and value
sliderData = get(handleObj, 'userdata');
scrollValue = get(handleObj, 'value');

% Initial position and axes handles
initialPosition = sliderData.initialPositions;
axesHandles = sliderData.axesHandles;

% Update axes positions
figPos = get(axesHandles, 'position');
% Convert cell to MAT
figPos = cell2mat(figPos);
figPos(:, 2) = initialPosition(:, 2) - scrollValue;
% Convert figpos to cell
figPos = mat2cell(figPos, ones(size(figPos, 1), 1), 4);

% Set position to all handles
set(axesHandles, {'position'}, figPos);

%% Refresh legend handles
if (figData.legendExists)
    for nH = 1:length(axesHandles)
        icatb_legend(axesHandles(nH));
    end
end

drawnow;

set(handles, 'pointer', 'arrow');


function buttonDownFcn(hObject, event_data, handles)
%% Execute button downfcn
%

set(hObject, 'pointer', 'watch');

% Identify selectionType
selectionType = get(hObject, 'SelectionType');

cH = gco;
objType = get(gco, 'Type');

if strcmpi(selectionType, 'open')
    if ~strcmpi(objType, 'axes')
        if strcmpi(get(get(cH, 'parent'), 'Type'), 'axes')
            cH = get(cH, 'parent');
        else
            set(hObject, 'pointer', 'arrow');
            return;
        end
    end
    graphicsHandle = icatb_getGraphics('Timecourse', 'normal', 'time_course');
    axisUserData = get(cH, 'userdata');
    new_handle = copyobj(cH, graphicsHandle);
    
    set(new_handle, 'position', [0.15 0.15 0.75 0.75]);
    set(graphicsHandle, 'name', get(get(new_handle, 'Title'), 'string'));
    
    if ~isempty(axisUserData)
        try
            icatb_legend(axisUserData.legend);
        catch
        end
    end
    
end

set(hObject, 'pointer', 'arrow');