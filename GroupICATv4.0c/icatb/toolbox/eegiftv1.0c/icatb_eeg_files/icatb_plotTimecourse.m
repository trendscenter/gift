function icatb_plotTimecourse(varargin)
% Plot the timecourse

icatb_defaults;
global FONT_COLOR;


titleColor = 'c';
axesTitle = '';
time_course_color = 'm';
yAxisLocation = 'left';
fontSize = 0.05;
legendString = {};

for ii = 1:2:nargin
    
    if strcmpi(varargin{ii}, 'data')
        data = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'parent')
        axesHandle = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'color')
        time_course_color = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'title')
        axesTitle = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'titlecolor')
        titleColor = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'yaxislocation')
        yAxisLocation = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'legendstring')
        legendString = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'xlabelstr')
        xlabelStr = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'ylabelstr')
        ylabelStr = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'sem')
        sem = varargin{ii + 1};
    end
    
end

if ~exist('data', 'var')
    error('Data variable must be present');
end

if size(data, 1) == 1 & size(data, 2) > 1
    data = data';
end

if size(data, 2) == 1
    xAxis = (1:length(data))';
    yAxis = data;
else
    xAxis = data(:, 1);
    yAxis = data(:, 2);
end

if exist('sem', 'var') && ~isempty(sem)
    plot(xAxis(:), (yAxis(:) + sem(:)), 'g:', 'parent', axesHandle, 'linewidth', 1.5);
    hold on;
    plot(xAxis(:), (yAxis(:) - sem(:)), 'g:', 'parent', axesHandle, 'linewidth', 1.5);
    hold on;
end


if ~exist('axesHandle', 'var')
    tag = 'Timecourse';
    [GraphicsHandle] = icatb_getGraphics(tag, 'timecourse', tag, 'on');
    plot_handle = plot(xAxis, yAxis, time_course_color, 'linewidth', 1.5);
    axesHandle = get(plot_handle, 'parent');
else
    plot(xAxis, yAxis, time_course_color, 'parent', axesHandle, 'linewidth', 1.5);
end

hold off;

% Prepare legend string
if ~isempty(legendString)
    icatb_legend(legendString{:});
end

xlabel(xlabelStr, 'parent', axesHandle);
ylabel(ylabelStr, 'parent', axesHandle);

%grid(axesHandle, 'on');
axis(axesHandle, 'tight');

title(axesTitle, 'color',  titleColor, 'parent', axesHandle);

set(axesHandle, 'YColor', FONT_COLOR, 'XColor', FONT_COLOR, 'YAxisLocation', yAxisLocation, 'box', 'on');