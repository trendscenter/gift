function varargout = icatb_plot_matrix(data, xTickLabel, yTickLabel, varargin)
%% Plot matrix using pcolor
%
% Inputs:
% 1. data - 2D Matrix
% 2.

icatb_defaults;
global FONT_COLOR;
global UI_FS;
global UI_FONTNAME;
global BG_COLOR;


titleStr = '';
tag = '';
axes_title = '';
xAxisLabel = '';
yAxisLabel = '';
ytickoff = 0;
colorbarOn = 1;
textOn = 0;
cmap = jet(64);
colorbar_title = '';

if ((~exist('xTickLabel', 'var') && ~exist('yTickLabel', 'var')) || (isempty(xTickLabel) && isempty(yTickLabel)))
    
    xTickLabel = cellstr(num2str((1:size(data, 1))'));
    yTickLabel = cellstr(num2str((1:size(data, 2))'));
end

for i = 1:2:length(varargin)
    if (strcmpi(varargin{i}, 'fig_title'))
        titleStr = varargin{i + 1};
    elseif (strcmpi(varargin{i}, 'axesh'))
        axesH = varargin{i + 1};
    elseif (strcmpi(varargin{i}, 'tag'))
        tag = varargin{i + 1};
    elseif (strcmpi(varargin{i}, 'clim'))
        clim = varargin{i + 1};
    elseif (strcmpi(varargin{i}, 'cmap'))
        cmap = varargin{i + 1};
    elseif (strcmpi(varargin{i}, 'title'))
        axes_title = varargin{i + 1};
    elseif (strcmpi(varargin{i}, 'xlabel'))
        xAxisLabel = varargin{i + 1};
    elseif (strcmpi(varargin{i}, 'ylabel'))
        yAxisLabel = varargin{i + 1};
    elseif (strcmpi(varargin{i}, 'ytickoff'))
        ytickoff = varargin{i + 1};
    elseif (strcmpi(varargin{i}, 'colorbar'))
        colorbarOn =  varargin{i + 1};
    elseif (strcmpi(varargin{i}, 'texton'))
        textOn =  varargin{i + 1};
    elseif (strcmpi(varargin{i}, 'colorbar_title'))
        colorbar_title =  varargin{i + 1};
    end
end

if (~exist('axesH', 'var'))
    figH = icatb_getGraphics(titleStr, 'graphics', tag, 'on');
    set(figH, 'resize', 'on');
    axesH = axes('parent', figH, 'units', 'normalized', 'position', [0.1, 0.1, 0.8, 0.8]);
end
data(end+1, :) = data(end, :);
data(:, end+1) = data(:, end);
M = pcolor(axesH,  1:length(xTickLabel)+1, 1:length(yTickLabel)+1, data');
set(M, 'EdgeColor', 0.8*[1 1 1]);

if (exist('clim', 'var'))
    set(axesH, 'Clim', clim);
end
colormap(cmap);

C = [];
xlabel(xAxisLabel, 'parent', axesH, 'Interpreter', 'tex');
if (colorbarOn)
    C = colorbar;
    ylabel(C, yAxisLabel, 'parent', C, 'Interpreter', 'tex');
end
set(axesH,'XTick',1.5:length(xTickLabel)+1,'XTickLabel', xTickLabel);
if (~ytickoff)
    set(axesH,'YTick',1.5:length(yTickLabel)+1,'YTickLabel', yTickLabel);
else
    set(axesH,'YTick', [], 'YTickLabel', []);
end

title(axes_title, 'parent', axesH, 'Interpreter', 'tex');
axis(axesH, 'ij');
set(axesH, 'TickLength', [0 0]);
set(axesH, 'YColor', FONT_COLOR, 'XColor', FONT_COLOR);
set(C, 'YColor', FONT_COLOR, 'XColor', FONT_COLOR);
try
    ylabel(colorbar_title, 'parent', C, 'color', FONT_COLOR);
catch
end

if (textOn)
    x = linspace(0, 1, size(data, 1));
    y =  linspace(0, 1, size(data, 2));
    for ii = 1:size(data, 1)
        for jj = 1:size(data, 2)
            th= text(x(ii),y(jj), num2str(data(ii, jj), '%0.2f'));
        end
    end
end
%set(axesH, 'color', BG_COLOR);

varargout{1} = get(axesH, 'parent');
if (colorbarOn)
    varargout{2} = C;
end
