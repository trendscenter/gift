function [axesHandle, colorbarHandle] = icatb_plotImage(varargin)
%% Plot image.

% Don't display colobar warnings for R2008a.
warning('off', 'MATLAB:colorbar:DeprecatedV6Argument');

icatb_defaults;
global FONT_COLOR;
global UI_FONTNAME;

titleText = '';
cDataMapping = 'scaled';
titleColor = 'c';

fontSize = 0.05;

axesHandle = [];
colorbarHandle = [];

% get input variables
for ii = 1:2:nargin
    if strcmpi(varargin{ii}, 'parent')
        axesHandle = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'data')
        data = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'cdatamapping')
        cDataMapping = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'axesclim')
        axesClim = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'colormap')
        cmap = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'titlecolor')
        titleColor = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'title')
        titleText = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'colorbarlim')
        colorbarLim = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'colorbarposition')
        colorbarPosition = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'colorbarminmaxtext')
        colorbarMinMaxText = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'text_left_right')
        text_left_right = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'xlabel')
        xlabelStr = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'ylabel')
        ylabelStr = varargin{ii + 1};
    end

end
% end for getting input variables

% Check the necessary vars
if ~exist('axesHandle', 'var')
    error('axes handle must be passed');
end


if ~exist('data', 'var')
    error('data must be present to plot the image');
end
% End for checking necessary vars

%%%%%%%%%%%%%%%%%%%%
% Plot image
%%%%%%%%%%%%%%%%%%%
ImageAxis = image(data, 'parent', axesHandle, 'CDataMapping', cDataMapping);

if ~exist('axesClim', 'var')
    minClim = min(colorbarLim);
    maxClim = 2*max(colorbarLim);
    %     minClim = min(data(:));
    %     maxClim = max(data(:));
    axesClim = [minClim, maxClim];
end

set(axesHandle, 'CLIM', axesClim);

if exist('cmap', 'var')
    handleFig = get(axesHandle, 'parent');
    set(handleFig, 'Colormap', cmap);
end

axis(axesHandle , 'off');

%axis(axesHandle, 'image');
title(titleText, 'color',  titleColor, 'parent', axesHandle);


if exist('text_left_right', 'var')
    if ~isempty(text_left_right)
        [axesPos] =   get(axesHandle, 'position');
        diff_xPos = colorbarPosition(1) + colorbarPosition(3) - axesPos(1) - axesPos(3);
        xPos = 1 + 0.1*diff_xPos;
        yPos = 0.5;
        text(xPos, yPos, text_left_right(1), 'units', 'normalized', 'parent', axesHandle, 'color', FONT_COLOR,  ...
            'fontsize', fontSize, 'HorizontalAlignment', 'left', 'verticalalignment', 'bottom', ...
            'FontName', UI_FONTNAME, 'fontunits', 'normalized',  'fontweight', 'bold');

        xPos = -0.05;
        text(xPos, yPos, text_left_right(2), 'units', 'normalized', 'parent', axesHandle, 'color', FONT_COLOR,  ...
            'fontsize', fontSize, 'HorizontalAlignment', 'left', 'verticalalignment', 'bottom', ...
            'FontName', UI_FONTNAME, 'fontunits', 'normalized', 'fontweight', 'bold');
    end

end

if exist('colorbarPosition', 'var')
    %%%%%%%%%%%%%%%%%%%%%
    % Plot colorbar
    %%%%%%%%%%%%%%%%%%%%%%
    colorbarHandle = colorbar('v6', 'peer', axesHandle);

    ChildH = get(colorbarHandle, 'Children');
    %set(ChildH, 'YData', axesClim);
    imInd = strmatch('image', lower(get(ChildH, 'Type')), 'exact');
    set(ChildH(imInd), 'YData', axesClim);

    set(colorbarHandle, 'YTickLabel', []);

    if exist('colorbarLim', 'var')
        set(colorbarHandle, 'YLim', colorbarLim);
    end

    set(colorbarHandle, 'YTick', []);

    set(colorbarHandle, 'position', colorbarPosition);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% Plot Text %%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if exist('colorbarMinMaxText', 'var')
        colorbarMin = deblank(colorbarMinMaxText(1, :));
        colorbarMax = deblank(colorbarMinMaxText(2, :));
        % Plot colorbar text

        % Maximum
        xPos = 0.01; yPos = 1.05;
        text(xPos, yPos, colorbarMax, 'units', 'normalized', 'parent', colorbarHandle, 'color', FONT_COLOR,  ...
            'fontsize', fontSize, 'HorizontalAlignment', 'left', ...
            'FontName', UI_FONTNAME, 'fontunits', 'normalized');

        % Minimum
        yPos = - 0.05;
        text(xPos, yPos, colorbarMin, 'units', 'normalized', 'parent', colorbarHandle, 'color', FONT_COLOR,  ...
            'fontsize', fontSize, 'HorizontalAlignment', 'left', ...
            'FontName', UI_FONTNAME, 'fontunits', 'normalized');
    end

end

fontUnits = get(axesHandle, 'fontUnits');
fontName = get(axesHandle, 'fontName');
fontSize = get(axesHandle, 'fontSize');

% Plot xlabel and ylabel
if exist('xlabelStr', 'var')
    text(0.42, -0.05, xlabelStr, 'units', 'normalized', 'Color', FONT_COLOR, 'parent', axesHandle, 'fontunits', fontUnits, ...
        'fontName', fontName, 'fontSize', fontSize);
    %xlabel(xlabelStr, 'parent', axesHandle);
end

if exist('ylabelStr', 'var')
    text(-0.05, 0.42, ylabelStr, 'units', 'normalized', 'rotation', 90, 'Color', FONT_COLOR, 'parent', axesHandle, 'fontunits', fontUnits, ...
        'fontName', fontName, 'fontSize', fontSize);
    %ylabel(ylabelStr, 'parent', axesHandle);
end