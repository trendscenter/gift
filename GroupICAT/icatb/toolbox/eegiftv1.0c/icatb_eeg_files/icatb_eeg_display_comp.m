function icatb_eeg_display_comp(varargin)
% Display EEG components

time_course_color = 'm';
titleColor = 'c';
number_per_figure = 4;

for ii = 1:2:nargin
    if strcmpi(varargin{ii}, 'plot_data')
        plotData = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'title_color')
        titleColor = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'color_map')
        cmap = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'topo_cmap')
        topo_cmap = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'time_course_color')
        time_course_color = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'number_per_figure')
        number_per_figure = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'eeglocs')
        EEGlocs = varargin{ii + 1};
    end
end

icatb_defaults;
global BG_COLOR;

% Number of components
numComponents = length(plotData.comp);

numCols = ceil(sqrt(number_per_figure));
numRows = round(number_per_figure / numCols);
% Number of figures possible
numFigures = ceil(numComponents / number_per_figure);

timeAxis = plotData.timeAxis;
xlabelStr = plotData.timeAxisXlabel;
ylabelStr = plotData.timeAxisYlabel;


% Offsets
yOffset = 0.01; xOffset = 0.08;

if number_per_figure == 1
    yOffset = 0.025;
    %xOffset = 0.08;
end

% Font size
fontSize = (0.05)*sqrt(numRows*numCols);

if fontSize > 1
    fontSize = 0.85;
end

% Box widths
boxWidth = ((1 - xOffset - numCols*xOffset) / numCols);

boxHeight = ((1 - 2*yOffset*numRows) / numRows);


startXOffset = 0.6*xOffset;
count = 0;
% Loop over  figures
for nFig = 1:numFigures
    
    % Get the graphics handle
    figHandle = icatb_getGraphics([plotData.figLabel, ' ', num2str(nFig)], 'graphics', ['Figure ', num2str(nFig)], 'on');
    
    % Loop over number of rows
    for getRow = 1:numRows
        
        boxXOrigin = startXOffset;
        % Loop over number of cols
        for getCol = 1:numCols
            
            count = count + 1;
            
            if count <= length(plotData.comp)
                
                boxYOrigin = (1 - getRow*boxHeight - getRow*yOffset);
                
                % Box Position
                boxPos = [boxXOrigin, boxYOrigin, boxWidth, boxHeight];
                
                % Draw Loading coefficients
                plotObjects(figHandle, boxPos, plotData.comp(count), timeAxis, EEGlocs, cmap, topo_cmap, titleColor, time_course_color, fontSize, ...
                    xlabelStr, ylabelStr);
            end
            
            boxXOrigin = boxXOrigin + boxWidth + xOffset;
            
        end
        % End loop over number of rows
    end
    % End loop over number of rows
    
    graphicsH(nFig).H = figHandle;
    graphicsH(nFig).userdata = get(figHandle, 'userdata');
    
    figureData.GraphicsHandle = graphicsH;
    set(figHandle, 'userdata', figureData);
    clear figureData;
    
    set(figHandle, 'color', BG_COLOR);
end
% End loop over figures


% Plot previous, next and exit push buttons
icatb_plotNextPreviousExitButtons(graphicsH);

for nFig = 1:numFigures
    set(graphicsH(nFig).H, 'WindowButtonDownFcn', @mouseClickCallback);
end

function plotObjects(figHandle, boxPos, dataStruct, timeAxis, EEGlocs, cmap, topo_cmap, titleColor, time_course_color, fontSize, xlabelStr, ylabelStr)
% Use boxposition and draw objects

icatb_defaults;

global AXES_COLOR;
global FONT_COLOR;
global UI_FONTNAME;

userdata = get(figHandle, 'userdata');

xOffset = 0.05;
yOffset = 0.06;
startYOffset = 0.03;

% numFeatures
numAxes = length(dataStruct.axes);

axes1Height = 0.5*(boxPos(4) - 2*yOffset - startYOffset);
axes2Height = 0.5*(boxPos(4) - 2*yOffset - startYOffset);

axesWidth = boxPos(3) - xOffset;
% Loading coefficient axes
axesOb(1).pos = [boxPos(1) + xOffset, boxPos(2) + boxPos(4) - startYOffset - axes1Height, axesWidth, axes1Height];

axesWidth = 0.45*axesOb(1).pos(3);
axesWidth = min([axesWidth, axes2Height]); axes2Height = axesWidth;
axesOb(2).pos = [axesOb(1).pos(1), boxPos(2) + yOffset, axesWidth, axes2Height];
axesOb(3).pos = [axesOb(2).pos(1) + axesOb(2).pos(3) + 1.85*xOffset, boxPos(2) + yOffset, axesWidth, axes2Height];

% Form a composite map
compositeMap = [topo_cmap; cmap];

countData = length(userdata);
for nA = 1:numAxes
    countData = countData + 1;
    axesH = axes('Parent', figHandle, 'units', 'normalized', 'position', axesOb(nA).pos, 'color', AXES_COLOR, ...
        'Fontname', UI_FONTNAME, 'fontunits', 'normalized', 'fontsize', fontSize);
    
    % Store user data
    userdata(countData).axesPosition = axesOb(nA).pos;
    userdata(countData).time_course_color = time_course_color;
    userdata(countData).titleColor = titleColor;
    userdata(countData).colorbarLIM = [];
    userdata(countData).axesCLIM = [];
    userdata(countData).cmap = cmap;
    userdata(countData).colorbarPos = [];
    userdata(countData).colorbarMinMaxText = [];
    userdata(countData).textLeftRight = [];
    userdata(countData).slice_plane = [];
    userdata(countData).offset = [xOffset, yOffset];
    userdata(countData).data = dataStruct.axes(nA).data;
    userdata(countData).EEGlocs = EEGlocs;
    userdata(countData).axesTitle = dataStruct.axes(nA).title;
    userdata(countData).topo_cmap = topo_cmap;
    userdata(countData).timeAxis = timeAxis;
    userdata(countData).sem = dataStruct.axes(nA).sem;
    userdata(countData).xlabelStr = xlabelStr;
    userdata(countData).ylabelStr = ylabelStr;
    
    if strcmpi(dataStruct.axes(nA).plotType, 'image')
        userdata(countData).plotType = 'image';
        % Plot Image
        axesPos = get(axesH, 'position');
        axesWidth = min([axesPos(3), axesPos(4)]);
        axesPos(3:4) = axesWidth;
        set(axesH, 'position', axesPos);
        
        % Colorbar position
        colorbarPos = get(axesH, 'position');
        colorbarPos(1) = colorbarPos(1) + colorbarPos(3) + 0.01*xOffset;
        colorbarPos(3) = 0.35*xOffset;
        
        colorbarLim = [dataStruct.minInterval dataStruct.maxInterval];
        
        % Colorbar min max text
        minClim = min(colorbarLim);
        maxClim = 2*max(colorbarLim);
        colorbarMinMaxText = str2mat(num2str(dataStruct.minICAIm), num2str(dataStruct.maxICAIm));
        
        axesCLIM = [minClim, maxClim];
        
        % Update user data
        userdata(countData).colorbarLIM = colorbarLim;
        userdata(countData).axesCLIM = axesCLIM;
        userdata(countData).cmap = cmap;
        userdata(countData).colorbarPos = colorbarPos;
        userdata(countData).colorbarMinMaxText = colorbarMinMaxText;
        
        ax2 = axesH;
        
        % Plot the image
        [axesH, colorbarHandle] = icatb_plotImage('parent', axesH, 'data', dataStruct.axes(nA).data, 'CDataMapping', 'scaled', 'axesCLIM', ...
            axesCLIM, 'colormap', cmap, 'title', dataStruct.axes(nA).title, 'titlecolor', titleColor, ...
            'colorbarPosition', colorbarPos, 'colorbarMinMaxText', ...
            colorbarMinMaxText); %, 'xlabel', 'Time Points', 'ylabel', 'Trials');
        axis(axesH, 'on');
        xlabel(xlabelStr, 'parent', axesH);
        ylabel('Trials', 'parent', axesH);
        set(axesH, 'YColor', FONT_COLOR, 'XColor', FONT_COLOR);
        if ~isempty(timeAxis)
            set(axesH, 'XTickLabel', num2str(timeAxis(get(axesH, 'XTick'))));
        end
        
        
    elseif strcmpi(dataStruct.axes(nA).plotType, 'timecourse')
        
        userdata(countData).plotType = 'timecourse';
        
        if ~isempty(timeAxis)
            %xlabelStr = 'Time (ms)';
            data = [timeAxis, dataStruct.axes(nA).data];
        else
            data = dataStruct.axes(nA).data;
            %xlabelStr = 'Time Points';
        end
        
        if ~isempty(dataStruct.axes(nA).sem)
            legendString = {'Mean + SEM', 'Mean - SEM', 'Mean'};
        else
            legendString{1} = 'Mean';
        end
        
        % Plot Timecourse
        icatb_plotTimecourse('parent', axesH, 'data', data, 'color', time_course_color, 'titleColor', ...
            titleColor, 'title', dataStruct.axes(nA).title, 'YAxisLocation', 'left', 'xlabelstr', xlabelStr, ...
            'ylabelstr', ylabelStr, 'sem', dataStruct.axes(nA).sem);
        
    elseif strcmpi(dataStruct.axes(nA).plotType, 'topoplot')
        % Colorbar min max text
        userdata(countData).plotType = 'topoplot';
        icatb_eeg_topoplot(dataStruct.axes(nA).data, EEGlocs, 'colormap', topo_cmap, 'electrodes', 'on');
        title(dataStruct.axes(nA).title, 'parent', axesH, 'color', titleColor);
        ax1 = axesH;
        
    end
    
end


%%%%% Set Axes CLim property %%%%%%%%%%%%%%%%

colormap(compositeMap);

CmLength   = size(compositeMap, 1);   % Colormap length
BeginSlot1 = 1;                  % Beginning slot
EndSlot1   = size(topo_cmap, 1); %length(cmap);    % Ending slot
BeginSlot2 = EndSlot1 + 1;
EndSlot2   = CmLength;
CLim1      = get(ax1, 'CLim');  % CLim values for each axis
CLim2      = get(ax2, 'CLim');

newCLIM1 = newclim(BeginSlot1, EndSlot1, CLim1(1), CLim1(2), CmLength);
newCLIM2 = newclim(BeginSlot2, EndSlot2, CLim2(1), CLim2(2), CmLength);

% Set new CLIM property
set(ax1, 'CLim', newCLIM1);
set(ax2, 'CLim', newCLIM2);

%%%%% End For Setting Axes CLim property %%%%%%%%%%%%%%%%


colorbarNewLim = get(colorbarHandle, 'YLim');
colorbarNewLim = [0.5*max(colorbarNewLim) + 1, max(colorbarNewLim)];
set(colorbarHandle, 'YLim', colorbarNewLim);

% Set userdata to figure
set(figHandle, 'userdata', userdata);


%%%%%%%%%%%%% Function Callbacks %%%%%%%%%%%%%%%

function mouseClickCallback(hObject, event_data, handles)

icatb_defaults;

% axes color
global AXES_COLOR;
global BG_COLOR;
global FONT_COLOR;

% FONT DEFAULTS
global UI_FONTNAME; % font name
global UI_FONTUNITS; % font units
global UI_FS;

% Identify selectionType
selectionType = get(hObject, 'SelectionType');

objTag = get(gco, 'Tag');

if strcmpi(objTag, 'legend')
    return;
end

% Object Type
objType = get(gco, 'Type');

% Identify double click
if (strcmpi(selectionType, 'normal') | strcmpi(selectionType, 'open')) & ...
        (strcmpi(objType, 'axes') | strcmpi(objType, 'image') | strcmpi(objType, 'line') | strcmpi(objType, 'surface') | ...
        strcmpi(objType, 'patch'))
    
    objData = get(hObject, 'userdata');
    
    figIndex = objData.index;
    
    userdata = objData.GraphicsHandle(figIndex).userdata;
    
    set(hObject, 'units', 'normalized');
    
    currentPoint = get(hObject, 'currentPoint');
    
    for nn = 1:length(userdata)
        
        if strcmpi(userdata(nn).plotType, 'image')
            colorbarPos = userdata(nn).colorbarPos;
            axesLocation = userdata(nn).axesPosition;
            % Include the axes until colorbar
            axesLocation(3) = colorbarPos(1) - axesLocation(1) + colorbarPos(3);
        else
            axesLocation = userdata(nn).axesPosition;
        end
        
        
        pointLiesInside = icatb_check_point_inside(currentPoint, axesLocation);
        
        if pointLiesInside
            
            data = userdata(nn).data; % data
            titleColor = userdata(nn).titleColor; % title color
            time_course_color = userdata(nn).time_course_color; % timecourse color
            axesTitle = userdata(nn).axesTitle; % axes title
            colorbarLim = userdata(nn).colorbarLIM;
            axesCLIM = userdata(nn).axesCLIM;
            cmap = userdata(nn).cmap;
            topo_cmap = userdata(nn).topo_cmap;
            colorbarPos = userdata(nn).colorbarPos;
            colorbarMinMaxText = userdata(nn).colorbarMinMaxText;
            offset = userdata(nn).offset;
            EEGlocs = userdata(nn).EEGlocs;
            timeAxis = userdata(nn).timeAxis;
            sem = userdata(nn).sem;
            xlabelStr = userdata(nn).xlabelStr;
            ylabelStr = userdata(nn).ylabelStr;
            
            % Open figure
            [graphicsHandle] = icatb_getGraphics(axesTitle, 'displaygui', axesTitle, 'on');
            
            axesH = axes('units', 'normalized', 'position', [0.15 0.15 0.75 0.75], 'color', AXES_COLOR, 'FontName', UI_FONTNAME, ...
                'fontunits', UI_FONTUNITS, 'fontsize', UI_FS);
            
            try
                set(hObject, 'pointer', 'watch');
                
                % Plot time course
                if strcmpi(userdata(nn).plotType, 'timecourse')
                    if ~isempty(timeAxis)
                        % xlabelStr = 'Time (ms)';
                        data = [timeAxis, data];
                    else
                        %xlabelStr = 'Time Points';
                    end
                    if ~isempty(sem)
                        legendString = {'Mean + SEM', 'Mean - SEM', 'Mean'};
                    else
                        legendString{1} = 'Mean';
                    end
                    % Plot timecourse
                    icatb_plotTimecourse('parent', axesH, 'data', data, 'color', time_course_color, 'titleColor', ...
                        titleColor, 'title', axesTitle, 'YAxisLocation', 'left', 'xlabelstr', xlabelStr, ...
                        'ylabelstr', ylabelStr, 'sem', sem, 'legendString', legendString);
                elseif strcmpi(userdata(nn).plotType, 'image')
                    % Colorbar position
                    colorbarPos = get(axesH, 'position');
                    colorbarPos(1) = colorbarPos(1) + colorbarPos(3) + 0.4*offset(1);
                    colorbarPos(3) = 0.5*offset(1);
                    % Plot the image
                    icatb_plotImage('parent', axesH, 'data', data, 'CDataMapping', 'scaled', 'axesCLIM', axesCLIM, ...
                        'colormap', cmap, 'title', axesTitle, 'titlecolor', titleColor, 'colorbarLIM', colorbarLim, ...
                        'colorbarPosition', colorbarPos, 'colorbarMinMaxText', colorbarMinMaxText); %, 'xlabel', 'Time Points', ...
                    %'ylabel', 'Trials');
                    axis(axesH, 'on');
                    xlabel(xlabelStr, 'parent', axesH);
                    ylabel('Trials', 'parent', axesH);
                    set(axesH, 'YColor', FONT_COLOR, 'XColor', FONT_COLOR);
                    if ~isempty(timeAxis)
                        set(axesH, 'XTickLabel', num2str(timeAxis(get(axesH, 'XTick'))));
                    end
                elseif strcmpi(userdata(nn).plotType, 'topoplot')
                    icatb_eeg_topoplot(data, EEGlocs, 'colormap', topo_cmap, 'electrodes', 'on');
                    title(axesTitle, 'parent', axesH, 'color', titleColor);
                    set(graphicsHandle, 'color', BG_COLOR);
                end
                set(hObject, 'pointer', 'arrow');
            catch
                set(hObject, 'pointer', 'arrow');
                icatb_displayErrorMsg;
            end
            
            return;
        end
    end
end



function CLim = newclim(BeginSlot,EndSlot,CDmin,CDmax,CmLength)
% Convert slot number and range to percent of colormap

PBeginSlot    = (BeginSlot - 1) / (CmLength - 1);
PEndSlot      = (EndSlot - 1) / (CmLength - 1);
PCmRange      = PEndSlot - PBeginSlot;

% Determine range and min and max of new CLim values
DataRange     = CDmax - CDmin;
ClimRange     = DataRange / PCmRange;
NewCmin       = CDmin - (PBeginSlot * ClimRange);
NewCmax       = CDmax + (1 - PEndSlot) * ClimRange;
CLim          = [NewCmin, NewCmax];
