function icatb_drawMComponents(parameters)
% draw components in multiple figures
% use function callbacks to replace the strings

% load defaults
icatb_defaults;
global BG_COLOR;
global FG_COLOR;
global AXES_COLOR;
global FONT_COLOR;
global ASPECT_RATIO_FOR_SQUARE_FIGURE;
global WS;
global FONT_COLOR;
global COLORLIST;
global UI_FONTNAME;
global UI_FONTUNITS;
global TEXT_DISPLAY_SLICES_IN_MM;
global PRINTTYPE_REGRESSORS;

% get fields from parameters
minInterval = parameters.minInterval; maxInterval = parameters.maxInterval;
plotCount = parameters.imagesperfigure; images = parameters.images;
numComp = parameters.numComp; figLabel = parameters.figLabel;
structHInfo = parameters.structHInfo; cm = parameters.cm;
minICAIM = parameters.minICAIM; maxICAIM = parameters.maxICAIM;
anatomicalPlane = parameters.anatomicalplane;

printTypeRegress = lower(PRINTTYPE_REGRESSORS);
if isempty(printTypeRegress)
    printTypeRegress = 'row_wise';
end
fontSizeText = round((8 / plotCount));
% get the slice range
slices_in_mm = parameters.slicerange;

if fontSizeText < 1
    fontSizeText = 1;
end

if isfield(parameters, 'inputPrefix')
    inputPrefix = parameters.inputPrefix;
else
    inputPrefix = '';
end

% necessary flag to distinguish between the different regressors for all
% data-sets and sessions
if isfield(parameters, 'spmMatFlag')
    spmMatFlag = lower(parameters.spmMatFlag);
else
    spmMatFlag = 'not_specified';
end

% remove the corresponding fields from the structure parameters
parameters = rmfield(parameters, {'minInterval', 'maxInterval', 'imagesperfigure', 'images', 'structHInfo', ...
    'cm', 'minICAIM', 'maxICAIM', 'imagevalues', 'convertToZ', 'thresholdvalue', 'anatomicalplane', ...
    'slicerange'});

outputDir = parameters.filesOutputDir;
modalityType = icatb_get_modality;
if strcmpi(modalityType, 'fmri')
    helpLabel = 'GIFT-Help';
else
    helpLabel = 'SBM-Help';
end

% get the fields sort parameters
if isfield(parameters, 'sortParameters')
    parameters.compLabels = parameters.sortParameters.compLabels;
    parameters.icaTimecourse = parameters.sortParameters.icaTimecourse;
    parameters.modelTimecourse = parameters.sortParameters.modelTimecourse;
    parameters.undetrendICA = parameters.sortParameters.undetrendICA;
    parameters.sortParameters = rmfield(parameters.sortParameters, {'icaTimecourse', 'modelTimecourse', ...
        'undetrendICA'});
end

sortingCriteria = []; sortingType = [];
if isfield(parameters, 'sortParameters')
    sortingCriteria = lower(deblank(parameters.sortParameters.sortingCriteria));
    sortingType = lower(deblank(parameters.sortParameters.sortingType));
end

% check for different time points
if isfield(parameters, 'sortParameters')
    diffTimePoints = parameters.sortParameters.diffTimePoints;
end
% end for checking different time points

% get images fields
DIM = [size(images, 2), size(images, 3), size(images, 4)];
firstColumn = size(images, 1);
% reshape images
images = reshape(images, [firstColumn, DIM(1), DIM(2), 1, DIM(3)]);

% get the corresponding fields from the obects
if isa(images, 'complex_data')
    
    % double the number of components
    numComp = 2*numComp; plotCount = 2*plotCount;
    
    % get complexInfo
    complexInfoWrite = parameters.complexInfoWrite;
    % get the naming of the strings
    if strcmpi(complexInfoWrite.complexType, 'real&imaginary')
        firstStringPart = '(Real Part) '; secondStringPart = '(Imag Part) ';
    else
        firstStringPart = '(Mag Part) '; secondStringPart = '(Phase Part) ';
    end
    
    % original component labels
    temp_compLabels =  parameters.compLabels;
    % replicate component labels
    compLabels = repmat(struct('string', []), 1, 2*length(temp_compLabels));
    
    % time courses of class complex_data
    allTimecourses = parameters.icaTimecourse;
    allTC_first = getfield(allTimecourses, 'firstField'); allTC_second = getfield(allTimecourses, 'secondField');
    clear allTimecourses;
    
    % time courses of class complex_data
    all_undetrendICATc = parameters.undetrendICA;
    TC_first = getfield(all_undetrendICATc, 'firstField'); TC_second = getfield(all_undetrendICATc, 'secondField');
    clear all_undetrendICATc;
    
    
    % get the interval range for the magnitude or real images
    minInterval = getfield(minInterval, 'firstField');
    maxInterval = getfield(maxInterval, 'firstField');
    
    % get maxICAIM and minICAIM
    maxICAIM_first = getfield(maxICAIM, 'firstField'); minICAIM_first = getfield(minICAIM, 'firstField');
    maxICAIM_second = getfield(maxICAIM, 'secondField'); minICAIM_second = getfield(minICAIM, 'secondField');
    clear maxICAIM; clear minICAIM;
    
    % handle images
    firstSetImages = getfield(images, 'firstField'); secondSetImages = getfield(images, 'secondField');
    clear images;
    % Initialise images to 5D array
    images = zeros(2*size(firstSetImages, 1), size(firstSetImages, 2), size(firstSetImages, 3), ...
        size(firstSetImages, 4), size(firstSetImages, 5));
    
    % Initialise
    maxICAIM = zeros(1, 2*length(maxICAIM_first)); minICAIM = maxICAIM;
    
    % Initialize to zeros
    icaTimecourses = zeros(size(allTC_first, 1), 2*size(allTC_first, 2));
    undetrendICATc = icaTimecourses;
    
    % concatenate arrays
    for ii = 1:length(maxICAIM_first)
        maxICAIM(2*ii-1) = maxICAIM_first(ii); minICAIM(2*ii-1) = minICAIM_first(ii);
        maxICAIM(2*ii) = maxICAIM_second(ii); minICAIM(2*ii) = minICAIM_second(ii);
        % handle images here
        images(2*ii-1, :, :, :, :) = firstSetImages(ii, :, :, :, :);
        images(2*ii, :, :, :, :) = secondSetImages(ii, :, :, :, :);
        % handle time courses here
        icaTimecourses(:, 2*ii-1) = allTC_first(:, ii); icaTimecourses(:, 2*ii) = allTC_second(:, ii);
        undetrendICATc(:, 2*ii-1) = TC_first(:, ii); undetrendICATc(:, 2*ii) = TC_second(:, ii);
        % handle component labels
        compLabels(2*ii-1).string = [firstStringPart, temp_compLabels(ii).string];
        compLabels(2*ii).string = [secondStringPart, temp_compLabels(ii).string];
    end
    % end for concatenating arrays
    clear firstSetImages; clear secondSetImages; clear allTC_first; clear allTC_second;
    clear TC_first; clear TC_second;
    clear temp_compLabels;
    clear maxICAIM_first; clear maxICAIM_second; clear minICAIM_first; clear minICAIM_second;
    % end for loop
    
    parameters.icaTimecourse =  icaTimecourses;
    clear icaTimecourses;
    parameters.undetrendICA = undetrendICATc;
    clear undetrendICATc;
    % update parameters with the appended comp labels
    parameters.compLabels = compLabels;
    clear compLabels;
    
end
% end for checking complex_data class


% Modifying the formula to calculate the number of rows and columns
numCol = ceil(sqrt(plotCount));
%numRow = numCol;
numRow = ceil(plotCount/numCol);
numberOfFigures = ceil(numComp/plotCount);

% axis climb property (check if two objects are concatenated)
CLIM = [minInterval 2*maxInterval];
colorList = COLORLIST;
colorLSize = size(COLORLIST,2);

% description:
% Draw components in groupings of 1, 4, 9, 16, 25 in multiple figures
% plot spatial map and time course
% use commands that force spatial maps or time courses to plot in the
% specified handle

% Loop over number of figures
for nFigs = 1:numberOfFigures
    
    %--Number of the components in each figure
    if(nFigs*plotCount < numComp)
        componentIndex = (nFigs - 1)*plotCount + 1 : nFigs*plotCount;
    else
        componentIndex = (nFigs - 1)*plotCount + 1 : numComp;
    end
    
    drawnow;
    %--Setup figure for imagess
    TitleFig = [figLabel,' ',num2str(nFigs)];
    % figure tag
    Tag = ['Explorer', num2str(nFigs)];
    % store the figure
    tags{nFigs} = Tag;
    % open figure
    GraphicsHandle(nFigs).H = icatb_getGraphics(TitleFig, 'Graphics', Tag);
    set(GraphicsHandle(nFigs).H, 'Toolbar', 'figure', 'menubar', 'figure');
    menu1H = uimenu('parent', GraphicsHandle(nFigs).H, 'label', helpLabel);
    if isfield(parameters, 'htmlFile')
        htmlFile = parameters.htmlFile;
    end
    if icatb_findstr(lower(htmlFile), 'subject')
        menu2H = uimenu(menu1H, 'label', 'Subject Explorer', 'callback', ...
            'icatb_openHTMLHelpFile(''icatb_subject_explorer.htm'');');
    else
        menu2H = uimenu(menu1H, 'label', 'Component Explorer', 'callback', ...
            'icatb_openHTMLHelpFile(''icatb_component_explorer.htm'');');
    end
    %--Variables that give the position of each image within the figure
    figurePixels = get(GraphicsHandle(nFigs).H, 'position');
    pixelsXDIM = figurePixels(3);
    pixelsYDIM = figurePixels(4)-50;
    spaceAbove = (pixelsXDIM*.2)/numCol;
    spaceRight = (pixelsYDIM*.05)/numRow;
    expandRight= (pixelsXDIM*.96)/(numCol);
    expandDown = (pixelsYDIM*.96)/(numRow);
    x0=[(pixelsXDIM*.05):expandRight:(numCol+1)*expandRight];
    x1 = x0-spaceRight;
    x1=x1(2:end);
    y0=[(pixelsYDIM*.05):expandDown:(numRow+1)*expandDown];
    y0=sort(y0*-1)*-1;
    y1 = y0+spaceAbove;
    y1=y1(2:end);
    
    %--Loop through each figure, adding the correct number of images to
    %that figure
    maxXPos = 0;
    %loop over the number of components in that figure
    for compIndex = 1:size(componentIndex, 2)
        col = mod(compIndex-1, numCol)+1;
        row = ceil(compIndex / numRow);
        if row > numRow
            row = numRow;
        end
        
        if col > numCol
            col = numCol;
        end
        
        plotColor = colorList(:, mod(compIndex-1, colorLSize) + 1);
        
        % --Draw image on figure
        fontSize = round( (1/log((plotCount)+1))*10) +6;
        axisPos = round([x0(col) y0(row+1) x1(col)-x0(col) y0(row)-y1(row)]);
        % axes for plotting spatial map
        SubHandle = axes('Parent', GraphicsHandle(nFigs).H, 'units', 'pixels', 'position', ...
            axisPos, 'color', [1 1 1]);
        
        % get the montage of the image
        [im, numImagesX, numImagesY, slices_in_mm_new] = icatb_returnMontage(images, componentIndex(compIndex), ...
            DIM, structHInfo.VOX, slices_in_mm);
        dim(1) = size(im,1);
        dim(2) = size(im,2);
        imFlat = reshape(im, [dim(1)*dim(2), 1]);
        indices = find(imFlat == 0);
        imFlat(indices) = maxInterval + (icatb_range(CLIM)/size(cm,1));
        im = reshape(imFlat, [dim(1), dim(2)]);
        clear temp;
        
        
        %         use image function instead of imagesc as figure properties can be passed as
        %         input in the image function
        %         plot the image to the specified axes
        ImageAxis = image(im, 'parent', SubHandle, 'CDataMapping', 'scaled');
        set(SubHandle, 'clim', CLIM); % set the axis positions to the specified
        axis(SubHandle, 'off');
        
        
        if strcmpi(TEXT_DISPLAY_SLICES_IN_MM, 'on')
            % name the text and place it in the correct order (slices in mm).
            textCount = 0;
            yPos = 1 + dim(1) / numImagesY;
            for nTextRows = 1:numImagesY
                xPos = 1;
                for nTextCols = 1:numImagesX
                    textCount = textCount + 1;
                    if textCount <= DIM(3)
                        txHandle(textCount) = text(xPos, yPos, num2str(round(slices_in_mm_new(textCount))), 'color', FONT_COLOR,  ...
                            'fontsize', fontSizeText, 'HorizontalAlignment', 'left', 'verticalalignment', 'bottom', ...
                            'FontName', UI_FONTNAME, 'parent', SubHandle);
                    end
                    xPos = xPos + (dim(2) / numImagesX);
                end
                % end for cols
                yPos = yPos + (dim(1) / numImagesY); % update the y position
            end
            % end for rows
            
        end
        
        %         use dimensions of y and x axis to keep aspect ratio
        %         ie makes axis square within allocated space
        imageAxisPos = get(SubHandle, 'position');
        yAxisRatio = imageAxisPos(4)/imageAxisPos(3);
        xAxisRatio = imageAxisPos(3)/imageAxisPos(4);
        if(yAxisRatio>1)
            yAxisRatio = 1;
        else
            xAxisRatio = 1;
        end
        R = ASPECT_RATIO_FOR_SQUARE_FIGURE;
        imageAxisPos = [imageAxisPos(1) imageAxisPos(2) imageAxisPos(3)*yAxisRatio imageAxisPos(4)*xAxisRatio];
        imageAxisPos = [imageAxisPos(1) imageAxisPos(2) imageAxisPos(3)*R(1) imageAxisPos(4)*R(2)];
        set(SubHandle, 'position', imageAxisPos);
        
        % use actual dimensions of image to keep aspect ratio
        imageAxisPos = get(SubHandle, 'position');
        imageXRatio = (dim(1)*structHInfo.VOX(2)) / (dim(2)*structHInfo.VOX(1));
        imageYRatio = (dim(2)*structHInfo.VOX(1)) / (dim(1)*structHInfo.VOX(2));
        if(imageXRatio>1)
            imageXRatio = 1;
        else
            imageYRatio = 1;
        end
        imageAxisPos = [imageAxisPos(1) imageAxisPos(2) imageAxisPos(3)*imageYRatio imageAxisPos(4)*imageXRatio];
        set(SubHandle, 'position', imageAxisPos);
        
        % axis limits
        xLimAxes = get(SubHandle, 'XLim');
        yLimAxes = get(SubHandle, 'YLim');
        
        %%%%%%%%% display the following if it is axial view
        % Plot text here and get the text position
        textFont = 12;
        
        if strcmpi(anatomicalPlane, 'axial')
            firstLRPos = [xLimAxes(2) + 2, 0.5*yLimAxes(2)];
            secondLRPos = [xLimAxes(1) - 2 - textFont, 0.5*yLimAxes(2)];
            text(firstLRPos(1), firstLRPos(2), parameters.text_left_right(1), ...
                'FontSize', textFont, 'parent', SubHandle, 'fontWeight', 'bold');
            text(secondLRPos(1), secondLRPos(2), parameters.text_left_right(2), ...
                'FontSize', textFont, 'parent', SubHandle, 'fontWeight', 'bold');
        end
        % end for plotting text
        
        % plot
        colormap(cm);
        cAxesPos =[imageAxisPos(1)+imageAxisPos(3)+(pixelsXDIM*.005) + 2*fontSize imageAxisPos(2) ...
            min([(pixelsYDIM*.08)/plotCount (pixelsYDIM*.04)]) imageAxisPos(4)];
        
        
        % colorbar label
        maxLabel = round(maxICAIM(componentIndex(compIndex))*10)/10;
        minLabel = round(minICAIM(componentIndex(compIndex))*10)/10;
        
        
        if ((icatb_get_matlab_version <= 2013) || strcmpi(version('-release'), '2014a'))
            ax = axes('Parent', GraphicsHandle(nFigs).H, 'position', cAxesPos, 'units', 'pixels', 'color', FONT_COLOR);
            axis(ax, 'off');
            ColorbarHandle = colorbar('peer', ax);
        else
            ColorbarHandle = colorbar;
        end
        
        set(ColorbarHandle, 'units', 'pixels');
        set(ColorbarHandle, 'position', cAxesPos);
        ChildH = get(ColorbarHandle,'Children');
        %set(ChildH,'YData',CLIM);
        imInd = strmatch('image', lower(get(ChildH, 'Type')), 'exact');
        set(ChildH(imInd), 'YData', CLIM);
        set(ColorbarHandle,'YLim',[minInterval maxInterval]);
        
        if ((icatb_get_matlab_version <= 2013) || strcmpi(version('-release'), '2014a'))
            set(ColorbarHandle,'YTick',[]);
            set(ColorbarHandle,'YTickLabel',[]);
            text(4, 2, num2str(minLabel), 'units', 'data', 'FontUnits', 'points', 'FontSize', ...
                fontSize, 'FontName', UI_FONTNAME, 'parent', ColorbarHandle);
            text(4, 98, num2str(maxLabel), 'units', 'data', 'FontUnits', 'points', 'FontSize', ...
                fontSize, 'FontName', UI_FONTNAME, 'parent', ColorbarHandle);
        else
            set(ColorbarHandle, 'color', FONT_COLOR);
            labelsH = get(ColorbarHandle, 'Label');
            set (labelsH, 'units', UI_FONTUNITS, 'FontSize', fontSize, 'FontName', UI_FONTNAME);            
            YTicks = get(ColorbarHandle, 'YTick');
            set(ColorbarHandle, 'YTick', [YTicks(1), YTicks(end)]);
            set(ColorbarHandle, 'YTickLabel', char(num2str(minLabel), num2str(maxLabel)));
        end
        
        
        set(SubHandle,'XTickLabel',[]);
        set(SubHandle,'YTickLabel',[]);
        timecourseBBox = [imageAxisPos(1) imageAxisPos(2)+imageAxisPos(4) (axisPos(3)-(pixelsXDIM*.04)) ...
            y0(row)-(imageAxisPos(2)+imageAxisPos(4))];
        subHeight = (timecourseBBox(4))*.6;
        timecourseBBox = [timecourseBBox(1) timecourseBBox(2)+(.5*subHeight) timecourseBBox(3) ...
            timecourseBBox(4)-subHeight];
        timecourseAxis = axes('Parent', GraphicsHandle(nFigs).H, 'units', 'pixels', 'position', timecourseBBox);
        
        GraphicsHandle(nFigs).timeCourseBoundary(compIndex,:,:,:,:) = timecourseBBox;
        GraphicsHandle(nFigs).timecourseIndex(compIndex) = componentIndex(compIndex);
        plotColor = colorList(:,mod(compIndex-1,colorLSize)+1);
        % Plot model Time courses if any
        if ~isempty(parameters.modelTimecourse)
            % Plot all the model timecourses
            [nrows, ncols] = size(parameters.modelTimecourse);
            % normalize model time courses
            parameters.modelTimecourse = icatb_normalizeTimecourse(parameters.modelTimecourse, parameters.icaTimecourse(:, componentIndex(compIndex)), ...
                diffTimePoints, parameters.num_DataSets);
            for numModels = 1:ncols
                plot(parameters.modelTimecourse(:, numModels), 'y:', 'parent', timecourseAxis);
                hold on;
            end
            
        end
        % Line fit and regression parameters
        if strcmp(sortingCriteria, 'multiple regression') & strcmp(sortingType, 'temporal')
            % plot regression line and store regression parameter
            regressionLineFit = parameters.sortParameters.regressionLineFit{componentIndex(compIndex)};
            plot(regressionLineFit, '-.c', 'parent', timecourseAxis);
            hold on;
        end
        % plot ICA Time course to the specified time course axes
        plot(parameters.icaTimecourse(:, componentIndex(compIndex)), 'm', 'parent', timecourseAxis);
        hold off;
        set(timecourseAxis,'XTickLabel',[]);
        set(timecourseAxis,'XTick',[]);
        if (strcmpi(modalityType, 'fmri'))
            xLabelstr = 'Scans';
        else
            xLabelstr = 'Subjects';
        end
        xlabel(xLabelstr, 'FontSize',[fontSize]);
        set(timecourseAxis,'YAxisLocation','right');
        fontSize = round( (1/log((plotCount)+1))*10) +4;
        set(timecourseAxis,'FontSize',[fontSize]);
        set(timecourseAxis,'Ycolor',FONT_COLOR);
        set(timecourseAxis,'Xcolor',FONT_COLOR);
        axis(timecourseAxis, 'tight');
        title([parameters.compLabels(componentIndex(compIndex)).string], 'color', plotColor, 'FontSize', [fontSize], ...
            'parent', timecourseAxis);
    end
    
end

% set parameters as figure data
figureData = parameters;
clear parameters;

for numFig = 1:numberOfFigures
    set(GraphicsHandle(numFig).H, 'userdata', figureData);
    set(GraphicsHandle(numFig).H, 'WindowButtonDownFcn', @mouseListenerCallback);
end

icatb_plotNextPreviousExitButtons(GraphicsHandle);


%%%%%%%%%%% Function callbacks %%%%%%%%%%%%%%%%%%%%%
% mouse listener callback
function mouseListenerCallback(handleObj, event_data)


set(handleObj, 'pointer', 'watch');

% get the associated figure Data
figureData = get(handleObj, 'userdata');
compLabels = figureData.compLabels;
figNum = figureData.index;
H = figureData.GraphicsHandle(figNum).H;
numPlots = size(figureData.GraphicsHandle(figNum).timeCourseBoundary,1);
set(H, 'units','pixels');
mPos = get(H, 'currentpoint');
drawnow;
for nPlots = 1:numPlots
    tPos = figureData.GraphicsHandle(figNum).timeCourseBoundary(nPlots,:,:,:,:);
    isBetweenX = mPos(1) > tPos(1) & mPos(1) < tPos(1)+tPos(3);
    isBetweenY = mPos(2) > tPos(2) & mPos(2) < tPos(2)+tPos(4);
    if(isBetweenX & isBetweenY)
        num = figureData.GraphicsHandle(figNum).timecourseIndex(nPlots);
        % input parameters to the figure
        inputParameters = figureData;
        inputParameters = rmfield(inputParameters, 'GraphicsHandle');
        inputParameters.icaTimecourse = inputParameters.icaTimecourse(:, num);
        inputParameters.modelTimecourse = inputParameters.modelTimecourse;
        % undetrended ICA Time course
        inputParameters.undetrendICA = inputParameters.undetrendICA(:, num);
        inputParameters.titleFig = inputParameters.compLabels(num).string;
        if isfield(inputParameters, 'sortParameters')
            sortParameters = inputParameters.sortParameters;
            inputParameters = rmfield(inputParameters, 'sortParameters');
            if isfield(sortParameters, 'regressionCoeff')
                sortParameters.regressionCoeff = sortParameters.regressionCoeff{num};
            end
            if isfield(sortParameters, 'regressionLineFit')
                sortParameters.regressionLineFit = sortParameters.regressionLineFit{num};
            end
            sortParameters.sortedValues = sortParameters.sortedValues(num);
            sortParameters.sortedIndices  = sortParameters.sortedIndices(num);
            inputParameters.sortParameters = sortParameters;
        end
        % using clean function to draw time courses
        icatb_drawTimecourses(inputParameters);
    end
end

drawnow;
set(handleObj, 'pointer', 'arrow');