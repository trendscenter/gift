function icatb_compositeViewer(parameters, minimumParametersStructure)

% composite viewer is used to display multiple components overlaid on each
% other

% get defaults
icatb_defaults;
global BG_COLOR;
global FG_COLOR;
global AXES_COLOR;
global FONT_COLOR;
global LEGENDCOLOR;
global ASPECT_RATIO_FOR_SQUARE_FIGURE;
global THRESHOLD_VALUE;
global ANATOMICAL_PLANE;
global DETRENDNUMBER;
global UI_FONTNAME;
global TEXT_DISPLAY_SLICES_IN_MM;
global SMOOTHPARA;
global SMOOTHINGVALUE;

fontSizeText = 8;

modalityType = icatb_get_modality;
if strcmpi(modalityType, 'fmri')
    helpLabel = 'GIFT-Help';
else
    helpLabel = 'SBM-Help';
end

if ~exist('parameters', 'var') | exist( 'minimumParametersStructure', 'var' )
    
    files_in_zip = {};
    
    if (exist( 'minimumParametersStructure', 'var' ) & ~isempty(  minimumParametersStructure ))
        parameters = minimumParametersStructure;
        clear minimumParametersStructure;
    else
        % open the input parameters for the composite viewer
        parameters = icatb_inputParameters_display('composite viewer');
        if isempty(parameters)
            error('figure window was quit');
        end
    end
    
    imageValues = parameters.returnValue; % image values
    slicePlane = parameters.anatomicalplane; % slice plane
    sliceRange = parameters.slicerange; % slice range
    numComp = parameters.numComp; % Number of components
    compFiles = parameters.compFiles; % component files
    
    if isfield(parameters, 'compNumbers')
        compNumbers = parameters.compNumbers; % component Numbers
    else
        error('compNumbers field must be passed in parameters structure variable');
    end
    
    if length(compNumbers) > 5
        compNumbers = compNumbers(1:5);
    end
    
    threshValue = parameters.thresholdvalue; % threshold
    convertToZ = parameters.convertToZ; % convert to z scores
    structuralFile = parameters.structFile; % structural file
    parameters.componentnumber = compNumbers; % store the component numbers
    
    if isfield(parameters, 'files_in_zip')
        % get files in zip
        files_in_zip = parameters.files_in_zip;
    end
    
    % load the component files here
    [icasig, A, structuralImage, coords, HInfo, text_left_right] = icatb_loadICAData('structFile', ...
        structuralFile, 'compFiles', compFiles, 'slicePlane', slicePlane, 'sliceRange', sliceRange, 'comp_numbers', ...
        compNumbers, 'convertToZ', convertToZ, 'returnValue', imageValues, 'threshValue', ...
        threshValue, 'dataType', 'real', 'complexInfo', []);
    
    %%%%%%%%% Apply structural image parameters %%%%%%%%
    
    % delete the unzipped files
    if ~isempty(files_in_zip)
        for ii = 1:length(files_in_zip)
            drawnow;
            delete(files_in_zip{ii});
        end
    end
    
    parameters.structHInfo = HInfo; structHInfo = HInfo;
    structDIM = HInfo.DIM; icaDIM = HInfo.DIM;
    parameters.structuralImage = structuralImage;
    
    undetrendICA = A; % show the user detrend and undetrend
    icaTimecourse = icatb_detrend(A, 1, size(A, 1), DETRENDNUMBER);
    modelTimecourse = [];
    figLabel = '';
    clear prefix;
else
    structuralImage = parameters.structuralImage; % structural image
    parameters = rmfield(parameters, 'structuralImage');
    structHInfo = parameters.structHInfo;
    parameters = rmfield(parameters, 'structHInfo');
    structDIM = structHInfo.DIM; % structural image dimensions
    icasig = parameters.icasig; % component map
    icaDIM = parameters.icaDIM; % ica dimensions
    icaTimecourse = parameters.icaTimecourse; % ICA Time course
    imageValues = parameters.imagevalues; % image values
    sliceRange = parameters.sliceRange; % slice range
    slicePlane = parameters.slicePlane; % slice plane
    text_left_right = parameters.text_left_right; % text for plotting
    compNumbers = parameters.componentnumber;
end


% smooth ica time course
if strcmpi(SMOOTHPARA, 'yes')
    icaTimecourse = icatb_gauss_smooth1D(icaTimecourse, SMOOTHINGVALUE);
end

numOfComp = size(icasig, 1);
xdim = size(icasig, 2);
ydim = size(icasig, 3);
zdim = size(icasig, 4);
slices_in_mm = sliceRange;

%overlay components on to structural
StatusHandle = icatb_statusBar('init', 20, 'Loading Data','','');
[im, maxICAIM, minICAIM, minInterval, maxInterval] = icatb_overlayImages(icasig, structuralImage, icaDIM, structDIM, ...
    imageValues);
clear structuralImage;
clear icasig;
imFlat = im;
im = reshape(im, structDIM(1), structDIM(2), structDIM(3));
StatusHandle = icatb_statusBar('increment', 10);


%--get image in correct plane
if icatb_findstr(slicePlane, 'sagittal')
    im = permute(im, [2 3 1]);
    structHInfo.DIM = [structHInfo.DIM(2) structHInfo.DIM(3) structHInfo.DIM(1)];
    structHInfo.VOX = [structHInfo.VOX(2) structHInfo.VOX(3) structHInfo.VOX(1)];
end
if icatb_findstr(slicePlane,'coronal')
    im = permute(im, [1 3 2]);
    structHInfo.DIM = [structHInfo.DIM(1) structHInfo.DIM(3) structHInfo.DIM(2)];
    structHInfo.VOX = [structHInfo.VOX(1) structHInfo.VOX(3) structHInfo.VOX(2)];
end

%%%%%%%%%%%%%% Commented below lines to display slices in mm %%%%%%%%%%
%--get slices using slice range parameter
%temp = im(:, :, sliceRange);
%im = temp;
%clear temp;

%get dimensions, they may have changed after changing planes
DIM(1) = size(im,1);
DIM(2) = size(im,2);
DIM(3) = size(im,3);


%get colormap
cm = icatb_getColormap(numOfComp,imageValues,1, 'composite');

%get montage of image;
% a1b = zeros(DIM(2),DIM(1),1,DIM(3));
% temp = reshape(im,DIM(1),DIM(2),1,DIM(3));
% for j=1:DIM(3)
%     a1b(:,:,1,end-j+1) = (temp(:,end:-1:1,1,j)');
% end
% clear temp;
% temp=reshape(a1b(:,:,:,:),DIM(2),DIM(1),1,DIM(3));
% % flip the slices
% slices_in_mm = slices_in_mm(end:-1:1);
% [montageImage, numImagesX, numImagesY] = icatb_get_montage(temp,structHInfo.VOX);

% montage image
[montageImage, numImagesX, numImagesY, slices_in_mm] = icatb_returnMontage(im, [], DIM, structHInfo.VOX, ...
    sliceRange);
StatusHandle = icatb_statusBar('increment',10);
icatb_statusBar('finished', StatusHandle);

titleStr = ['Composite Viewer: Components: '];
for nn = 1:length(compNumbers)
    titleStr = [titleStr, ' ', num2str(compNumbers(nn))];
end
%--Setup figure for imagess
Tag = ['CompositeViewer'];
GraphicsHandle =icatb_getGraphics(titleStr,'graphics',Tag);
menu1H = uimenu('parent', GraphicsHandle, 'label', helpLabel);
menu2H = uimenu(menu1H, 'label', 'Composite Viewer', 'callback', 'icatb_openHTMLHelpFile(''icatb_composite_viewer.htm'');');
set(GraphicsHandle, 'Toolbar', 'figure', 'menubar', 'figure');
%drawnow;

%--Image
colormap(cm);
CLIM = [0 (numOfComp+1)*icatb_range([minInterval maxInterval])];
%axes('units','normalized','position',[.01 .01 .7 .9]);
componentAxis = axes('Parent', GraphicsHandle, 'units', 'normalized', 'position', [.05 .01 .7 .9]);
ImageAxis = image(montageImage, 'parent', componentAxis, 'CDataMapping', 'scaled');
set(componentAxis, 'clim', CLIM); % set the axis positions to the specified
%ImageAxis = imagesc(montageImage,CLIM);

dim(1) =size(montageImage,1);
dim(2) = size(montageImage,2);

if strcmpi(TEXT_DISPLAY_SLICES_IN_MM, 'on')
    
    % name the text and place it in the correct order (slices in mm).
    textCount = 0;
    yPos = 1 + dim(1) / numImagesY;
    for nTextRows = 1:numImagesY
        xPos = 1;
        for nTextCols = 1:numImagesX
            textCount = textCount + 1;
            if textCount <= DIM(3)
                txHandle(textCount) = text(xPos, yPos, num2str(round(slices_in_mm(textCount))), 'color', FONT_COLOR,  ...
                    'fontsize', fontSizeText, 'HorizontalAlignment', 'left', 'verticalalignment', 'bottom', ...
                    'FontName', UI_FONTNAME, 'parent', componentAxis);
            end
            xPos = xPos + (dim(2) / numImagesX);
        end
        % end for cols
        yPos = yPos + (dim(1) / numImagesY); % update the y position
    end
    % end for rows
    
end


%use dimensions of y and x axis to keep aspect ratio
% ie makes axis square within allocated space
imageAxisPos = get(get(ImageAxis,'parent'),'position');
yAxisRatio = imageAxisPos(4)/imageAxisPos(3);
xAxisRatio = imageAxisPos(3)/imageAxisPos(4);
if(yAxisRatio>1)
    yAxisRatio = 1;
else
    xAxisRatio = 1;
end
R= ASPECT_RATIO_FOR_SQUARE_FIGURE;
imageAxisPos = [imageAxisPos(1) imageAxisPos(2) imageAxisPos(3)*yAxisRatio imageAxisPos(4)*xAxisRatio];
imageAxisPos = [imageAxisPos(1) imageAxisPos(2) imageAxisPos(3)*R(1) imageAxisPos(4)*R(2)];
set(get(ImageAxis,'Parent'),'position',imageAxisPos);

%use actual dimensions of image to keep aspect ratio
imageAxisPos = get(get(ImageAxis,'parent'),'position');
imageXRatio = (dim(1)*structHInfo.VOX(2)) / (dim(2)*structHInfo.VOX(1));
imageYRatio = (dim(2)*structHInfo.VOX(1)) / (dim(1)*structHInfo.VOX(2));
if(imageXRatio>1)
    imageXRatio = 1;
else
    imageYRatio = 1;
end
imageAxisPos = [imageAxisPos(1) imageAxisPos(2) imageAxisPos(3)*imageYRatio imageAxisPos(4)*imageXRatio];
set(get(ImageAxis,'Parent'),'position',imageAxisPos);
axis off;

%--Colorbar
imagePos = get(get(ImageAxis,'parent'),'position');
pos = [imagePos(1)+imagePos(3)+.04 imagePos(2)+.1 .03 .13];
pos_colorbar = [0 0 0 0];

% Plot left right on figure
LRtextFont = 12;
xLimAxes = get(componentAxis, 'XLim');
yLimAxes = get(componentAxis, 'YLim');

if strcmpi(slicePlane, 'axial')
    text(xLimAxes(2) + 4, 0.5*yLimAxes(2), text_left_right(1), 'parent', componentAxis, ...
        'FontSize', LRtextFont, 'fontWeight', 'bold');
    text(xLimAxes(1) + 4 - LRtextFont, 0.5*yLimAxes(2), text_left_right(2), 'parent', componentAxis, ...
        'FontSize', LRtextFont, 'fontWeight', 'bold');
end
% end for plotting left right text

xOffset = 0.005;
for i=1:numOfComp
    
    if ((icatb_get_matlab_version <= 2013) || strcmpi(version('-release'), '2014a'))
        ax = axes('Parent', GraphicsHandle, 'position', [pos(1)+(i-1)*pos(3) pos(2) pos(3) pos(4)]);
        axis off;
        ColorbarHandle = colorbar('peer',ax);
    else
        ColorbarHandle = colorbar;
    end
    set(ColorbarHandle,'position',[pos(1)+(i-1)*pos(3)+(i-1)*xOffset pos(2) pos(3) pos(4)]);
    %pos_colorbar(1) = pos(3)*0.20 + pos(1) + (i-1)*pos(3);
    %pos_colorbar(2:4) = pos(2:4);
    %set(ColorbarHandle,'position', pos_colorbar);
    ChildH = get(ColorbarHandle,'Children');
    imInd = strmatch('image', lower(get(ChildH, 'Type')), 'exact');
    set(ChildH(imInd), 'YData', CLIM);
    YTick = [minInterval+(i-1)*maxInterval i*maxInterval];
    
    maxLabel = round(maxICAIM(i)*10)/10;
    minLabel = round(minICAIM(i)*10)/10;
    set(ColorbarHandle,'YLim',YTick);
    
    set(ColorbarHandle,'YTick',[]);
    
    if ((icatb_get_matlab_version <= 2013) || strcmpi(version('-release'), '2014a'))
        
        title(maxLabel, 'parent', ColorbarHandle, 'color', FONT_COLOR, 'Horizontalalignment', 'center');
        xlabel(minLabel, 'parent', ColorbarHandle, 'color', FONT_COLOR, 'Horizontalalignment', 'center');
        %text(.1,1.1,num2str(maxLabel));
        %text(0.1 + (i-1)*0.5*pos(3), 1.1, num2str(maxLabel));
        %text(.1,-.1,num2str(minLabel));
        %text(0.1 + (i-1)*0.5*pos(3), -0.1, num2str(minLabel));
        %drawnow;
        
    else
        
        labelsH = get(ColorbarHandle, 'label');
        set(ColorbarHandle, 'color', FONT_COLOR);
        set(labelsH, 'units', 'normalized');
        title(maxLabel, 'parent', ColorbarHandle, 'color', FONT_COLOR, 'Horizontalalignment', 'center');
        % set(labelsH, 'string', maxLabel, 'rotation', 0);
        % set(labelsH, 'position', [0.45, 1.2, 0]);
        set(labelsH, 'string', minLabel, 'rotation',0);
        set(labelsH, 'position', [0.45, -.12, 0]);
        
    end
    
    set(ColorbarHandle,'YTickLabel',[]);
    
end


%--Timecourses
y0 = imagePos(2)+imagePos(4)+.05;
deltaY = (.95-y0)/(numOfComp);
pos = [imagePos(1)+.05 y0 .8 deltaY];

meanPos = 0;
for i=1:numOfComp
    %calculate mean of icaTimecourse;
    meanA = zeros(1,size(icaTimecourse(:, i), 1)) + mean(icaTimecourse(:, i));
    
    ax = axes('Parent', GraphicsHandle, 'position', [pos(1) pos(2)+(i-1)*pos(4) pos(3) pos(4)]);
    
    meanPos = meanPos + pos(4);
    
    P = plot(icaTimecourse(:, i));
    
    % Name the components in the figure
    if exist('parameters', 'var')
        if isfield(parameters, 'componentnumber')
            icatb_legend(['IC ', num2str(parameters.componentnumber(i))]);
        end
    end
    interval=floor(size(cm,1)/(numOfComp+1));
    plotColor = cm((interval*(i-1))+(interval/2),:);
    set(P,'color',plotColor);
    hold on;
    P=plot(meanA,':');
    set(P,'color',[1 1 1]);
    hold off;
    set(ax,'Ycolor',FONT_COLOR);
    set(ax,'Xcolor',FONT_COLOR);
    %turn x labels off if plot is not first plot
    if(i~=1)
        set(ax,'XTickLabel',[]);
        set(ax,'XTick',[]);
    else
        if (strcmpi(modalityType, 'fmri'))
            xlabel('Scans');
            ylabel('Signal Units');
        else
            xlabel('Subjects');
            ylabel('ICA Loadings');
        end
    end
    
    %alternate sides of y axis
    if(mod(i,2) == 1)
        set(ax,'YAxisLocation','right');
    end
    axis('tight');
end