function icatb_orthoViewer(parameters, minimumParametersStructure)

icatb_defaults;
global DETRENDNUMBER;
global COLORMAP_FILE;
global THRESHOLD_VALUE;
global ANATOMICAL_PLANE;
global SMOOTHPARA;
global SMOOTHINGVALUE;

complexInfoRead = [];
complexInfoWrite = [];
dataType = 'real';

zipContents.zipFiles = {};
zipContents.files_in_zip(1).name = {};
outputDir = pwd;

files_in_zip = {};

modalityType = icatb_get_modality;
if strcmpi(modalityType, 'fmri')
    helpLabel = 'GIFT-Help';
else
    helpLabel = 'SBM-Help';
end

% parameters structure is used
if ~exist('parameters', 'var') | exist( 'minimumParametersStructure', 'var' )
    
    files_in_zip = {};
    
    if (exist( 'minimumParametersStructure', 'var' ) & ~isempty(  minimumParametersStructure ))
        parameters = minimumParametersStructure;
        clear minimumParametersStructure;
    else
        % open the input parameters for the Ortho Viewer
        parameters = icatb_inputParameters_display('orthogonal viewer');
        if isempty(parameters)
            error('figure window was quit');
        end
    end
    
    imageValues = parameters.returnValue; % image values
    threshValue = parameters.thresholdvalue; % threshold
    compFiles = parameters.compFiles; % component file
    compNumber = parameters.compNumbers; % component number
    compNumber = compNumber(1);
    convertToZ = parameters.convertToZ; % convert to z scores
    
    if isfield(parameters, 'fmriFiles')
        fmriFiles(1).ses(1).name = parameters.fmriFiles; % fmri data
    else
        error('fmriFiles field is not passed in parameters structure variable');
    end
    
    diffTimePoints = icatb_get_countTimePoints(fmriFiles(1).ses(1).name);
    structuralFile = parameters.structFile;
    [path1, comp_name] = fileparts(deblank(compFiles(1, :)));
    underScorePos = icatb_findstr(comp_name, '_');
    parameters.title = [comp_name(1:underScorePos(end)), num2str(compNumber)];
    figLabel = parameters.title;
    if isfield(parameters, 'files_in_zip')
        % get files in zip
        files_in_zip = parameters.files_in_zip;
    end
    % output directory
    if isfield(parameters, 'filesOutputDir')
        outputDir = parameters.filesOutputDir;
    end
    clear parameters;
    
    % load the component files here
    [icasig, A, structuralImage, real_world_coords, HInfo, text_left_right] = icatb_loadICAData('structFile', ...
        structuralFile, 'compFiles', compFiles, 'comp_numbers', compNumber, 'convertToZ', convertToZ, ...
        'returnValue', imageValues, 'threshValue', threshValue, 'dataType', dataType, 'complexInfo', ...
        complexInfoWrite);
    
    icaHInfo = HInfo; structHInfo = HInfo;
    icasig = reshape(icasig, HInfo.DIM(1), HInfo.DIM(2), HInfo.DIM(3));
    icaHInfo = HInfo;
    icaFiles = compFiles;
    clear HInfo;
    fmriHInfo = icaHInfo;
    undetrendICA = A; % show the user detrend and undetrend
    A = icatb_detrend(A, 1, size(A, 1), DETRENDNUMBER); %detrend(A); % detrended time course
    icaTimecourse = A;
    modelTimecourse = [];
else
    % get the corresponding parameters
    structuralImage = parameters.structuralImage; % anatomical image
    parameters = rmfield(parameters, 'structuralImage'); % remove the structural image
    structHInfo = parameters.structHInfo; % header information of the anatomical image
    icasig = parameters.icasig; % component image
    parameters = rmfield(parameters, 'icasig'); % remove the icasig
    icaHInfo = parameters.icaHInfo; % ICA image header information
    A = parameters.icaTimecourse; % ICA timecourse
    imageValues = parameters.imagevalues; % image values
    icaFiles = parameters.icaFiles; % ICA files
    fmriFiles = parameters.fmriFiles; % functional files
    fmriHInfo = parameters.fmriHInfo;
    compNumber = parameters.compNum;
    figLabel = parameters.title;
    % get structural file
    structuralFile = parameters.structFile;
    % defaults
    convertToZ = parameters.convertToZ;
    threshValue = parameters.thresholdvalue;
    complexInfoRead = parameters.complexInfoRead;
    complexInfoWrite = parameters.complexInfoWrite;
    dataType = parameters.dataType;
    if isfield(parameters, 'zipContents')
        zipContents = parameters.zipContents;
    end
    if isfield(parameters, 'outputDir')
        outputDir = parameters.outputDir;
    end
    diffTimePoints = parameters.diffTimePoints;
    real_world_coords = parameters.real_world_coords;
    % get files in zip
    %files_in_zip = parameters.files_in_zip;
    clear parameters;
end

checkTp = find(diffTimePoints ~= diffTimePoints(1));

if ~isempty(checkTp)
    flagTimePoints = 'different_time_points';
else
    flagTimePoints = 'same_time_points';
end

% display maximum voxel
[maxVoxel, maxVoxPos] = max(icasig(:));

[maxVoxX, maxVoxY, maxVoxZ] = ind2sub(size(icasig), maxVoxPos);

[maximumVoxelRealPos, maximumVoxelPos] = getCoord(fmriHInfo.V(1), real_world_coords, [maxVoxX, maxVoxY, maxVoxZ]);


% maximumVoxelRealPos = squeeze(real_world_coords(maxVoxX, maxVoxY, maxVoxZ, :));
% maximumVoxelRealPos = maximumVoxelRealPos(:)';
%
% maximumVoxelPos = icatb_real_to_voxel(fmriHInfo.V(1), maximumVoxelRealPos);

%maximumVoxelPos = [maxVoxX, maxVoxY, maxVoxZ];

disp(['Maximum voxel is at ', num2str(maximumVoxelPos), ' and its value is ', num2str(maxVoxel)]);
disp(['Real world coordinates are ', num2str(maximumVoxelRealPos)]);
fprintf('\n');

% make A a column vector
if size(A, 1) == 1
    A = A';
end

% smooth ica time course
if strcmpi(SMOOTHPARA, 'yes')
    A = icatb_gauss_smooth1D(A, SMOOTHINGVALUE);
end

%% Minimum voxel position
[minVoxel, minVoxPos] = min(icasig(:));
maxAbsValue = max([abs(minVoxel), maxVoxel]);
minICATC = (minVoxel*A) / maxAbsValue;


% show only the first field
if isa(icasig, 'complex_data')
    icasig = getfield(icasig, 'firstField');
end

if isa(A, 'complex_data')
    A = getfield(A, 'firstField');
end

if isa(icasig, 'complex_data')
    icasig = getfield(icasig, 'firstField');
end

load(COLORMAP_FILE);

%  save pre interpolation ica dim
origICADIM(1) = icaHInfo.V(1).dim(1);
origICADIM(2) = icaHInfo.V(1).dim(2);
origICADIM(3) = icaHInfo.V(1).dim(3);

% overlay functional image on top of structural
StatusHandle = icatb_statusBar('init', 100, 'Loading Data For Ortho Viewer', '', '');
icaDIM = icaHInfo.DIM;
structDIM = structHInfo.DIM;
tempIC = reshape(icasig, 1, icaDIM(1), icaDIM(2), icaDIM(3));
[im, maxICAIM, minICAIM, minInterval, maxInterval] = icatb_overlayImages(tempIC, structuralImage, icaDIM, ...
    structDIM, imageValues);
clear tempIC;
StatusHandle = icatb_statusBar('increment', 100/3);

% get unit color
unitColor = (icatb_range(im))/256;

% reshape imDIM
imDIM = structDIM;
im = squeeze(reshape(im, imDIM(1), imDIM(2), imDIM(3)));

% get max value for scaling
maxICA = maxICAIM;
minICA = minICAIM;
if minICA >= 0
    icaCLIM = [minICA maxICA];
    % Special case when the maxICA is equal to zero
elseif maxICA == 10^-8 & minICA < 0
    icaCLIM = [minICA 0];
else
    tempMax = max([abs(maxICA) abs(minICA)]);
    icaCLIM = [-tempMax tempMax];
end

% reshape image so that it doesn't look "squashed"
im = reshape(im, structDIM(1), structDIM(2), structDIM(3));
imOrigin = [0 0 0];
StatusHandle = icatb_statusBar('increment', 100/3);

% get starting voxel location
% icaPos = [round(icaHInfo.DIM(1)/2) round(icaHInfo.DIM(2)/2) round(icaHInfo.DIM(3)/2)];
% realPos(1) = [round(icaPos(1)*(fmriHInfo.DIM(1)/icaHInfo.DIM(1)))];
% realPos(2) = [round(icaPos(2)*(fmriHInfo.DIM(2)/icaHInfo.DIM(2)))];
% realPos(3) = [round(icaPos(3)*(fmriHInfo.DIM(3)/icaHInfo.DIM(3)))];


% Calculate voxel origin
calcOrigin = (fmriHInfo.V(1).mat\[0 0 0 1]')';
voxelOrigin = calcOrigin(1:3);

icaPos = round(voxelOrigin);
realPos = round(voxelOrigin);

[realWorldPos, voxelCoord] = getCoord(fmriHInfo.V(1), real_world_coords, realPos);

% get starting timecourse
P = fmriFiles(1).ses(1).name;
[timecourse] = icatb_getVoxelValues(realPos, P, structuralFile, dataType, complexInfoRead, 'read');
StatusHandle = icatb_statusBar('increment', 100/3);
BOLDTimecourse = timecourse;
if (~strcmpi(icatb_get_modality, 'smri'))
    BOLDTimecourseLabel = ['fMRI Data for sub ', num2str(1), ' sess ', num2str(1), ' at Voxel ', num2str(voxelCoord)];
else
    BOLDTimecourseLabel = ['sMRI Data at Voxel ', num2str(voxelCoord)];
end
clear temp;

icatb_statusBar('finished', StatusHandle);

% figure handle
GraphicsHandle = icatb_getGraphics(figLabel, 'graphics', 'orthoviewer');
set(GraphicsHandle, 'units', 'normalized', 'Toolbar', 'figure', 'menubar', 'figure');

optionsH = uimenu('parent', GraphicsHandle, 'label', 'Options');
optionsSubH = uimenu(optionsH, 'label', 'Set Voxel Position', 'callback', {@setVoxelPosCallback, GraphicsHandle});


menu1H = uimenu('parent', GraphicsHandle, 'label', helpLabel);
menu2H = uimenu(menu1H, 'label', 'Orthogonal Viewer', 'callback', 'icatb_openHTMLHelpFile(''icatb_ortho_viewer.htm'');');


cm = icatb_getColormap(1,imageValues,1);

% save data into figure
data = struct('icaFiles', icaFiles, 'BOLDTimecourse', BOLDTimecourse, 'BOLDTimecourseLabel', BOLDTimecourseLabel, ...
    'im', im, 'imOrigin', imOrigin, 'imDIM', imDIM, 'realPos', realPos, 'icaPos', icaPos, 'voxelOrigin', voxelOrigin, ...
    'A', A, 'icasig', icasig, 'icaHInfo', icaHInfo, 'structHInfo', structHInfo, 'fmriFiles', fmriFiles, 'fmriHInfo', ...
    fmriHInfo, 'icaCLIM', icaCLIM, 'minInterval', minInterval, 'maxInterval', maxInterval, 'unitColor', unitColor, ...
    'cm', cm, 'imagevalues', imageValues, 'convertToZ', convertToZ, 'structFile', structuralFile, 'threshValue', ...
    threshValue, 'dataType', dataType, 'complexInfoRead', complexInfoRead, 'complexInfoWrite', complexInfoWrite, ...
    'maximumVoxelPos', maximumVoxelPos, 'zipContents', zipContents, 'outputDir', outputDir, 'minICATC', minICATC, ...
    'diffTimePoints', diffTimePoints, 'flagTimePoints', flagTimePoints, 'real_world_coords', real_world_coords);

data.files_to_be_deleted = files_in_zip ;
clear cm;

%save data in figure
set(GraphicsHandle, 'userdata', data);

set(GraphicsHandle, 'CloseRequestFcn', @figCloseCallback);

clear data;

redrawCallback(GraphicsHandle);

function redrawCallback(handles)
% redraw callback

%-------------------------------------------------------------------
%CODE TO Draw figure
%-------------------------------------------------------------------

%get global defaults
icatb_defaults;
global COLORLIST;
% Screen Color Defaults
global BG_COLOR;
global BG2_COLOR;
global BUTTON_COLOR;
global BUTTON_FONT_COLOR;
global FONT_COLOR;
global AXES_COLOR;
global FONT_COLOR;

% Fonts
global UI_FONTNAME;
global UI_FONTUNITS;
global UI_FS;
global LEGENDCOLOR;
global ASPECT_RATIO_FOR_SQUARE_FIGURE;
global DETRENDNUMBER;

%--load data and grab variables
data = get(handles, 'userdata');
minInterval = data.minInterval;
maxInterval = data.maxInterval;
real_world_coords = data.real_world_coords;
cm = data.cm;
A=data.A;
minICATC = data.minICATC;
%GraphicsHandle = data.GraphicsHandle;
realPos = data.realPos;
icaPos = data.icaPos;
% voxel origin
voxelOrigin = data.voxelOrigin;
unitColor = data.unitColor;
im = data.im;
fmriHInfo = data.fmriHInfo;
fmriFiles = data.fmriFiles;


numOfSub = size(fmriFiles,2);
numOfSess = size(fmriFiles(1).ses,2);


% real world position
%realWorldPos = (realPos - voxelOrigin).*(fmriHInfo.VOX);
%real_world_coords
%realWorldPos = voxel_to_real(fmriHInfo.V(1), realPos);

[realWorldPos, voxelCoord] = getCoord(fmriHInfo.V(1), real_world_coords, realPos);

imDIM = data.imDIM;
icasig = squeeze(data.icasig);
icaCLIM = data.icaCLIM;
BOLDTimecourse = data.BOLDTimecourse;
BOLDTimecourseLabel = data.BOLDTimecourseLabel;
set(handles, 'Units', 'normalized');

voxelIntensity = icasig(realPos(1), realPos(2), realPos(3));


%display current position
disp(['Current Voxel Position ', num2str(voxelCoord), ' [x y z]. Value is ', num2str(voxelIntensity)]);

% Print real world coordinates
disp(['Real world coordinate is ', num2str(realWorldPos(1)), ' ', num2str(realWorldPos(2)), ' ', num2str(realWorldPos(3))]);


%--save image
xdim = size(im,1);
ydim = size(im,2);
zdim = size(im,3);


%--get max and min value from ica image for scaling;
maxICA = icaCLIM(2);
minICA = icaCLIM(1);

%--get max and min of image
maxIM = max(max(max(im)));
minIM = min(min(min(im)));

signal = im(icaPos(1),icaPos(2),icaPos(3));
diffInterval = icatb_range([maxInterval abs(minInterval)]);
diffInterval = diffInterval / 2;
if(diffInterval - signal <=0)
    voxelColor=minInterval+unitColor;
else
    voxelColor=maxInterval-unitColor;
end

im(icaPos(1), icaPos(2),icaPos(3))=voxelColor;

sliceXY=reshape(im(:,:,icaPos(3)),size(im,1),size(im,2));
sliceXZ=reshape(im(:,icaPos(2),:),size(im,1),size(im,3));
sliceYZ=reshape(im(icaPos(1),:,:),size(im,2),size(im,3));


%scaledA from icasig;
maxAbsValue = max([max(max(max(icasig))) abs(min(min(min(icasig))))]);
A=A;
% Scaling with positive number
%scaledVal = A * (abs(icasig(icaPos(1),icaPos(2),icaPos(3)))/maxAbsValue);
scaledVal = A * (icasig(icaPos(1),icaPos(2),icaPos(3))/maxAbsValue);


%figure out max of y-axis and min y-axis for plotting ica timecourse
maxA=max(A);
minA=min(A);


newColor = cm;

%----------------------------------------------------
%Draw and scale image in Quad 1

%use dimensions of y and x axis to keep aspect ratio
% ie makes axis square within allocated space
VOX = data.structHInfo.VOX;
axQuad1Pos = [.05 .70 .45 .25];
dim = size(sliceXZ);
ImageAxis = axes('Parent', handles, 'units', 'normalized', 'position', axQuad1Pos, 'color', [1 1 1]);
axis(ImageAxis, 'off');
% using image function as it is general function than imagesc
origData.axisHandle1 = image(rot90(sliceXZ), 'parent', ImageAxis, 'CDataMapping', 'scaled');
set(ImageAxis, 'clim', [minInterval 2*maxInterval]);
%origData.axisHandle1 = imagesc(rot90(sliceXZ), [minInterval 2*maxInterval]); % replace with image function
imageAxisPos = get(ImageAxis, 'position');
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
set(ImageAxis, 'position', imageAxisPos);


%use actual dimensions of image to keep aspect ratio
imageAxisPos = get(ImageAxis,'position');
imageXRatio = dim(2)/dim(1);
imageYRatio = dim(1)/dim(2);
if(imageXRatio>1)
    imageXRatio = 1;
else
    imageYRatio = 1;
end
imageAxisPos = [imageAxisPos(1) imageAxisPos(2) imageAxisPos(3)*imageYRatio imageAxisPos(4)*imageXRatio];
set(ImageAxis, 'position', imageAxisPos);


%use voxel dimensions of image to keep aspect ratio
voxYAxisRatio = VOX(1)/VOX(3);
voxXAxisRatio = VOX(3)/VOX(1);
if(voxYAxisRatio > 1)
    voxYAxisRatio =1;
else
    voxXAxisRatio =1;
end
R=ASPECT_RATIO_FOR_SQUARE_FIGURE;
imageAxisPos = [imageAxisPos(1) imageAxisPos(2) imageAxisPos(3)*voxYAxisRatio imageAxisPos(4)*voxXAxisRatio];
set(ImageAxis, 'position', imageAxisPos);

%keeping the same aspect ratio as calculated above enlarge image to fit
% in to its allocated space
enlargeScalar = min([axQuad1Pos(3) axQuad1Pos(4)])/max([imageAxisPos(3) imageAxisPos(4)]);
imageAxisPos = [imageAxisPos(1) imageAxisPos(2) imageAxisPos(3)*enlargeScalar imageAxisPos(4)*enlargeScalar];
set(ImageAxis,'position',imageAxisPos);
colormap(newColor);
quad(1).xCord =[ imageAxisPos(1), imageAxisPos(1)+imageAxisPos(3),imageAxisPos(1)+imageAxisPos(3),imageAxisPos(1)];
quad(1).yCord =[ imageAxisPos(2)+imageAxisPos(4), imageAxisPos(2)+imageAxisPos(4),imageAxisPos(2),imageAxisPos(2)];
set(ImageAxis,'XTickLabel',[]);
set(ImageAxis,'YTickLabel',[]);

%--------------------------------------------------------
%Draw and scale image in Quad 2


%use dimensions of y and x axis to keep aspect ratio
% ie makes axis square within allocated space
axQuad2Pos = [.5 .70 .45 .25];
dim = size(sliceYZ);
ImageAxis = axes('Parent', handles, 'units', 'normalized', 'position', axQuad2Pos, 'color', [1 1 1]);
axis(ImageAxis, 'off');
%origData.axisHandle1=imagesc(rot90(sliceYZ),[minInterval 2*maxInterval]);
% using image as it is more general than the imagesc
origData.axisHandle1 = image(rot90(sliceYZ), 'parent', ImageAxis, 'CDataMapping', 'scaled');
set(ImageAxis, 'clim', [minInterval 2*maxInterval]); % set the axis clim property
imageAxisPos = get(ImageAxis, 'position');
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
set(ImageAxis,'position',imageAxisPos);
axQuad2Pos = imageAxisPos;

%use actual dimensions of image to keep aspect ratio
imageAxisPos = get(ImageAxis, 'position');
imageXRatio = dim(2)/dim(1);
imageYRatio = dim(1)/dim(2);
if(imageXRatio>1)
    imageXRatio = 1;
else
    imageYRatio = 1;
end
imageAxisPos = [imageAxisPos(1) imageAxisPos(2) imageAxisPos(3)*imageYRatio imageAxisPos(4)*imageXRatio];
set(ImageAxis,'position',imageAxisPos);

%use voxel dimensions of image to keep aspect ratio
voxYAxisRatio = VOX(2)/VOX(3);
voxXAxisRatio = VOX(3)/VOX(2);
if(voxYAxisRatio > 1)
    voxYAxisRatio =1;
else
    voxXAxisRatio =1;
end
R=ASPECT_RATIO_FOR_SQUARE_FIGURE;
imageAxisPos = [imageAxisPos(1) imageAxisPos(2) imageAxisPos(3)*voxYAxisRatio imageAxisPos(4)*voxXAxisRatio];
set(ImageAxis,'position',imageAxisPos);

%keeping the same aspect ratio as calculated above enlarge image to fit
% in to its allocated space
enlargeScalar = min([axQuad2Pos(3) axQuad2Pos(4)])/max([imageAxisPos(3) imageAxisPos(4)]);
imageAxisPos = [imageAxisPos(1) imageAxisPos(2) imageAxisPos(3)*enlargeScalar imageAxisPos(4)*enlargeScalar];
set(ImageAxis,'position',imageAxisPos);
colormap(newColor);
quad(2).xCord =[ imageAxisPos(1), imageAxisPos(1)+imageAxisPos(3),imageAxisPos(1)+imageAxisPos(3),imageAxisPos(1)];
quad(2).yCord =[ imageAxisPos(2)+imageAxisPos(4), imageAxisPos(2)+imageAxisPos(4),imageAxisPos(2),imageAxisPos(2)];
set(ImageAxis,'XTickLabel',[]);
set(ImageAxis,'YTickLabel',[]);


%Quad 3

%use dimensions of y and x axis to keep aspect ratio
% ie makes axis square within allocated space
axQuad3Pos = [.05 .4 .45 .25];
dim = size(sliceXY);
ImageAxis = axes('Parent', handles, 'units','normalized','position',axQuad3Pos,'color',[1 1 1]);
axis(ImageAxis, 'off');
% origData.axisHandle1=imagesc(rot90(sliceXY),[minInterval 2*maxInterval]);
% use image as it is more general than imagesc
origData.axisHandle1 = image(rot90(sliceXY), 'parent', ImageAxis, 'CDataMapping', 'scaled');
set(ImageAxis, 'clim', [minInterval 2*maxInterval]); % set the axis clim property
imageAxisPos = get(ImageAxis,'position');
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
set(ImageAxis,'position',imageAxisPos);
axQuadPos = imageAxisPos;

%use actual dimensions of image to keep aspect ratio
imageAxisPos = get(ImageAxis,'position');
imageXRatio = dim(2)/dim(1);
imageYRatio = dim(1)/dim(2);
if(imageXRatio>1)
    imageXRatio = 1;
else
    imageYRatio = 1;
end
imageAxisPos = [imageAxisPos(1) imageAxisPos(2) imageAxisPos(3)*imageYRatio imageAxisPos(4)*imageXRatio];
set(ImageAxis,'position',imageAxisPos);


%use voxel dimensions of image to keep aspect ratio
voxYAxisRatio = VOX(2)/VOX(1);
voxXAxisRatio = VOX(1)/VOX(2);
if(voxYAxisRatio > 1)
    voxYAxisRatio =1;
else
    voxXAxisRatio =1;
end
R=ASPECT_RATIO_FOR_SQUARE_FIGURE;
imageAxisPos = [imageAxisPos(1) imageAxisPos(2) imageAxisPos(3)*voxYAxisRatio imageAxisPos(4)*voxXAxisRatio];
set(ImageAxis,'position',imageAxisPos);

%keeping the same aspect ratio as calculated above enlarge image to fit
% in to its allocated space
enlargeScalar = min([axQuad1Pos(3) axQuad1Pos(4)])/max([imageAxisPos(3) imageAxisPos(4)]);
imageAxisPos = [imageAxisPos(1) imageAxisPos(2) imageAxisPos(3)*enlargeScalar imageAxisPos(4)*enlargeScalar];
set(ImageAxis,'position',imageAxisPos);
colormap(newColor);
quad(3).xCord =[ imageAxisPos(1), imageAxisPos(1)+imageAxisPos(3),imageAxisPos(1)+imageAxisPos(3),imageAxisPos(1)];
quad(3).yCord =[ imageAxisPos(2)+imageAxisPos(4), imageAxisPos(2)+imageAxisPos(4),imageAxisPos(2),imageAxisPos(2)];
set(ImageAxis,'XTickLabel',[]);
set(ImageAxis,'YTickLabel',[]);

%colorbar
colorPos = [0.85 0.4 .045 .25];
drawnow;
%cbAxis = axes('Parent', handles, 'units','normalized', 'position', colorbarPos, 'color', [1 1 1]);
ColorBarHandle = colorbar; %(cbAxis);
set(ColorBarHandle, 'position', colorPos);
%colorPos = get(ColorBarHandle, 'position');
%axis off;

childH = get(ColorBarHandle,'Children');
%set(childH, 'YData', [minInterval 2*maxInterval]);
imInd = strmatch('image', lower(get(childH, 'Type')), 'exact');
set(childH(imInd), 'YData', [minInterval 2*maxInterval]);
set(ColorBarHandle, 'YLim', [minInterval maxInterval]);
set(ColorBarHandle, 'YColor', FONT_COLOR);
set(ColorBarHandle, 'YTick', [minInterval maxInterval]);
icaCLIM = round(icaCLIM*10)/10;
set(ColorBarHandle, 'YTickLabel', [icaCLIM(1) icaCLIM(2)]);
%colorPos=get(ColorBarHandle,'position');
set(ColorBarHandle, 'position', [colorPos(1)-.04,colorPos(2)+.03,colorPos(3),colorPos(4)]);
%     title('z score values','color',FONT_COLOR);


%setup origData timecourse axes and plot BOLD signal
origData.S4Handle = subplot(3, 2, 5);
cL = COLORLIST;
cLSize = size(cL,2);

%     for i=1:size(BOLDTimecourse,1)
%         c = cL(mod(i,cLSize));
%   end
if(size(BOLDTimecourse,1)==1)
    origData.axisHandle4 = plot(BOLDTimecourse(1,:), 'r', 'parent', origData.S4Handle);
else
    origData.axisHandle4 = plot(BOLDTimecourse(1,:), 'r', 'parent', origData.S4Handle);
    hold on;
    origData.axisHandle4 = plot(BOLDTimecourse(2,:), 'y:', 'parent', origData.S4Handle);
    hold on;
    origData.axisHandle4 = plot(BOLDTimecourse(3,:), 'y:', 'parent', origData.S4Handle);
end
hold off;
axis(origData.S4Handle, 'tight');
set(origData.S4Handle, 'YColor', [1 1 1])
set(origData.S4Handle, 'XColor', [1 1 1])
title(BOLDTimecourseLabel, 'parent', origData.S4Handle);
if (strcmpi(icatb_get_modality, 'fmri'))
    ylabel('BOLD Signal', 'color', [1 1 1], 'parent', origData.S4Handle);
else
    ylabel('Original Data', 'color', [1 1 1], 'parent', origData.S4Handle);
end



%setup ica timecourse axes and plot ica timecourse
origData.S5Handle = subplot(3,2,6);
detrendA = icatb_detrend(A, 1, length(A), DETRENDNUMBER);
detrendMinICA = icatb_detrend(minICATC, 1, length(minICATC), DETRENDNUMBER);

plot((1:length(detrendA)), detrendA, 'r', (1:length(detrendA)), detrendMinICA, 'r:', ...
    'parent', origData.S5Handle);

axis(origData.S5Handle, 'tight');
%label = ['ICA Timecourse(', [num2str(realWorldPos(1)), ',', num2str(realWorldPos(2)), ',', num2str(realWorldPos(3))], ')'];
if (strcmpi(icatb_get_modality, 'fmri'))
    label = ['ICA Timecourse(', [num2str(realWorldPos(1)), ',', num2str(realWorldPos(2)), ',', ...
        num2str(realWorldPos(3))], ')', ' (DETRENDNUMBER = ', num2str(DETRENDNUMBER), ')'];
else
    label = ['ICA Loadings (', [num2str(realWorldPos(1)), ',', num2str(realWorldPos(2)), ',', num2str(realWorldPos(3))], ')'];
end
title(label, 'parent', origData.S5Handle);
clear label;
if (strcmpi(icatb_get_modality, 'fmri'))
    ylabel('Signal Units', 'color', [1 1 1], 'parent', origData.S5Handle);
    xlabel('Scans', 'color', [1 1 1], 'parent', origData.S5Handle);
else
    ylabel('ICA Loadings', 'color', [1 1 1], 'parent', origData.S5Handle);
    xlabel('Subjects', 'color', [1 1 1], 'parent', origData.S5Handle);
end

% move BOLD timecourse position
plotPos = get(get(origData.axisHandle4, 'Parent'), 'position');
plotPos = [plotPos(1),plotPos(2)+.03,(plotPos(4)+.05)*2+.1,(plotPos(4)+.05)/2.4];
set(get(origData.axisHandle4, 'Parent'), 'position', [plotPos(1),plotPos(2)+.1,plotPos(3), plotPos(4)] );

%BOLD timecourse popup
pos = [plotPos(1)+plotPos(3)+.01,plotPos(2)+.19,.18,.05];
label = 'Plot fMRI data';
if (strcmpi(icatb_get_modality, 'smri'))
    label = 'Plot sMRI data';
end
timePlotLabel = icatb_uicontrol('parent', handles, 'style', 'text', 'units', 'normalized', 'string', label, ...
    'position', pos, 'visible', 'on', 'tag', 'fmriLabel', 'fontunits', UI_FONTUNITS, 'fontname', UI_FONTNAME, ...
    'FontSize', UI_FS - 2, 'ForegroundColor', FONT_COLOR, 'BackgroundColor', BG2_COLOR);

[outString, newPos] = textwrap(timePlotLabel, {get(timePlotLabel, 'string')});
pos(4) = newPos(4);
set(timePlotLabel, 'position', pos, 'string', outString);
%timePlotLabel = icatb_getUILabel(handles, label, pos, 'on', 'fmriLabel');

label = ['At voxel ', num2str(voxelCoord)];

pos = [pos(1),pos(2)-.025,pos(3),pos(4)];
%timePlotLabel = icatb_getUILabel(handles, label, pos, 'on', 'VoxelLabel');
timePlotLabel = icatb_uicontrol('parent', handles, 'style', 'text', 'units', 'normalized', 'string', label, ...
    'position', pos, 'visible', 'on', 'tag', 'VoxelLabel', 'fontunits', UI_FONTUNITS, 'fontname', UI_FONTNAME, ...
    'FontSize', UI_FS - 2, 'ForegroundColor', FONT_COLOR, 'BackgroundColor', BG2_COLOR);

[outString, newPos] = textwrap(timePlotLabel, {get(timePlotLabel, 'string')});
pos(4) = newPos(4);
set(timePlotLabel, 'position', pos, 'string', outString);
pos=[plotPos(1)+plotPos(3)+.04,plotPos(2)+.1,.15,.05];
if (~strcmpi(icatb_get_modality, 'smri'))
    for i=1:size(fmriFiles,2)
        for k=1:size(fmriFiles(1).ses,2)
            if(i==1 & k==1)
                labels = str2mat(['Subject ',num2str(i),' Session ',num2str(k)]);
            else
                labels = str2mat(labels,['Subject ',num2str(i),' Session ',num2str(k)]);
            end
        end
    end
else
    labels = 'All Subjects';
end

if numOfSub*numOfSess > 1
    labels = str2mat(labels,'All(Mean & SEM)');
end

popupH = findobj(handles, 'tag', 'fMRIPopUp');
if isempty(popupH)
    data.timePlotPopUp = icatb_uicontrol('parent', handles, 'units', 'normalized', 'style', 'popup', 'position', ...
        pos, 'string', labels, 'BackgroundColor', BG2_COLOR, 'ForegroundColor', FONT_COLOR, 'fontunits', ...
        UI_FONTUNITS, 'fontname', UI_FONTNAME, 'fontsize', UI_FS - 2, 'tag', 'fMRIPopUp', 'callback', ...
        {@fMRITimecourseCallback, handles});
end
%plot scaled ica timecourse on top of timecourse
set(origData.S5Handle, 'position', [plotPos(1),plotPos(2)-.08,plotPos(3),plotPos(4)]);
hold on;
detrendScaledVal = icatb_detrend(scaledVal, 1, length(scaledVal), DETRENDNUMBER);
plot((1:length(detrendScaledVal)), detrendScaledVal, 'g', 'parent', origData.S5Handle);
set(origData.S5Handle, 'YColor', [1 1 1])
set(origData.S5Handle, 'XColor', [1 1 1])
hold off;
%axis([0 size(A,1) minA maxA]);
axis(origData.S4Handle, 'tight');
axis(origData.S5Handle, 'tight');
if (~strcmpi(icatb_get_modality, 'smri'))
    icatb_legend(origData.S5Handle, 'ICA TC at Max Voxel', 'ICA TC at Min Voxel', 'IC TC at Current Voxel');
else
    icatb_legend(origData.S5Handle, 'ICA Loadings at Max Voxel', 'ICA Loadings at Min Voxel', 'IC Loadings at Current Voxel');
end
% legend2H = legend(origData.S5Handle, 'ICA TC at Max Voxel', 'ICA TC at Min Voxel', 'IC TC at Current Voxel', 0);
% if whichVersion > 13
%     set(legend2H, 'textColor', LEGENDCOLOR);
% end
%'Voxel Fit To ICA TC',0)


%setup top  ica timecourse button
pos = [plotPos(1)+plotPos(3)+.01,plotPos(2),.2,.02];
%timePlotLabel = icatb_getUILabel(handles, 'Top 5 component timecourses', pos, 'on', 'fmriPlot');
if (~strcmpi(icatb_get_modality, 'smri'))
    timeStr = 'Top 5 component timecourses';
else
    timeStr = 'Top 5 component loadings';
end
timePlotLabel = icatb_uicontrol('parent', handles, 'style', 'text', 'units', 'normalized', 'string', ...
    timeStr, 'position', pos, 'visible', 'on', 'tag', 'fmriPlot', 'fontunits', ...
    UI_FONTUNITS, 'fontname', UI_FONTNAME, 'FontSize', UI_FS - 2, 'ForegroundColor', FONT_COLOR, ...
    'BackgroundColor', BG2_COLOR);
[outString, newPos] = textwrap(timePlotLabel, {get(timePlotLabel, 'string')});
pos(4) = newPos(4);
set(timePlotLabel, 'position', pos, 'string', outString);

plotTag = findobj(handles, 'Tag', 'topTimePlotButton');
if isempty(plotTag)
    topTimePlotButton = icatb_uicontrol('parent', handles, 'Tag', 'topTimePlotButton', 'string', 'Plot',...
        'units', 'normalized', 'position', [plotPos(1)+plotPos(3)+.08,plotPos(2)-.08,.05,.05],...
        'style','pushbutton', 'Visible', 'on', ...
        'backgroundcolor', BUTTON_COLOR, 'foregroundcolor', BUTTON_FONT_COLOR, 'fontunits', UI_FONTUNITS, ...
        'fontname', UI_FONTNAME, 'FontSize', UI_FS);
    extentPos = get(topTimePlotButton, 'Extent');
    extentPos(3) = extentPos(3) + 0.01; extentPos(4) = extentPos(4) + 0.005;
    plotButtonPos = get(topTimePlotButton, 'position');
    plotButtonPos(3:4) = extentPos(3:4);
    set(topTimePlotButton, 'position', plotButtonPos);
    set(topTimePlotButton, 'callback', {@plotTopTimecoursesCallback, handles});
end

%save data in figure
data.quad = quad;
set(handles, 'userdata', data);
clear data;

% set up mouseClick listener
set(handles, 'WindowButtonDownFcn', @getPointsCallback);

%%%%%%%%%% Function callbacks %%%%%%%%%%%%

% using function callbacks
function getPointsCallback(hObject, event_data)
% get points
% purpose:
%get handles to each figure
data = get(hObject, 'userdata');
figureName = get(hObject, 'name');
quad = data.quad;
structHInfo = data.structHInfo;
imDIM = data.imDIM;
fmriHInfo = data.fmriHInfo;
realPos = data.realPos;
icaPos = data.icaPos;
imOrigin = data.imOrigin;
set(hObject, 'Units', 'normalized');


p = get(hObject, 'currentpoint');

%find out if the click was in a quad
isQuad1= (p(1) > quad(1).xCord(1)) & (p(1) < quad(1).xCord(2)) & (p(2) < quad(1).yCord(1)) & (p(2) >quad(1).yCord(3));
isQuad2= (p(1) > quad(2).xCord(1)) & (p(1) < quad(2).xCord(2)) & (p(2) < quad(2).yCord(1)) & (p(2) >quad(2).yCord(3));
isQuad3= (p(1) > quad(3).xCord(1)) & (p(1) < quad(3).xCord(2)) & (p(2) < quad(3).yCord(1)) & (p(2) >quad(3).yCord(3));


if(~isQuad1 & ~isQuad2 & ~isQuad3)
    newX=icaPos(1);
    newY=icaPos(2);
    newZ=icaPos(3);
else
    %get points from mouse click and decide which quad they clicked in
    %quad 1 is upper left, quad 2 is upper right, quad 3 is lower
    %left
    if(isQuad1)
        %quad 1
        
        %convert units(from mouse click) to image dimensions
        [x y]=ginput(1);
        newX = ceil(x);%round(p2(1)*DIM(1)/maxP2(1));
        newZ = ceil(imDIM(3)-y);%round(p2(2)*DIM(3)/maxP2(2));
        newY = icaPos(2);
        
    end
    if(isQuad2)
        %quad 2
        %Y DIM
        %convert units(from mouse click) to image dimensions
        [x y]=ginput(1);
        
        newY= ceil(x);%round((p2(1)*DIM(2))/maxP2(1));
        newZ = ceil(imDIM(3)-y);%round(p2(2)*DIM(3)/maxP2(2));
        newX = icaPos(1);
        
        
    end
    
    if(isQuad3)
        %quad 3
        
        %             %X DIM
        %convert units(from mouse click) to image dimensions
        [x y] = ginput(1);
        
        newX = ceil(x);%round(p2(1)*DIM(1)/maxP2(1));
        newY = ceil(imDIM(2)-y);%round(p2(2)*DIM(2)/maxP2(2));
        newZ = icaPos(3);
        
    end
    
    icaPos(1)=newX;
    icaPos(2)=newY;
    icaPos(3)=newZ;
    
    
    %     %get the real positiosn
    realPos(1) = [round(icaPos(1)*(fmriHInfo.DIM(1)/imDIM(1)))];
    realPos(2) = [round(icaPos(2)*(fmriHInfo.DIM(2)/imDIM(2)))];
    realPos(3) = [round(icaPos(3)*(fmriHInfo.DIM(3)/imDIM(3)))];
    
    
    realPos = icaPos;
    
    
    %save data
    data.realPos = realPos;
    data.icaPos = icaPos;
    set(hObject, 'userdata', data);
    clear data;
    redrawCallback(hObject);
    if (~strcmpi(icatb_get_modality, 'smri'))
        disp('If you want to look at the fMRI time course for the current voxel, select the subject under popup box.');
    else
        disp('If you want to look at the ICA loadings for the current voxel, select the all subjects under popup box.');
    end
    
    fprintf('\n');
    %popUpHandle = findobj(hObject, 'tag', 'fMRIPopUp');
    % update the bold plot
    %fMRITimecourseCallback(popUpHandle, [], hObject);
end


function fMRITimecourseCallback(hObject, event_data, handles)

set(handles, 'pointer', 'watch');

icatb_defaults;
global DETRENDNUMBER;

drawnow;

data = get(handles, 'userdata');

real_world_coords = data.real_world_coords;
fmriFiles = data.fmriFiles;
realPos = data.realPos;
icaPos = data.icaPos;
complexInfoRead = data.complexInfoRead;
dataType = data.dataType;
diffTimePoints = data.diffTimePoints;
flagTimePoints = data.flagTimePoints;
fmriHInfo = data.fmriHInfo;

popUpValue = get(hObject, 'value');
numOfSub = size(fmriFiles,2);
numOfSess = size(fmriFiles(1).ses,2);
answerQuestion = 1;


[realWorldPos, voxelCoord] = getCoord(fmriHInfo.V(1), real_world_coords, realPos);

%check to make sure they want to look at mean
if(popUpValue > numOfSub*numOfSess)
    string1 = ['Do you want to look at the mean and standard error of the mean for all ', num2str(numOfSub*numOfSess), ' data sets.'];
    %answer = questdlg(string);
    [answerQuestion] = icatb_questionDialog('title', 'Calculate Mean', 'textbody', string1);
    
end

if answerQuestion == 1
    %--Mean and Standard Error of Mean for all Subject's Timecourse at a
    %--particular voxel
    if(popUpValue > numOfSub*numOfSess)
        %timecourse=zeros(numOfSub,size(fmriFiles(1).ses(1).name,1));
        tp =  min(diffTimePoints);
        timecourse = zeros(numOfSub, tp);
        
        %         if strcmpi(flagTimePoints, 'different_time_points')
        %             disp('Not calculating mean as functional data has different time dimensions over subjects');
        %             return;
        %         end
        
        icatb_statusBar('init',numOfSub*numOfSess,'Calculating Mean','','');
        counter =1;
        for i = 1:numOfSub
            for k=1:numOfSess
                P=fmriFiles(i).ses(k).name;
                tmpTC = icatb_getVoxelValues(realPos, P, data.structFile, ...
                    dataType, complexInfoRead, 'read');
                timecourse(counter,:) = tmpTC(1:tp);
                %timecourse(counter,:) = 100*detrend(timecourse(counter,:),0)/mean(timecourse(counter,:));
                timecourse(counter,:) = 100*(icatb_detrend(timecourse(counter, :), 1, length(timecourse(counter,:)), ...
                    DETRENDNUMBER)') / mean(timecourse(counter, :));
                counter=counter+1;
                icatb_statusBar('increment',1);
                disp(['Subject ',num2str(i),' Session ',num2str(k), ' timecourse loaded']);
            end
        end
        meanTimecourse = mean(timecourse);
        SEMTimecourse = std(timecourse);
        SEMTimecourse = SEMTimecourse./sqrt(numOfSub*numOfSess - 1);
        data.BOLDTimecourse = [meanTimecourse; SEMTimecourse + meanTimecourse; meanTimecourse - SEMTimecourse];
        data.BOLDTimecourseLabel = ['fMRI Data for all Subjects at Voxel ',num2str(voxelCoord),'(Mean & SEM)'];
        icatb_statusBar('finished');
    else
        %--Single subject's timecourse at a particular voxel
        
        subNum = ceil(popUpValue/numOfSess);
        sesNum = mod(popUpValue-1,numOfSess)+1;
        P=fmriFiles(subNum).ses(sesNum).name;
        
        numPlots = size(data.BOLDTimecourse,1);
        timecourse = icatb_getVoxelValues(realPos, P, data.structFile, dataType, complexInfoRead, 'read');
        data.BOLDTimecourse = timecourse;
        if (~strcmpi(icatb_get_modality, 'smri'))
            data.BOLDTimecourseLabel = ['fMRI Data for sub ',num2str(subNum), ' sess ', num2str(sesNum), ' at Voxel ', num2str(voxelCoord)];
        else
            data.BOLDTimecourseLabel = ['sMRI Data at Voxel ', num2str(voxelCoord)];
        end
        
    end
    
    set(handles, 'userdata', data);
    clear data;
    
    redrawCallback(handles);
end

clear data;

set(handles, 'pointer', 'arrow');

function plotTopTimecoursesCallback(hObject, event_data, handles)
% get points
% purpose:

icatb_defaults;
global COLORLIST
global LEGENDCOLOR;
global DETRENDNUMBER;
global SMOOTHPARA;
global SMOOTHINGVALUE;

% get data
data = get(handles, 'userdata');
icaFiles = data.icaFiles;
realPos = data.realPos;
structFile = data.structFile;
dataType = data.dataType;
complexInfo = data.complexInfoWrite;
zipContents = data.zipContents;
fmriHInfo = data.fmriHInfo;
real_world_coords = data.real_world_coords;

zipFileName = [];
files_in_zip = {};

[realWorldPos, voxelCoord] = getCoord(fmriHInfo.V(1), real_world_coords, realPos);

% get the files from the zip file
[zipFileName, files_in_zip] = icatb_getViewingSet_zip(icaFiles, complexInfo, dataType, zipContents);

icaFiles = char(strrep(cellstr(icaFiles), data.outputDir, ''));

% full file path
icaFiles = icatb_fullFile('directory', data.outputDir, 'files', icaFiles);

% Unzip files
if ~isempty(zipFileName)
    icatb_unzip(regexprep(zipFileName, ['.*\', filesep], ''), fileparts(icaFiles(1, :)));
end

% maximum number of components
maxComp = icatb_get_countTimePoints(icaFiles);


if maxComp >= 5
    if (strcmpi(icatb_get_modality, 'fmri'))
        figLabel = 'Components Making Up Majority of BOLD Signal';
    else
        figLabel = 'Top 5 Components';
    end
    Tag = 'TopComponents';
    [GraphicsHandle] = icatb_getGraphics(figLabel, 'timecourse', Tag);
    [A, values, topComponents] = icatb_getTopTimecourses(realPos, icaFiles, structFile, dataType, complexInfo);
    
    % smooth ica time course
    if strcmpi(SMOOTHPARA, 'yes')
        A = icatb_gauss_smooth1D(A, SMOOTHINGVALUE);
    end
    
    % detrend ICA time course
    A = icatb_detrend(A, 1, size(A, 1), DETRENDNUMBER);
    
    disp(['Showing top 5 components at voxel location ', num2str(voxelCoord)]);
    for i=5:-1:1
        disp(['Component ', num2str(topComponents(i)),' value ',num2str(values(i))]);
    end
    
    figure(GraphicsHandle);
    axisH = axes('parent', GraphicsHandle);
    
    % plot top five components
    values = abs(values);
    for ii = 1:5
        plot(values(ii)*A(:, ii), COLORLIST(ii), 'parent', axisH);
        hold on;
    end
    axis(axisH, 'tight');
    icatb_legend(['Component ',num2str(topComponents(1))],['Component ',num2str(topComponents(2))],...
        ['Component ',num2str(topComponents(3))],['Component ',num2str(topComponents(4))],['Component ',num2str(topComponents(5))]);
    hold off;
    if (strcmpi(icatb_get_modality, 'fmri'))
        ylabel('Signal Units', 'parent', axisH);
        xlabel('Scans', 'parent', axisH);
    else
        ylabel('ICA Loadings', 'parent', axisH);
        xlabel('Subjects', 'parent', axisH);
    end
    % show in the title
    title(['Top five components (DETRENDNUMBER = ', num2str(DETRENDNUMBER), ')'], 'parent', axisH);
    set(axisH, 'YColor', [1 1 1]);
    set(axisH, 'XColor', [1 1 1]);
    
else
    disp('No. of component files should be greater than or equal to five');
end

% % zip the images back
if ~isempty(zipFileName)
    fileN = str2mat(files_in_zip);
    %fileN = icatb_fullFile('directory', data.outputDir, 'files', fileN);
    icatb_delete_file_pattern(fileN, data.outputDir);
end

function setVoxelPosCallback(hObject, event_data, handles)
% set voxel position

data = get(handles, 'userdata');
real_world_coords = data.real_world_coords;

voxelOrigin = data.voxelOrigin;
fmriHInfo = data.fmriHInfo;

% open input dialog box
prompt = {'Enter voxel position in real world coordinates (mm)'};
dlg_title = 'Set voxel position';
num_lines = 1;
[realWorldPos, voxelCoord] = getCoord(fmriHInfo.V(1), real_world_coords, data.realPos);
def = {num2str(realWorldPos)};

% get voxel coords
voxelPos = icatb_inputdlg2(prompt, dlg_title, num_lines, def);

drawnow;

% set the data
if ~isempty(voxelPos)
    voxelPos = str2num(voxelPos{1});
    %voxelPos = round(voxelPos);
    
    if length(voxelPos) ~= 3
        error('Please enter the three coords correctly');
    end
    
    voxelPos = round(voxelOrigin + (voxelPos./fmriHInfo.VOX));
    
    %checkIndex = find(voxelPos == data.realPos);
    
    
    % check if the position is the same
    %if length(checkIndex) ~= 3
    data.realPos = voxelPos;
    data.icaPos = voxelPos;
    
    set(handles, 'userdata', data);
    popupH = findobj(handles, 'tag', 'fMRIPopUp');
    % plot bold time course
    %fMRITimecourseCallback(popupH, [], handles);
    redrawCallback(handles);
    if (~strcmpi(icatb_get_modality, 'smri'))
        disp('If you want to look at the fMRI time course for the current voxel, select the subject under popup box.');
    else
        disp('If you want to look at the ICA loadings for the current voxel, select the all subjects under popup box.');
    end
    fprintf('\n');
    %end
end


function figCloseCallback(hObject, event_data, handles)
% figure close callback

data = get(hObject, 'userdata');

files_in_zip = data.files_to_be_deleted;


% delete the unzipped files
if ~isempty(files_in_zip)
    icatb_delete_file_pattern(files_in_zip);
end

delete(hObject);


function [realWorldPos, voxelCoord] = getCoord(V, real_world_coords, pixelPos)
% Get real world and voxel coordinates

realWorldPos = squeeze(real_world_coords(pixelPos(1), pixelPos(2), pixelPos(3), :));
realWorldPos = realWorldPos(:)';
voxelCoord = icatb_real_to_voxel(V, realWorldPos);