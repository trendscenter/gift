function varargout = icatb_orth_views(file, varargin)
%% Orthogonal viewer
%
% Inputs:
% 1. file - Input file
% 2. varargin - Variable no. of arguments
%
% Outputs:
% slices
image_values = 'positive and negative';
cmap = [];
threshold = 0;
convert_to_zscores = 'no';
coords = [0, 0, 0];
structFile = '';
isInteractive = 1;
set_to_max_voxel = 0;
fig_title = 'OrthoViews';
returnInterpdata = 0;


icatb_defaults;
global UI_FS;
global UI_FONTNAME;

%% Loop over variable no. of arguments
for n = 1:2:length(varargin)
    if (strcmpi(varargin{n}, 'structfile'))
        structFile = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'coords'))
        coords = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'image_values'))
        image_values = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'cmap'))
        cmap = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'threshold'))
        threshold = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'convert_to_zscores'))
        convert_to_zscores = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'labels'))
        labels = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'isinteractive'))
        isInteractive = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'fig_title'))
        fig_title = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'tc'))
        tc = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'set_to_max_voxel'))
        set_to_max_voxel = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'get_interp_data'))
        returnInterpdata = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'position'))
        fig_pos = varargin{n + 1};
    end
end

if (length(coords) == numel(coords))
    coords = coords(:)';
end

if (size(coords, 2) ~= 3)
    error('No. of columns of coords is not equal to 3. coords must contain x,y,z coordsinates');
end

% Turn off interactive if multiple coordsinates are plotted
if (length(coords) ~= numel(coords))
    isInteractive = 0;
end

if (~exist('labels', 'var'))
    labels = repmat({''}, 1, size(coords, 1));
end

labels = cellstr(labels);
convertToZ = strcmpi(convert_to_zscores, 'yes');

returnValue = strmatch(lower(image_values), {'positive and negative', 'positive', 'absolute value', 'negative'}, 'exact');

%% Interpolate selected image to structural dimensions. Get the interactive
%% input?
[pathstr, file, extn] = fileparts(file);
if isempty(pathstr)
    pathstr = pwd;
end

file = fullfile(pathstr, [file, extn]);

if (~isempty(structFile))
    if (ischar(structFile))
        [pathstr, structFile, extn] = fileparts(structFile);
        if isempty(pathstr)
            pathstr = pwd;
        end
        structFile = fullfile(pathstr, [structFile, extn]);
    end
else
    structFile = file;
end

[structFile, n] = icatb_parseExtn(structFile);
structVol = icatb_spm_vol_nifti(structFile, n);
structVol = structVol(1);

[images, real_world_coords] = icatb_resizeImage(structVol, file, 'axial', [], 1);
%data = reshape(images(2, :, :, :), structVol.dim(1:3));
%structData = reshape(images(1, :, :, :), structVol.dim(1:3));
data = reshape(images(2, :, :, :), [size(images, 2), size(images, 3), size(images, 4)]);
structData = reshape(images(1, :, :, :), [size(images, 2), size(images, 3), size(images, 4)]);
clear images;

stdData = 1;
meanData = 0;
if (convertToZ)
    meanData = mean(data(data~=0));
    stdData = std(data(data~=0));
    data(data~=0) =  (data(data ~= 0) - meanData)/stdData;
end

data = applyDispParams(data, returnValue, threshold);

good_inds = find(abs(data) > eps);

minICAIM = min(data(good_inds));
maxICAIM = max(data(good_inds));

maxICAIM = max([maxICAIM, abs(threshold)]);

if (isempty(minICAIM))
    minICAIM = eps;
end

if (isempty(maxICAIM))
    maxICAIM = eps;
end

if (returnValue == 1)    
    maxICAIM = max(abs([minICAIM, maxICAIM]));
    minICAIM = -maxICAIM;
end

dataRange = [minICAIM, maxICAIM];

[dd, inds] = max(abs(data(:)));
[x, y, z] = ind2sub(size(data), inds);
maxVoxelPos = squeeze(real_world_coords(x, y, z, :))';
if (~set_to_max_voxel)
    pixelPos = round(icatb_real_to_voxel(structVol, coords(1, :)));
else
    pixelPos = [x, y, z];
end


maxICAIM = round(maxICAIM*10)/10;
minICAIM = round(minICAIM*10)/10;

minInterval = 0;
maxInterval = 100;

%% Get colormap associated with the image values
if (isempty(cmap))
    cmap = icatb_getColormap(1, returnValue, 1);
    if (returnValue == 1)
        colorLength = size(cmap, 1)/2;
        load icatb_colors coldhot;
        nToSkip = ceil(size(coldhot, 1)/colorLength);
        cmap(1:colorLength, :) = coldhot(1:nToSkip:end, :);
    end
else
    cmap = [cmap; gray(size(cmap, 1))];
end

%composite_data = reshape(composite_data, structVol.dim(1:3));

% % Convert to voxel coordinates
%for nC = 1:size(coords, 1)
%coords(nC, :) = round(icatb_real_to_voxel(structVol, coords(nC, :)));
%end

CLIM = [minInterval, 2*maxInterval];

if (returnInterpdata)   
    composite_data = getCompositeData(structData, data, returnValue, minInterval, maxInterval, dataRange);
    handles_data.data = composite_data;
else
    handles_data.dispPrefs = struct('returnValue', returnValue, 'stdData', stdData, 'meanData', meanData, 'imgVol', ...
        icatb_returnHInfo(file), 'threshold', threshold);
end

handles_data.CLIM = CLIM;
handles_data.range = dataRange;
handles_data.minInterval = minInterval;
handles_data.maxInterval = maxInterval;
handles_data.pixelPos = pixelPos;
handles_data.coords = coords;
handles_data.origin = [0, 0, 0];
handles_data.voxelOrigin = icatb_real_to_voxel(structVol(1), handles_data.origin);
handles_data.VOX = abs([structVol(1).mat(1, 1), structVol(1).mat(2, 2), structVol(1).mat(3, 3)]);
handles_data.structVol = structVol;
handles_data.minICAIM = minICAIM;
handles_data.maxICAIM = maxICAIM;
handles_data.labels = labels;
handles_data.maxVoxelPos = maxVoxelPos;
handles_data.isInteractive = isInteractive;
handles_data.cmap = cmap;

if (returnInterpdata)
    handles_data.maxVoxelPos = handles_data.pixelPos;
    varargout{1} = handles_data;
    return;
end

if exist('tc', 'var')
    handles_data.numRowsToPlot = (size(coords, 1) + 1);
else
    handles_data.numRowsToPlot = (size(coords, 1));
end

%% Figure and Menus
gH = icatb_getGraphics(fig_title, 'graphics', 'orthoviews', 'off');

if (exist('fig_pos', 'var'))
    set(gH, 'position', fig_pos);
end

if (isInteractive)
    cmenu = uicontextmenu;
    set(gH, 'UIContextMenu', cmenu);
    uimenu(cmenu, 'Label', 'Set Voxel Position (mm)', 'Callback', {@setVoxelPos, gH});
    maxVoxelH = uimenu(cmenu, 'Label', 'Go To Max Voxel', 'Callback', {@setMaxVoxel, gH});
end

colormap(cmap);

handles_data.currentFigure = gH;
% if (~set_to_max_voxel)
%     handles_data.centre = voxelOrigin;
% else
%     handles_data.centre = maxVoxelPos;
% end
%else
%    handles_data.numRowsToPlot = size(coords, 1);
%end
set(gH, 'userdata', handles_data);

drawSlices(handles_data);

if (set_to_max_voxel)
    setMaxVoxel(maxVoxelH, [], gH);
end


if (exist('tc', 'var'))
    axesH = subplot( handles_data.numRowsToPlot, 3, (handles_data.numRowsToPlot*3-2:handles_data.numRowsToPlot*3));
    icatb_plot_spectra(axesH, tc);
end

axisH = axes('Parent', gH, 'position', [0 0 1 1], 'visible', 'off');
xPos = 0.5; yPos = 0.97;
titleColor = 'c';
text(xPos, yPos, fig_title, 'color', titleColor, 'fontweight', 'bold', 'fontsize', UI_FS + 1, 'HorizontalAlignment', 'center', 'FontName', UI_FONTNAME, 'parent', axisH);

set(gH, 'visible', 'on');
drawnow;

if (nargout == 1)
    varargout{1} = gH;
end


function drawSlices(handles_data)
%% Draw Slices
%

icatb_defaults;
global FONT_COLOR;

CLIM = handles_data.CLIM;
minInterval = handles_data.minInterval;
maxInterval = handles_data.maxInterval;
coords = handles_data.coords;
structVol = handles_data.structVol;
labels = handles_data.labels;
minICAIM = handles_data.minICAIM;
maxICAIM = handles_data.maxICAIM;
gH = handles_data.currentFigure;

allAxes = zeros(numel(coords), 1);

countA = 0;
for ncoords = 1:size(coords, 1)
    
    realCoords = coords(ncoords, :);
    
    %realCoords = (icatb_voxel_to_real(structVol, coords(ncoords, :)));
    
    [sliceXY, sliceXZ, sliceYZ] = returnSlices(handles_data, realCoords);
    
    countA = countA + 1;
    sh = subplot( handles_data.numRowsToPlot, 3, countA);
    set(sh, 'units', 'normalized');
    plotImage(sh, sliceYZ, CLIM);
    xlabel(['X = ', num2str(realCoords(1)), ' mm'], 'parent', sh);
    firstPos = get(sh, 'position');
    allAxes(countA) = sh;
    
    countA = countA + 1;
    sh = subplot( handles_data.numRowsToPlot, 3, countA);
    set(sh, 'units', 'normalized');
    plotImage(sh, sliceXZ, CLIM);
    xlabel(['Y = ', num2str(realCoords(2)), ' mm'], 'parent', sh);
    
    middleAxes = sh;
    allAxes(countA) = sh;
    
    countA = countA + 1;
    sh = subplot( handles_data.numRowsToPlot, 3, countA);
    set(sh, 'units', 'normalized');
    plotImage(sh, sliceXY, CLIM);
    xlabel(['Z = ', num2str(realCoords(3)), ' mm'], 'parent', sh);
    lastPos = get(sh, 'position');
    allAxes(countA) = sh;
    
    % Plot label
    title(labels{ncoords}, 'parent', middleAxes);
    width = lastPos(1) + lastPos(3) - firstPos(1);
    %         xPos = -width/2 - 0.05/2;
    %         yPos = 0.05/2;
    % text(xPos, yPos, labels{ncoordss}, 'parent', sh, 'units', 'normalized');
    
    % Plot horizontal colorbar
    if (ncoords == size(coords, 1))
        colorBarWidth = 0.3;
        colorBarHeight = 0.016;
        pos = [firstPos(1) + width/2 - colorBarWidth/2, lastPos(2) - 0.08, colorBarWidth, colorBarHeight];
        ch = colorbar('horiz');
        set(ch, 'units', 'normalized');
        set(ch, 'position', pos);
        %ch = colorbar(pos);
        set(ch, 'xLim', [minInterval, maxInterval]);
        xTicks = get(ch, 'xTick');
        set(ch, 'yTick', []);
        set(ch, 'xTick', [xTicks(1), xTicks(end)]);
        set(ch, 'xTicklabel', cellstr(num2str([minICAIM;maxICAIM])));
        set(ch, 'XColor', FONT_COLOR, 'YColor', FONT_COLOR);
        set(ch, 'tag', 'colorbar');
        set(ch, 'units', 'normalized');
        handles_data.ch = ch;
    end
    
end

set(allAxes, 'XColor', get(gH, 'color'), 'YColor', get(gH, 'color'));

xlabelH = get(allAxes, 'XLabel');
% if (iscell(xlabelH))
%     xlabelH = cell2mat(xlabelH);
% end

childAxesH = get(findobj(gH, 'Type', 'Axes'), 'children');
% if (iscell(childAxesH))
%     childAxesH = cell2mat(childAxesH);
% end

for nX = 1:length(xlabelH)
    if (iscell(xlabelH))
        set(xlabelH{nX}, 'color', FONT_COLOR);
    else
        set(xlabelH(nX), 'color', FONT_COLOR);
    end
end

if (handles_data.isInteractive)
    for nX = 1:length(childAxesH)
        set(childAxesH{nX}, 'hittest', 'off');
    end
    %set(childAxesH, 'hittest', 'off');
    set(allAxes(1), 'ButtonDownFcn', {@quad1Callback, gH});
    set(allAxes(2), 'ButtonDownFcn', {@quad2Callback, gH});
    set(allAxes(3), 'ButtonDownFcn', {@quad3Callback, gH});
end

handles_data.allAxes = allAxes;
set(gH, 'userdata', handles_data);

set(gH, 'visible', 'on');

function subH = plotImage(subH, data, CLIM)
%% Function to plot the image at the specified position
%

imageAxisPos = get(subH, 'position');
image(rot90(data), 'parent', subH, 'CDataMapping', 'scaled');
set(subH, 'units', 'normalized');
oldWidth = imageAxisPos(3);
oldHeight = imageAxisPos(4);
yAxisRatio = imageAxisPos(4)/imageAxisPos(3);
xAxisRatio = imageAxisPos(3)/imageAxisPos(4);
if(yAxisRatio>1)
    yAxisRatio = 1;
else
    xAxisRatio = 1;
end
newWidth = imageAxisPos(3)*yAxisRatio;
newHeight = imageAxisPos(4)*xAxisRatio;

if (newHeight < oldHeight)
    imageAxisPos(2) = imageAxisPos(2) + 0.5*(oldHeight - newHeight);
end

imageAxisPos = [imageAxisPos(1), imageAxisPos(2), newWidth, newHeight];
set(subH, 'position', imageAxisPos);
%setImagePos(subH, imagePos);
set(subH, 'clim', CLIM); % set the axis positions to the specified
%axis(subH, 'square');
set(subH, 'XTick', []);
set(subH, 'XTickLabel', []);
set(subH, 'YTick', []);
set(subH, 'YTickLabel', []);


function [sliceXY, sliceXZ, sliceYZ] = returnSlices(handles_data, voxelcoords)
%% Slices (XY, XZ, YZ)


%% XY
[sliceXY, structData] = getData(handles_data.structVol, handles_data.dispPrefs.imgVol, 'axial', voxelcoords(end));
sliceXY = applyDispParams((sliceXY - handles_data.dispPrefs.meanData) / handles_data.dispPrefs.stdData, ...
    handles_data.dispPrefs.returnValue, handles_data.dispPrefs.threshold);
sliceXY = getCompositeData(structData, sliceXY, handles_data.dispPrefs.returnValue, handles_data.minInterval, handles_data.maxInterval, handles_data.range);

%% XZ
[sliceXZ, structData] = getData(handles_data.structVol, handles_data.dispPrefs.imgVol, 'coronal', voxelcoords(2));
sliceXZ = applyDispParams((sliceXZ - handles_data.dispPrefs.meanData) / handles_data.dispPrefs.stdData, ...
    handles_data.dispPrefs.returnValue, handles_data.dispPrefs.threshold);
sliceXZ = getCompositeData(structData, sliceXZ, handles_data.dispPrefs.returnValue, handles_data.minInterval, handles_data.maxInterval, handles_data.range);

%% YZ
[sliceYZ, structData] = getData(handles_data.structVol, handles_data.dispPrefs.imgVol, 'sagittal', voxelcoords(1));
sliceYZ = applyDispParams((sliceYZ - handles_data.dispPrefs.meanData) / handles_data.dispPrefs.stdData, ...
    handles_data.dispPrefs.returnValue, handles_data.dispPrefs.threshold);
sliceYZ = getCompositeData(structData, sliceYZ, handles_data.dispPrefs.returnValue, handles_data.minInterval, handles_data.maxInterval, handles_data.range);

% sliceXY = reshape(data(:, :, voxelcoords(end)), size(data, 1), size(data, 2));
% sliceXZ = reshape(data(:, voxelcoords(2), :), size(data, 1),size(data, 3));
% sliceYZ = reshape(data(voxelcoords(1), :, :), size(data, 2), size(data, 3));


function scaledData = scaleIm(tmin, tmax, imageData, returnValue, data_range)
%% Scale images
%

if (~exist('returnValue', 'var'))
    returnValue = 2;
end

if (~exist('data_range', 'var'))
    minVal = min(imageData);
    maxVal = max(imageData);
else
    minVal = min(data_range);
    maxVal = max(data_range);
end

if (returnValue == 1)
    maxVal = max(abs([minVal, maxVal]));
    minVal = -maxVal;
end

rangeVal = (maxVal-minVal) + eps;
trange = tmax-tmin;
if (rangeVal == 0)
    rangeVal = eps;
end
scaledData = (((imageData-minVal)./rangeVal)./(1/trange))+tmin;

function updateSlices(handles_data)
%% Update slices and labels

%realCoords = (icatb_voxel_to_real(handles_data.structVol, handles_data.coords(1, :)));
realCoords = handles_data.coords(1, :);
[sliceXY, sliceXZ, sliceYZ] = returnSlices(handles_data, realCoords);

%% X
set(findobj(handles_data.allAxes(1), 'Type', 'Image'), 'CData', rot90(sliceYZ));
set(get(handles_data.allAxes(1), 'XLabel'), 'string', ['X = ', num2str(realCoords(1)), ' mm']);

%% Y
set(findobj(handles_data.allAxes(2), 'Type', 'Image'), 'CData', rot90(sliceXZ));
set(get(handles_data.allAxes(2), 'XLabel'), 'string', ['Y = ', num2str(realCoords(2)), ' mm']);

%% Z
set(findobj(handles_data.allAxes(3), 'Type', 'Image'), 'CData', rot90(sliceXY));
set(get(handles_data.allAxes(3), 'XLabel'), 'string', ['Z = ', num2str(realCoords(3)), ' mm']);

set(handles_data.allAxes, 'CLIM', handles_data.CLIM);

set(handles_data.ch, 'xLim', [handles_data.minInterval, handles_data.maxInterval]);
xTicks = get(handles_data.ch, 'xTick');
set(handles_data.ch, 'xTick', [xTicks(1), xTicks(end)]);
set(handles_data.ch, 'xTicklabel', cellstr(num2str([handles_data.minICAIM;handles_data.maxICAIM])));


function setVoxelPos(hObject, event_data, handles)
%% Set Voxel Position

handles_data = get(handles, 'userdata');

% open input dialog box
prompt = {'Enter voxel position in real world coordinates (mm)'};
dlg_title = 'Set voxel position';
num_lines = 1;
def = {num2str(handles_data.coords(1, :))};
%def = {num2str((icatb_voxel_to_real(handles_data.structVol, handles_data.coords)))};

% get voxel coords
voxelPos = icatb_inputdlg2(prompt, dlg_title, num_lines, def);
drawnow;

if (~isempty(voxelPos))
    voxelPos = str2num(voxelPos{1});
    if length(voxelPos) ~= 3
        error('Please enter the three coords correctly');
    end
    %handles_data.coords = round(icatb_real_to_voxel(handles_data.structVol, voxelPos));
    handles_data.coords = voxelPos;
    %voxelOrigin = icatb_real_to_voxel(handles_data.structVol(1), handles_data.origin);
    voxelOrigin = handles_data.voxelOrigin;
    handles_data.pixelPos = round(voxelOrigin + (voxelPos./handles_data.VOX));
    set(handles, 'userdata', handles_data);
    updateSlices(handles_data);
end

function setMaxVoxel(hObject, event_data, handles)
%% Set Max voxel

handles_data = get(handles, 'userdata');
handles_data.coords = handles_data.maxVoxelPos;
handles_data.pixelPos = round(handles_data.voxelOrigin + (handles_data.coords./handles_data.VOX));
%handles_data.centre = handles_data.maxVoxelPos;
set(handles, 'userdata', handles_data);
updateSlices(handles_data);

function [data, structData] = getData(structVol, imgVol, anatomicalView, slices)
%% Get data
%

parameters = icatb_get_slice_def(structVol, 'axial');
transform = parameters.transform;
slicedef = parameters.slicedef;
%
X = 1;Y = 2; Z = 3;
dims = slicedef;
xmm = dims(X,1):dims(X,2):dims(X,3);
ymm = dims(Y,1):dims(Y,2):dims(Y,3);
zmm = parameters.slices;
if (strcmpi(anatomicalView, 'axial'))
    zmm = slices;
end
if (strcmpi(anatomicalView, 'coronal'))
    ymm = slices;
end
if (strcmpi(anatomicalView, 'sagittal'))
    xmm = slices;
end

% if (strcmpi(anatomicalView, 'saggital'))
%     xmm = xmm(indices);
%     tmp = xmm;
%     xmm = ymm;
%     ymm = zmm;
%     zmm = tmp;
%     clear tmp;
% end
%
% if (strcmpi(anatomicalView, 'coronal'))
%     ymm = ymm(indices);
%     tmp = zmm;
%     zmm = ymm;
%     ymm = tmp;
%     clear tmp;
% end
%
% if (strcmpi(anatomicalView, 'axial'))
%     zmm = zmm(indices);
% end

[y, x] = meshgrid(ymm, xmm');
vdims = [length(xmm), length(ymm), length(zmm)];
nslices = vdims(Z);
nvox = prod(vdims(1:2));

structData = zeros(vdims);
data = zeros(vdims);
for ii = 1:nslices
    ixyzmm = [x(:)'; y(:)'; ones(1,nvox)*zmm(ii); ones(1, nvox)];
    vixyz = (transform*imgVol(1).mat) \ ixyzmm;
    i1 = icatb_spm_sample_vol(imgVol(1), vixyz(X, :), vixyz(Y, :), vixyz(Z, :), [1, nan]);
    data(:, :, ii) = reshape(i1, vdims(1:2));
    vixyz = (transform*structVol(1).mat) \ ixyzmm;
    i1 = icatb_spm_sample_vol(structVol(1), vixyz(X, :), vixyz(Y, :), vixyz(Z, :), [1, nan]);
    structData(:, :, ii) = reshape(i1, vdims(1:2));
end

data(isfinite(data) == 0) = 0;
structData(isfinite(structData) == 0) = 0;

data = squeeze(data);
structData = squeeze(structData);

function data = applyDispParams(data, returnValue, threshold)
%% Apply display parameters

threshold = abs(threshold);

if (returnValue == 2)
    % Positive
    data(data < 0) = 0;
elseif (returnValue == 3)
    % Absolute value
    data = abs(data);
elseif (returnValue == 4)
    % Negative value
    data(data > 0) = 0;
end

if (length(threshold) > 1)
    
    if ((returnValue == 2) || (returnValue == 3))
        % Positive or absolute
        data(data < min(threshold)) = 0;
        data(data > max(threshold)) = max(threshold);
    elseif (returnValue == 4)
        % Negative
        data(abs(data) < min(threshold)) = 0;
        data(abs(data) > max(threshold)) = -max(threshold);
    else
        % Positive and Negative
        neg_threshold = -threshold;
        
        pos_inds = (data > 0);
        neg_inds = (data < 0);
        
        % Handle positive range
        tmp1 = data(pos_inds);
        tmp1(tmp1 < min(threshold)) = 0;
        tmp1(tmp1 > max(threshold)) = max(threshold);
        
        % handle negative range
        tmp2 = data(neg_inds);
        tmp2(tmp2 > max(neg_threshold)) = 0;
        tmp2(tmp2 < min(neg_threshold)) = min(neg_threshold);
        
        data(pos_inds) = tmp1;
        data(neg_inds) = tmp2;
    end
    
    %data(abs(data) < min(threshold)) = 0;
    %data(abs(data) > max(threshold)) = 0;
else
    data(abs(data) < threshold) = 0;
end

function composite_data = getCompositeData(structData, data, returnValue, minInterval, maxInterval, data_range)
%% Get composite data

good_inds = find(abs(data) > eps);

unitColor = icatb_range(data(:))/64;

if (exist('data_range', 'var'))
    data(data == data_range(2)) = data(data == data_range(2)) - unitColor;
end

% Overlay components
composite_data = scaleIm(maxInterval + 1, 2*maxInterval, structData(:));
if (~isempty(good_inds))
    if (~exist('data_range', 'var'))
        composite_data(good_inds) = scaleIm(minInterval, maxInterval, data(good_inds), returnValue);
    else
        composite_data(good_inds) = scaleIm(minInterval, maxInterval, data(good_inds), returnValue, data_range);
    end
end
clear data;

composite_data = reshape(composite_data, size(structData));

function quad1Callback(hObject, event_data, handles)
%% Quad1

handles_data = get(handles, 'userdata');
pixelPos = handles_data.pixelPos;
childH = get(hObject, 'children');
data = get(childH, 'CData')';
imDIM = size(data);
points = get(hObject, 'currentPoint');
x = points(1, 1);
y = points(1, 2);
newY = ceil(x);
newZ = ceil(imDIM(end) - y);
if ((newY < 1) || (newY > imDIM(1)))
    return;
end
if ((newZ < 1) || (newZ > imDIM(end)))
    return;
end
pixelPos(2:3) = [newY, newZ];
handles_data.coords = (pixelPos - handles_data.voxelOrigin).*handles_data.VOX;
handles_data.pixelPos = pixelPos;
set(handles, 'userdata', handles_data);
updateSlices(handles_data);

function quad2Callback(hObject, event_data, handles)
%% Quad2

handles_data = get(handles, 'userdata');
pixelPos = handles_data.pixelPos;
childH = get(hObject, 'children');
data = get(childH, 'CData')';
imDIM = size(data);
points = get(hObject, 'currentPoint');
x = points(1, 1);
y = points(1, 2);
newX = ceil(x);
newZ = ceil(imDIM(end)-y);

if ((newX < 1) || (newX > imDIM(1)))
    return;
end

if ((newZ < 1) || (newZ > imDIM(end)))
    return;
end

pixelPos([1, 3]) = [newX, newZ];
handles_data.coords = (pixelPos - handles_data.voxelOrigin).*handles_data.VOX;
handles_data.pixelPos = pixelPos;
set(handles, 'userdata', handles_data);
updateSlices(handles_data);

function quad3Callback(hObject, event_data, handles)
%% Quad3

handles_data = get(handles, 'userdata');
pixelPos = handles_data.pixelPos;
childH = get(hObject, 'children');
data = get(childH, 'CData')';
imDIM = size(data);
points = get(hObject, 'currentPoint');
x = points(1, 1);
y = points(1, 2);
newX = ceil(x);
newY = ceil(imDIM(2)-y);

if ((newX < 1) || (newX > imDIM(1)))
    return;
end

if ((newY < 1) || (newY > imDIM(end)))
    return;
end

pixelPos(1:2) = [newX, newY];
handles_data.coords = (pixelPos - handles_data.voxelOrigin).*handles_data.VOX;
handles_data.pixelPos = pixelPos;
set(handles, 'userdata', handles_data);
updateSlices(handles_data);