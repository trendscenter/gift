function  varargout = icatb_display_composite(files, varargin)
%% Make composite plots. If no output parameter is specified, figure is displayed with the orthogonal sections.
%
% Inputs:
% 1. files - Cell array of strings
% 2. varargin - Arguments passed in pairs:
%   a. convert_to_zscores - Options are 'yes' or 'no'
%   b. image_values - Options are 'positive', 'positive and negative', 'absolute
%   value', 'negative'
%   c. anatomical_file - Anatomical file full path.
%   d. threshold - Threshold value
%   e. colorbar_label - Labels for each component in a cell array
%   d. title - Title of the ortho-slices
%

icatb_defaults;
global FONT_COLOR;
global TEXT_DISPLAY_SLICES_IN_MM;
global UI_FONTNAME;
global UI_FS;
global USE_UNIFORM_COLOR_COMPOSITE;


%% Parse params
threshold = 1.5;
convert_to_zscores = 'yes';
image_values = 'positive';
template_file = fullfile(fileparts(which('groupica.m')), 'icatb_templates', 'ch2bet.nii');
slice_plane = 'axial';
slices_in_mm = [];
axesTitle = '';
display_type = 'ortho';
uniformColor = USE_UNIFORM_COLOR_COMPOSITE;
set_to_max_voxel = 0;

for nF = 1:2:length(varargin)
    if (strcmpi(varargin{nF}, 'convert_to_zscores') || strcmpi(varargin{nF}, 'convert_to_z'))
        convert_to_zscores = varargin{nF + 1};
    elseif (strcmpi(varargin{nF}, 'image_values'))
        image_values = varargin{nF + 1};
    elseif (strcmpi(varargin{nF}, 'anatomical_file'))
        template_file = varargin{nF + 1};
    elseif (strcmpi(varargin{nF}, 'anatomical_view'))
        slice_plane = varargin{nF + 1};
    elseif (strcmpi(varargin{nF}, 'threshold'))
        threshold = varargin{nF + 1};
    elseif (strcmpi(varargin{nF}, 'colorbar_label'))
        labels = varargin{nF + 1};
    elseif (strcmpi(varargin{nF}, 'title'))
        axesTitle = varargin{nF + 1};
    elseif (strcmpi(varargin{nF}, 'display_type'))
        display_type = varargin{nF + 1};
    elseif (strcmpi(varargin{nF}, 'slices_in_mm'))
        slices_in_mm = varargin{nF + 1};
    elseif (strcmpi(varargin{nF}, 'coords'))
        coordsIn = varargin{nF + 1};
    elseif (strcmpi(varargin{nF}, 'cmap'))
        cmap = varargin{nF + 1};
    end
end


if (strcmpi(display_type, 'render'))
    global prevrend;
    rendfile = fullfile(fileparts(which('groupica.m')), 'icatb_templates', 'render_single_subj.mat');
    prevrend.rendfile = rendfile;
    prevrend.col = [1 0 0; 0 1 0; 0 0 1];
    prevrend.brt = NaN;
end

if (isempty(uniformColor))
    uniformColor = 0;
end

if (~exist('files', 'var') || isempty(files))
    files = icatb_selectEntry('title', 'Select images ...', 'typeSelection', 'multiple', 'typeEntity', 'file', 'filter', '*img;*nii', 'fileType', 'image');
    drawnow;
end

%% Parse params
files = icatb_rename_4d_file(files);
files = cellstr(files);

if (~exist('labels', 'var'))
    labels = cellstr(num2str((1:length(files))'));
end

if (~strcmpi(display_type, 'montage'))
    slices_in_mm = [];
end

structVol = icatb_spm_vol(template_file);

returnValue = strmatch(lower(image_values), {'positive and negative', 'positive', 'absolute value', 'negative'}, 'exact');

colorLen = 64;

% get colormap
if (exist('cmap', 'var'))
    
    if (iscell(cmap))
        cmap = cat(1, cmap{:});
    end
    
    if (length(files)*colorLen ~= size(cmap, 1))
        error('Specify colormap of length 64 for each file');
    end
    
    imagesCmap = [cmap; gray(colorLen)];
    
else
    imagesCmap = getCmap(length(files), colorLen, uniformColor);
end

if (~strcmpi(display_type, 'render'))
    
    topInds = [];
    
    maxCoords = zeros(length(files), 3);
    
    minMaxLabels = repmat({''}, length(files), 2);
    
    for nComp = 1:length(files)
        
        fileName = files{nComp};
        
        [tmpF, tmpN] = icatb_parseExtn(fileName);
        
        [imageTmp, coordsAll, HInfo, slices_in_mm] = icatb_resizeImage(structVol, tmpF, slice_plane, slices_in_mm, tmpN);
        
        
        structDIM = [size(imageTmp, 2), size(imageTmp, 3), size(imageTmp, 4)];
        structData = reshape(imageTmp(1, :, :, :), structDIM);
        imageTmp = reshape(imageTmp(end, :, :, :), structDIM);
        
        permuteDIM = [1, 2, 3];
        
        if strcmpi(slice_plane, 'sagittal')
            permuteDIM = [2, 3, 1];
        end
        
        if strcmpi(slice_plane, 'coronal')
            permuteDIM = [1, 3, 2];
        end
        
        structDIM = structDIM(permuteDIM);
        structData = permute(structData, permuteDIM);
        imageTmp = permute(imageTmp, permuteDIM);
        HInfo.VOX = HInfo.VOX(permuteDIM);
        HInfo.DIM = structDIM;
        coordsAll = permute(coordsAll, [permuteDIM, 4]);
        
        if (nComp == 1)
            compositeMap = scaleIm(colorLen*length(files) + 1, colorLen*(length(files) + 1 ), structData(:));
        end
        
        imageTmp = icatb_applyDispParameters_comp(imageTmp(:)', strcmpi(convert_to_zscores, 'yes'), returnValue, threshold);
        
        mask_inds = (abs(imageTmp) > eps);
        
        mask_inds = find(mask_inds == 1);
        
        
        
        %masks(:, nComp) = mask_inds;
        
        if (~isempty(mask_inds))
            if (returnValue == 1)
                minMaxTmp = [-max(abs(imageTmp(mask_inds))), max(abs(imageTmp(mask_inds)))];
            else
                minMaxTmp = [min(imageTmp(mask_inds)), max(imageTmp(mask_inds))];
            end
            minMaxLabels{nComp, 1} = num2str(minMaxTmp(1), '%0.1f');
            minMaxLabels{nComp, 2} = num2str(minMaxTmp(2), '%0.1f');
            compositeMap(mask_inds) = scaleIm(colorLen*(nComp - 1) + 1, colorLen*nComp, imageTmp(mask_inds), returnValue);
            
            [dd, tmpTopInds] = sort(abs(imageTmp(mask_inds)));
            tmpTopInds = tmpTopInds(end:-1:1);
            tmpTopInds = mask_inds(tmpTopInds);
            %ddMin = min([length(tmpTopInds), 50]);
            %tmpTopInds = tmpTopInds(1:ddMin);
            tmpTopInds = tmpTopInds(:)';
            topInds = [topInds, tmpTopInds];
            [xx,yy,zz]=ind2sub(structDIM, tmpTopInds(1));
            maxCoords(nComp, :) = squeeze(coordsAll(xx, yy, zz, :));
            
        end
        
    end
    
    handles_data.colorLen = colorLen;
    handles_data.real_world_coords = coordsAll;
    handles_data.VOX = HInfo.VOX;
    handles_data.origin = [0, 0, 0];
    
    voxelOrigin = icatb_real_to_voxel(structVol(1), handles_data.origin );
    coords = handles_data.origin;
    pixelPos = voxelOrigin;
    if (~isempty(topInds))
        %% figure out which voxels have most components and return the first in the list
        uniqueInds = unique(topInds);
        NN = histc(topInds, uniqueInds);
        
        [dd, inds] = max(NN);
        [xx,yy,zz]=ind2sub(structDIM, uniqueInds(inds));
        
        coords = squeeze(coordsAll(xx, yy, zz, :))';
        coords = coords(:)';
        pixelPos = [xx, yy, zz];
    end
    
    if (length(files) == 1)
        coords = maxCoords;
        %pixelPos = round(voxelOrigin + (coords./handles_data.VOX));
    end
    
    if (exist('coordsIn', 'var'))
        coords = coordsIn;
    end
    
    pixelPos = round(voxelOrigin + (coords(1, :)./handles_data.VOX));
    
    compositeMap = reshape(ceil(compositeMap), structDIM);
    maxVoxelPos = coords;
    
    % Store handles data
    handles_data.coords = coords;
    handles_data.voxelOrigin = voxelOrigin;
    handles_data.structVol = structVol;
    handles_data.pixelPos = pixelPos;
    handles_data.maxVoxelPos = maxVoxelPos;
    handles_data.data = compositeMap;
    
    compMenuLabels = cell(1, length(files));
    for nComp = 1:length(files)
        compMenuLabels{nComp} = ['Max ', labels{nComp}];
    end
    
    handles_data.compMenuLabels = compMenuLabels;
    handles_data.cmap = imagesCmap;
    handles_data.title = axesTitle;
    handles_data.maxCoords = maxCoords;
    handles_data.minMaxLabels = minMaxLabels;
    
    
    if (strcmpi(display_type, 'montage'))
        [im, numImagesX, numImagesY, slices_in_mm] = icatb_returnMontage(compositeMap, [], size(compositeMap), handles_data.VOX, slices_in_mm);
        im = reshape(imagesCmap(im(:), :), [size(im), 3]);
    end
    
    if (nargout == 1)
        if (strcmpi(display_type, 'ortho'))
            
            [stackedData, xlabels] = returnStackedData(handles_data, coords);
            handles_data.slices = stackedData;
            handles_data.xlabels = xlabels;
        elseif (strcmpi(display_type, 'render'))
            handles_data.slices = rgb;
        else
            handles_data.slices = im;
            handles_data.slices_in_mm = slices_in_mm;
        end
        varargout{1} = handles_data;
        return;
    end
    
    gH = icatb_getGraphics('Composite Maps', 'graphics', 'composite_maps', 'on');
    set(gH, 'resize', 'on');
    
    if (exist('fig_pos', 'var'))
        set(gH, 'position', fig_pos);
    end
    
    colormap(imagesCmap);
    
    handles_data.currentFigure = gH;
    
    
    if (strcmpi(display_type, 'ortho'))
        
        cmenu = uicontextmenu;
        set(gH, 'UIContextMenu', cmenu);
        uimenu(cmenu, 'Label', 'Set Voxel Position (mm)', 'Callback', {@setVoxelPos, gH});
        uimenu(cmenu, 'Label', 'Composite', 'Callback', {@setMaxVoxel, gH});
        for nComp = 1:length(files)
            uimenu(cmenu, 'Label', compMenuLabels{nComp}, 'Callback', {@setMaxVoxel, gH});
        end
        
        handles_data.compMenuLabels = compMenuLabels;
        
        set(gH, 'userdata', handles_data);
        
        % Draw slices and plot colorbar
        sh = drawSlices(gH);
        
    else
        
        axesPos = [0.15, 0.25, 0.7, 0.7];
        
        sh = axes('parent', gH, 'units', 'normalized', 'position', axesPos, 'tag', 'image_axes');
        imagesc(im);
        colormap(imagesCmap);
        set(sh, 'CLIM', [1, size(imagesCmap, 1)]);
        axis(sh,'image');
        axis(sh, 'off');
        
        if strcmpi(TEXT_DISPLAY_SLICES_IN_MM, 'on')
            dim(1) = size(im, 1);
            dim(2) = size(im, 2);
            % name the text and place it in the correct order (slices in mm).
            textCount = 0;
            yPos = 1 + dim(1) / numImagesY;
            for nTextRows = 1:numImagesY
                xPos = 1;
                for nTextCols = 1:numImagesX
                    textCount = textCount + 1;
                    if textCount <= structDIM(3)
                        txHandle(textCount) = text(xPos, yPos, num2str(round(slices_in_mm(textCount))), 'color', FONT_COLOR,  ...
                            'fontsize', UI_FS-4, 'HorizontalAlignment', 'left', 'verticalalignment', 'bottom', ...
                            'FontName', UI_FONTNAME, 'parent', sh);
                    end
                    xPos = xPos + (dim(2) / numImagesX);
                end
                % end for cols
                yPos = yPos + (dim(1) / numImagesY); % update the y position
            end
            % end for rows
            
        end
        
        %     nrows = ceil(sqrt(length(slices_in_mm)));
        %     ncols = ceil(length(slices_in_mm)./nrows);
        
        %     count = 0;
        %     for nr = 1:nrows
        %         for nc = 1:ncols
        %             count = count + 1;
        %             if (count > length(slices_in_mm))
        %                 break;
        %             end
        %             sh = subplot(nrows, ncols, count);
        %             dd = rot90(squeeze(compositeMap(:, :, count)));
        %             dd = reshape(imagesCmap(dd(:), :), [size(dd), 3]);
        %             imagesc(dd);
        %             xlabel(num2str(slices_in_mm(count)), 'parent', sh, 'horizontalAlignment', 'left', 'color', FONT_COLOR);
        %             axis image;
        %             %axis off;
        %             set(sh, 'XColor', get(gH, 'color'));
        %             set(sh, 'YColor', get(gH, 'color'));
        %             set(sh, 'XTickLabel', []);
        %             set(sh, 'YTickLabel', []);
        %             set(sh, 'box', 'off');
        %             set(sh, 'color', 'none');
        %             set(sh, 'CLIM', [1, size(imagesCmap, 1)]);
        %         end
        %     end
        
        %colormap(imagesCmap);
        
    end
    
else
    % render
    
    minMaxLabels = repmat({''}, length(files), 2);
    rendMat = cell(1, length(files));
    MinMaxVals = zeros(length(files), 2);
    
    for nV = 1:length(files)
        
        [V, HInfo] = icatb_returnHInfo(files{nV});
        [R, C, P] = ndgrid(1:V(1).dim(1),1:V(1).dim(2),1:V(1).dim(3));
        RCP = [R(:)';C(:)';P(:)'];
        data = icatb_spm_read_vols(V);
        data(isfinite(data) == 0) = 0;
        data = icatb_applyDispParameters(data(:)', strcmpi(convert_to_zscores, 'yes'), returnValue, threshold, V.dim(1:3), HInfo);
        mask = find(abs(data) > eps);
        
        %if (~isempty(mask))
        
        minMaxTmp = [min(data(mask)), max(data(mask))];
        minMaxLabels{nV, 1} = num2str(minMaxTmp(1), '%0.1f');
        minMaxLabels{nV, 2} = num2str(minMaxTmp(2), '%0.1f');
        
        MinMaxVals(nV, :) = minMaxTmp;
        
        XYZ = RCP(:, mask);
        Z = data(mask); %ones(length(mask), 1);
        dat   = struct(	'XYZ',	XYZ,...
            't',	Z, ...
            'mat',	V.mat, ...
            'dim',	V.dim(1:3)');
        
        [blah, rendMat{nV}] = icatb_spm_render(dat, NaN, rendfile);
        
    end
    
    countVoxelsRend = zeros(length(rendMat), length(rendMat{1}));
    rgb = cell(1, length(rendMat{1}));
    for i = 1:length(rendMat{1})
        
        for nV = 1:length(rendMat)
            
            if (any(MinMaxVals(nV, :) < 0))
                CLIM = [-max(abs(MinMaxVals(nV, :))), max(abs(MinMaxVals(nV, :)))];
            else
                CLIM = [min(MinMaxVals(nV, :)), max(MinMaxVals(nV, :))];
            end
            
            XD = rendMat{nV}{i}.data{1};
            
            if (nV == 1)
                ren = rendMat{nV}{i}.ren;
                minCm = (colorLen*length(files)) + 1;
                maxCm = colorLen*(length(files) + 1);
                ren = reshape(scaleIm(minCm, maxCm, ren(:), 2), size(ren));
            end
            
            msk = find(abs(XD) > eps);
            countVoxelsRend(nV, i) = length(msk);
            minCm = (colorLen*(nV - 1)) + 1;
            maxCm = colorLen*nV;
            
            %             if (returnValue == 1)
            %                 minCmA = minCm;
            %                 minCmB = ceil(minCm + (colorLen/2) - 1);
            %                 minCmC = minCmB + 1;
            %                 minCmD = maxCm;
            %             else
            %                 minCmA = minCm;
            %                 minCmB = maxCm;
            %                 minCmC = minCm;
            %                 minCmD = maxCm;
            %             end
            
            if (~isempty(msk))
                ren(msk) = scaleIm(minCm, maxCm, XD(msk), returnValue);
            end
            
            %             msk = find(XD > 0);
            %             if (~isempty(msk))
            %                 ren(msk) = scaleIm(minCmC, minCmD, XD(msk), returnValue, [0, CLIM(2)]);
            %             end
            
        end
        
        if (i > 2)
            ren = round(flipud(ren));
        else
            ren = round(rot90(ren));
        end
        
        ren = reshape(imagesCmap(ren(:), :), [size(ren), 3]);
        
        rgb{i} = ren;
        
        clear ren;
        
    end
    
    slices_rgb = rgb;
    rgb = concatDat(rgb);
    
    if (nargout == 1)
        handles_data.slices = rgb;
        handles_data.slices_all = slices_rgb;
        handles_data.cmap = imagesCmap;
        handles_data.minMaxLabels = minMaxLabels;
        handles_data.colorLen = colorLen;
        handles_data.countVoxelsRend = countVoxelsRend;
        varargout{1} = handles_data;
        return;
    end
    
    axesPos = [0.15, 0.25, 0.7, 0.7];
    gH = icatb_getGraphics('Composite Renderer', 'graphics', 'Image Viewer', 'on');
    sh = axes('parent', gH, 'units', 'normalized', 'position', axesPos);
    colormap(imagesCmap);
    
    ImageAxis = image(rgb, 'parent', sh);
    axis(sh, 'image');
    axis(sh, 'off');
    set(sh, 'CLIM', [1, size(imagesCmap, 1)]);
    
end


axesPos = get(sh, 'position');

offset = 0;
maxWidth = 0.9;
cwdiths = (maxWidth-offset*length(files))/(length(files));

if (cwdiths > 0.2)
    cwdiths = 0.2;
end

pos(1) = 0.05;

pos(3) = cwdiths;
pos(2) = axesPos(2) - 0.16;
pos(4) = 0.035;

pos(1) = 0.5 - 0.5*(length(files)*pos(3));

if (pos(1) < 0.05)
    pos(1) = 0.05;
end

cMax = 0;
allColorbars = zeros(1, length(files));
for n = 1:length(files)
    cMin = cMax + 1;
    cMax = cMax + colorLen;
    ch = colorbar('horiz');
    set(ch, 'tag', ['Colorbar', num2str(n)]);
    set(ch, 'position', pos);
    set(ch, 'xlim', [1, size(imagesCmap, 1)]);
    pos(1) = pos(1) + pos(3) + offset;
    set(ch, 'xlim', [cMin, cMax]);
    set(ch, 'xtick', []);
    set(ch, 'color', FONT_COLOR);
    xlabel(char(minMaxLabels{n, 1}, '', labels{n}), 'parent', ch, 'fontsize', UI_FS-4);
    title(minMaxLabels{n, 2}, 'parent', ch, 'fontsize', UI_FS-4, 'color', FONT_COLOR);
    set(ch, 'color', FONT_COLOR);
    allColorbars(n) = ch;
end

handles_data = get(gH, 'userdata');

handles_data.allColorbars = allColorbars;
set(gH, 'userdata', handles_data);

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



function quad1Callback(hObject, event_data, handles)
%% Quad1

handles_data = get(handles, 'userdata');
pixelPos = handles_data.pixelPos;
childH = get(hObject, 'children');
%ata = get(childH, 'CData')';
imDIM = [size(handles_data.data, 2), size(handles_data.data, 3)];
%imDIM = size(data);
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
[handles_data.coords, voxelCoord] = getCoord(handles_data.structVol, handles_data.real_world_coords, pixelPos);
%handles_data.coords = (pixelPos - handles_data.voxelOrigin).*handles_data.VOX;
handles_data.pixelPos = pixelPos;
set(handles, 'userdata', handles_data);
updateSlices(handles_data);

function quad2Callback(hObject, event_data, handles)
%% Quad2

handles_data = get(handles, 'userdata');
pixelPos = handles_data.pixelPos;
childH = get(hObject, 'children');
%data = get(childH, 'CData');
%imDIM = size(data);
imDIM = [size(handles_data.data, 1), size(handles_data.data, 3)];
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
%handles_data.coords = (pixelPos - handles_data.voxelOrigin).*handles_data.VOX;
[handles_data.coords, voxelCoord] = getCoord(handles_data.structVol, handles_data.real_world_coords, pixelPos);
handles_data.pixelPos = pixelPos;
set(handles, 'userdata', handles_data);
updateSlices(handles_data);

function quad3Callback(hObject, event_data, handles)
%% Quad3

handles_data = get(handles, 'userdata');
pixelPos = handles_data.pixelPos;
childH = get(hObject, 'children');
%data = get(childH, 'CData')';
imDIM = [size(handles_data.data, 1), size(handles_data.data, 2)];
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
[handles_data.coords, voxelCoord] = getCoord(handles_data.structVol, handles_data.real_world_coords, pixelPos);
%handles_data.coords = (pixelPos - handles_data.voxelOrigin).*handles_data.VOX;
handles_data.pixelPos = pixelPos;
set(handles, 'userdata', handles_data);
updateSlices(handles_data);


function sh = drawSlices(gH)
%% Draw Slices
%

icatb_defaults;
global FONT_COLOR;

handles_data = get(gH, 'userdata');
colorLen = handles_data.colorLen;
coords = handles_data.coords;
gH = handles_data.currentFigure;
icaPos = handles_data.pixelPos;
cmap = handles_data.cmap;

allAxes = zeros(3, 1);

countA = 0;
ncoords = 1;
%for ncoords = 1:size(coords, 1)

realCoords = coords(ncoords, :);
%voxelCoord = round(icatb_real_to_voxel(handles_data.structVol, realCoords));
voxelCoord = handles_data.pixelPos;
sliceXY = rot90(squeeze(handles_data.data(:,:,voxelCoord(3))));
sliceXZ = rot90(squeeze(handles_data.data(:,voxelCoord(2), :)));
sliceYZ = rot90(squeeze(handles_data.data(voxelCoord(1),:,:)));

sliceXY = reshape(cmap(sliceXY(:), :), [size(sliceXY), 3]);
sliceXZ = reshape(cmap(sliceXZ(:), :), [size(sliceXZ), 3]);
sliceYZ = reshape(cmap(sliceYZ(:), :), [size(sliceYZ), 3]);


%realCoords = (icatb_voxel_to_real(structVol, coords(ncoords, :)));

%[sliceXY, sliceXZ, sliceYZ] = returnSlices(handles_data, realCoords);

countA = countA + 1;
sh = subplot(1, 3, countA);
set(sh, 'units', 'normalized');
plotImage(sh, sliceYZ);
set(sh, 'CLIM', [1, size(cmap, 1)]);
%imagesc(sliceYZ);
%plotImage(sh, sliceYZ, CLIM);
xlabel(['X = ', num2str(realCoords(1)), ' mm'], 'parent', sh, 'color', FONT_COLOR);
firstPos = get(sh, 'position');
allAxes(countA) = sh;

countA = countA + 1;
sh = subplot(1, 3, countA);
set(sh, 'units', 'normalized');
plotImage(sh, sliceXZ);
set(sh, 'CLIM', [1, size(cmap, 1)]);
%plotImage(sh, sliceXZ, CLIM);
xlabel(['Y = ', num2str(realCoords(2)), ' mm'], 'parent', sh, 'color', FONT_COLOR);

middleAxes = sh;
allAxes(countA) = sh;

countA = countA + 1;
sh = subplot( 1, 3, countA);
set(sh, 'units', 'normalized');
plotImage(sh, sliceXY);
set(sh, 'CLIM', [1, size(cmap, 1)]);
%plotImage(sh, sliceXY, CLIM);
xlabel(['Z = ', num2str(realCoords(3)), ' mm'], 'parent', sh, 'color', FONT_COLOR);
lastPos = get(sh, 'position');
allAxes(countA) = sh;

% Plot label
%  title(labels{ncoords}, 'parent', middleAxes);
width = lastPos(1) + lastPos(3) - firstPos(1);
%         xPos = -width/2 - 0.05/2;
%         yPos = 0.05/2;
% text(xPos, yPos, labels{ncoordss}, 'parent', sh, 'units', 'normalized');

%end

%set(allAxes, 'XColor', get(gH, 'color'), 'YColor', get(gH, 'color'));

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

%if (handles_data.isInteractive)
for nX = 1:length(childAxesH)
    set(childAxesH{nX}, 'hittest', 'off');
end
%set(childAxesH, 'hittest', 'off');
set(allAxes(1), 'ButtonDownFcn', {@quad1Callback, gH});
set(allAxes(2), 'ButtonDownFcn', {@quad2Callback, gH});
set(allAxes(3), 'ButtonDownFcn', {@quad3Callback, gH});
%end

handles_data.allAxes = allAxes;
set(gH, 'userdata', handles_data);
title(handles_data.title, 'parent', allAxes(2), 'fontsize', 12, 'color', FONT_COLOR);

set(gH, 'visible', 'on');



% VIVID Creates a Personalized Vivid Colormap
%
%  VIVID(M,...) Creates a colormap with M colors
%  VIVID(MINMAX,...) Creates a colormap with a custom intensity range
%  VIVID(CLRS,...) Creates a colormap with custom colors
%  CMAP = VIVID(...) Exports the vivid colormap to a Mx3 matrix
%
%   Inputs:
%       M - (optional) an integer between 1 and 256 specifying the number
%           of colors in the colormap. Default is 128.
%       MINMAX - (optional) is a 1x2 vector with values between 0 and 1
%           representing the intensity range for the colors, which correspond
%           to black and white, respectively. Default is [0.15 0.85].
%       CLRS - (optional) either a Nx3 matrix of values between 0 and 1
%           representing the desired colors in the colormap
%               -or-
%           a string of characters that includes any combination of the
%           following letters:
%               'r' = red           'g' = green         'b' = blue
%               'y' = yellow        'c' = cyan          'm' = magenta
%               'o' = orange        'l' = lime green    'a' = aquamarine
%               's' = sky blue      'v' = violet        'p' = pink
%               'n' = navy blue     'f' = forest green
%               'k' or 'w' = black/white/grayscale
%
%   Outputs:
%       CMAP - an Mx3 colormap matrix
%
%   Example:
%       % Default Colormap
%       imagesc(sort(rand(200),'descend'));
%       colormap(vivid); colorbar
%
%   Example:
%       % Mapping With 256 Colors
%       imagesc(peaks(500))
%       colormap(vivid(256)); colorbar
%
%   Example:
%       % Mapping With Full Intensity Range
%       mesh(peaks(500))
%       colormap(vivid([0 1])); colorbar
%
%   Example:
%       % Mapping With Light Colors
%       mesh(peaks(500))
%       colormap(vivid([.5 1])); colorbar
%
%   Example:
%       % Mapping With Dark Colors
%       mesh(peaks(500))
%       colormap(vivid([0 .5])); colorbar
%
%   Example:
%       % Mapping With Custom Color Matrix
%       imagesc(peaks(500))
%       clrs = [.5 0 1; 0 .5 1; 0 1 .5; .5 1 0; 1 .5 0; 1 0 .5;];
%       colormap(vivid(clrs)); colorbar
%
%   Example:
%       % Mapping With Color String
%       imagesc(peaks(500))
%       colormap(vivid('pmvbscaglyor')); colorbar
%
%   Example:
%       % Colormap With Multiple Custom Settings
%       imagesc(sort(rand(300,100),'descend'));
%       colormap(vivid(64,[.1 .9],'bwr')); colorbar
%
%   Example:
%       % Topo Colormap
%       load topo;
%       imagesc(topo); axis xy; caxis([-6000 6000])
%       colormap(vivid('bg')); colorbar
%
% See also: jet, hsv, gray, hot, cold, copper, bone, fireice
%
% Author: Joseph Kirk
% Email: jdkirk630@gmail.com
% Release: 1.1
% Date: 12/09/11
function cmap = vivid(varargin)


% Default Color Spectrum
clrs = [1 0 1;.5 0 1;0 0 1;0 1 1    % Magenta, Violet, Blue, Cyan
    0 1 0;1 1 0;1 .5 0;1 0 0];      % Green, Yellow, Orange, Red

% Default Min/Max Intensity Range
minmax = [0.15 0.85];

% Default Colormap Size
m = 256;

% Process Inputs
for var = varargin
    input = var{1};
    if ischar(input)
        nColors = length(input(:));
        colorMat = zeros(nColors,3);
        c = 0;
        for k = 1:nColors
            c = c + 1;
            switch lower(input(k))
                case 'r', colorMat(c,:) = [1 0 0];  % red
                case 'g', colorMat(c,:) = [0 1 0];  % green
                case 'b', colorMat(c,:) = [0 0 1];  % blue
                case 'y', colorMat(c,:) = [1 1 0];  % yellow
                case 'c', colorMat(c,:) = [0 1 1];  % cyan
                case 'm', colorMat(c,:) = [1 0 1];  % magenta
                case 'p', colorMat(c,:) = [1 0 .5]; % pink
                case 'o', colorMat(c,:) = [1 .5 0]; % orange
                case 'l', colorMat(c,:) = [.5 1 0]; % lime green
                case 'a', colorMat(c,:) = [0 1 .5]; % aquamarine
                case 's', colorMat(c,:) = [0 .5 1]; % sky blue
                case 'v', colorMat(c,:) = [.5 0 1]; % violet
                case 'f', colorMat(c,:) = [0 .5 0]; % forest green
                case 'n', colorMat(c,:) = [0 0 .5]; % navy
                case {'k','w'}, colorMat(c,:) = [.5 .5 .5]; % grayscale
                otherwise, c = c - 1;
                    fprintf('Warning: Input character [%s] is not a recognized color ...\n',input(k));
            end
        end
        colorMat = colorMat(1:c,:);
        if ~isempty(colorMat)
            clrs = colorMat;
        end
    elseif isnumeric(input)
        if isscalar(input)
            m = max(1,min(256,round(real(input))));
        elseif size(input,2) == 3
            clrs = max(0,min(1,real(input)));
        elseif length(input) == 2
            minmax = max(0,min(1,real(input)));
        end
    end
end

% Calculate Parameters
nc = size(clrs,1);  % number of spectrum colors
ns = ceil(m/nc);    % number of shades per color
n = nc*ns;
d = n - m;

% Scale Intensity
sup = 2*minmax;
sub = 2*minmax - 1;
if ns == 1
    high = repmat(min(1,mean(sup)),[1 nc 3]);
    low = repmat(max(0,mean(sub)),[1 nc 3]);
else
    high = repmat(min(1,linspace(sup(1),sup(2),ns))',[1 nc 3]);
    low = repmat(max(0,linspace(sub(1),sub(2),ns))',[1 nc 3]);
end

% Determine Color Spectrum
rgb = repmat(reshape(clrs,1,nc,3),ns,1);
map = rgb.*high + (1-rgb).*low;

% Obtain Color Map
cmap = reshape(map,n,3,1);
cmap(1:ns:d*ns,:) = [];



function cmap = getCmap(numComp, colorLen, uniformColor)

if (~exist('uniformColor', 'var'))
    uniformColor = 0;
end

cmaps = cell(1, numComp);

if (~uniformColor)
    load icatb_colors;
    colors = {redV, blueV, greenV, pinkV, yellowV, 'c', 'o', 'p', 'v', 'y', 'm', 'g', 'r', 'l', 'a', 'b', 's'};
    
    for nC = 1:length(cmaps)
        if (nC <= 5)
            tmp = colors{nC};
            colorsToSkip = ceil(size(tmp, 1)/colorLen);
            cmaps{nC} = tmp(1:colorsToSkip:end, :);
        elseif ((nC > 5) && (nC <= length(colors)))
            cmaps{nC} = vivid(colorLen, [0.08, 0.8], colors{nC});
        else
            tmp = rand(1, 3);
            cmaps{nC} = vivid(colorLen, [0.08, 0.8], tmp);
        end
    end
    
else
    
    %colorLen = 64;
    frame_colors = [166, 206, 227; 31,120,180; 178,223,138; 51,160,44; 251,154,153;
        227,26,28; 253,191,111; 255,127,0; 202,178,214; 106,61,154; 255,255,153; 177,89,40];
    frame_colors = frame_colors/256;
    frame_colors = [frame_colors(2:2:end,:);frame_colors(1:2:end,:)];
    frame_colors = num2cell(frame_colors, 2)';
    colors = {[1, 1, 0], ... % yellow
        [1, 0, 1], ... % Magenta
        [1, 0, 0], ... % red
        [0, 1, 0], ... % Green
        [0, 1, 1], ... %cyan
        [0, 0, 1], ... % blue
        [1, .5, 0], ... % orange
        [1, 0, .5], ... % pink
        [.5, 0, 1], ... % violet
        [0.5, 1, 0], ... % Lime green
        [0, 1, 0.5], ...
        [0, .5, 1]};
    colors = [colors, frame_colors];
    
    for nC = 1:length(cmaps)
        if (nC <= length(colors))
            tmp = colors{nC};
        else
            tmp = rand(1, 3);
        end
        cmaps{nC} = repmat(tmp, colorLen, 1);
    end
    
end

cmap = cat(1, cmaps{:});
cmap = [cmap; gray(colorLen)];


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
myLabel = get(hObject,'Label');
if (strcmpi(myLabel, 'composite'))
    coords = handles_data.maxVoxelPos;
else
    chk = strmatch(lower(myLabel), lower(handles_data.compMenuLabels), 'exact');
    coords = handles_data.maxCoords(chk(1), :);
end
handles_data.coords = coords(1, :);
handles_data.pixelPos = round(handles_data.voxelOrigin + (handles_data.coords./handles_data.VOX));
%handles_data.centre = handles_data.maxVoxelPos;
set(handles, 'userdata', handles_data);
updateSlices(handles_data);


function subH = plotImage(subH, data)
%% Function to plot the image at the specified position
%

imageAxisPos = get(subH, 'position');
image((data), 'parent', subH);
set(subH, 'units', 'normalized');
set(subH, 'tag', 'image_axes');
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
%axis(subH, 'square');
set(subH, 'XTick', []);
set(subH, 'XTickLabel', []);
set(subH, 'YTick', []);
set(subH, 'YTickLabel', []);



function updateSlices(handles_data)
%% Update slices and labels

%realCoords = (icatb_voxel_to_real(handles_data.structVol, handles_data.coords(1, :)));
colorLen = handles_data.colorLen;
realCoords = handles_data.coords(1, :);
cmap = handles_data.cmap;

voxelCoord = handles_data.pixelPos;
%voxelCoord = round(icatb_real_to_voxel(handles_data.structVol, realCoords));
sliceXY = rot90(squeeze(handles_data.data(:,:,voxelCoord(3))));
sliceXZ = rot90(squeeze(handles_data.data(:,voxelCoord(2), :)));
sliceYZ = rot90(squeeze(handles_data.data(voxelCoord(1),:,:)));

sliceXY = reshape(cmap(sliceXY(:), :), [size(sliceXY), 3]);
sliceXZ = reshape(cmap(sliceXZ(:), :), [size(sliceXZ), 3]);
sliceYZ = reshape(cmap(sliceYZ(:), :), [size(sliceYZ), 3]);

%% X
set(findobj(handles_data.allAxes(1), 'Type', 'Image'), 'CData', (sliceYZ));
set(get(handles_data.allAxes(1), 'XLabel'), 'string', ['X = ', num2str(realCoords(1)), ' mm']);

%% Y
set(findobj(handles_data.allAxes(2), 'Type', 'Image'), 'CData', (sliceXZ));
set(get(handles_data.allAxes(2), 'XLabel'), 'string', ['Y = ', num2str(realCoords(2)), ' mm']);

%% Z
set(findobj(handles_data.allAxes(3), 'Type', 'Image'), 'CData', (sliceXY));
set(get(handles_data.allAxes(3), 'XLabel'), 'string', ['Z = ', num2str(realCoords(3)), ' mm']);
set(handles_data.allAxes, 'CLIM', [1, size(cmap, 1)]);


allColorbars = handles_data.allColorbars;
cMax = 0;
for n = 1:length(allColorbars)
    cMin = cMax + 1;
    cMax = cMax + colorLen;
    set(allColorbars(n), 'Xlim', [cMin, cMax]);
    set(allColorbars(n), 'Xtick', []);
    set(allColorbars(n), 'Ytick', []);
end


function [realWorldPos, voxelCoord] = getCoord(V, real_world_coords, pixelPos)
% Get real world and voxel coordinates

realWorldPos = squeeze(real_world_coords(pixelPos(1), pixelPos(2), pixelPos(3), :));
realWorldPos = realWorldPos(:)';
voxelCoord = icatb_real_to_voxel(V, realWorldPos);


function MAT = concatDat(dat)
%% Concatenate rendered data
%

nCols = 2;
nLoops = ceil(length(dat)/nCols);

row_size = 0;
col_size = [];
for nLoop = 1:nLoops
    dims = size(dat{2*nLoop-1});
    row_size = row_size + dims(1);
    col_size = max([col_size, dims(2)]);
end

col_size = 2*col_size;

MAT = zeros(row_size, col_size, 3);

er = 0;
loopNum = 0;

for i = 1:nLoops
    
    loopNum = loopNum + 1;
    ec = 0;
    if (loopNum <= length(dat))
        tmp = dat{loopNum};
        
        current_dims = size(tmp);
        
        sr = er + 1;
        er = er + current_dims(1);
        
        sc = ec + round((col_size/2 - current_dims(2))/2)+1;
        ec = sc + current_dims(2)-1;
        
        MAT(sr:er, sc:ec, :) = tmp;
    end
    
    loopNum = loopNum + 1;
    
    if (loopNum <= length(dat))
        tmp = dat{loopNum};
        
        current_dims = size(tmp);
        
        sc = round(col_size/2 + (col_size/2 - current_dims(2))/2)+1;
        ec = sc + current_dims(2)-1;
        
        MAT(sr:er, sc:ec, :) = tmp;
    end
    
end



function data = stackData(slices, minVal)
%% Stack slices in columns
%

% tmp = slices{1}; tmp = [tmp;minVal*ones(3, size(tmp, 2))];
% slices{1} = tmp;
%
% tmp = slices{2}; tmp = [tmp;minVal*ones(3, size(tmp, 2))];
% slices{2} = tmp;
%
% tmp = slices{3}; tmp = [tmp;minVal*ones(3, size(tmp, 2))];
% slices{3} = tmp;

ddN = cellfun(@size, slices, 'UniformOutput', false);
ddN = cat(1, ddN{:});
maxSizeX = max(ddN(:, 1));
maxSizeY = max(ddN(:, 2));

% m1 = size(slices{1}, 1); n1 = size(slices{1}, 2);
% m2 = size(slices{2}, 1); n2 = size(slices{2}, 2);
% m3 = size(slices{3}, 1); n3 = size(slices{3}, 2);

%maxSizeX = max([m1, m2, m3]);
%maxSizeY = max([n1, n2, n3]);

data = minVal*ones(maxSizeX, maxSizeY, 3);


e = 0;
for nS = 1:length(slices)
    tmp = slices{nS};
    xinda = ceil(maxSizeX/2) - ceil(size(tmp, 1)/2) + 1;
    xindb = xinda + size(tmp, 1) - 1;
    s = e + 1;
    e = e + size(tmp, 2);
    % inda = ceil(maxSizeY/2) - ceil(size(tmp, 2)/2) + 1;
    % indb = inda + size(tmp, 2) - 1;
    data(xinda:xindb, s:e, :) = tmp;
end


function [stackedData, labels] = returnStackedData(handles_data, coords)

stackedData = {};

imagesCmap = handles_data.cmap;

coordsN = coords;
for n = 1:size(coords, 2)
    [dd, xInds] = unique(coordsN(:, n));
    xInds = sort(xInds);
    inds = 1:length(coordsN(:, n));
    inds(xInds) = [];
    coordsN(inds, n) = NaN;
end

labels = [];

for nCoords = 1:size(coords, 1)
    
    pixelPos = round(handles_data.voxelOrigin + (coords(nCoords, :)./handles_data.VOX));
    % Return slices
    sliceXY = rot90(squeeze(handles_data.data(:,:,pixelPos(3))));
    sliceXZ = rot90(squeeze(handles_data.data(:,pixelPos(2), :)));
    sliceYZ = rot90(squeeze(handles_data.data(pixelPos(1),:,:)));
    
    sliceXY = reshape(imagesCmap(sliceXY(:), :), [size(sliceXY), 3]);
    sliceXZ = reshape(imagesCmap(sliceXZ(:), :), [size(sliceXZ), 3]);
    sliceYZ = reshape(imagesCmap(sliceYZ(:), :), [size(sliceYZ), 3]);
    
    if ~isnan(coordsN(nCoords, 1))
        stackedData{end + 1} = sliceYZ;
        labels = [labels, ' X = ', num2str(coords(nCoords, 1), '%0.1f')];
    end
    
    
    if ~isnan(coordsN(nCoords, 2))
        stackedData{end + 1} = sliceXZ;
        labels = [labels, ' Y = ', num2str(coords(nCoords, 2), '%0.1f')];
    end
    
    if ~isnan(coordsN(nCoords, 3))
        stackedData{end + 1} = sliceXY;
        labels = [labels, ' Z = ', num2str(coords(nCoords, 3), '%0.1f')];
    end
    
end

stackedData = stackData(stackedData, 0);