function icatb_groupNetworks(file_names, network_names, network_vals, varargin)
%% Group components by network names and display slices as a ortho plot.
% 1. file_names - File names.
% 2. structFile - Anatomical file.
% 3. netowrk_names - Network names in a cell array.
% 4. network_vals - File numbers for each network in a cell array.
% varargin - arguments must be passed in pairs
%   a. image_values - Options are 'positive and negative', 'positive',
%   'absolute value' and 'negative'.
%   b. threshold - Threshold
%   c. structfile - Anatomical file to overlay
%   d. prefix - Prefix for plotting title. Default is IC
%   e. convert_to_z - Options are 'yes' and 'no'.
%

icatb_defaults;
global UI_FS;
global IMAGE_VALUES;
global CONVERT_Z;
global THRESHOLD_VALUE;
global FONT_COLOR;
global BG_COLOR;

%% Initialise vars
image_values = 'positive';
thresholds = 1;
structFile = fullfile(fileparts(which('gift.m')), 'icatb_templates', 'ch2bet.nii');
convert_to_z = 'yes';
prefix = 'IC';
interMediatePlot = 1;
colorbarPlot = 0;

%% Parse arguments
for i = 1:2:length(varargin)
    if (strcmpi(varargin{i}, 'image_values'))
        image_values = varargin{i + 1};
    elseif (strcmpi(varargin{i}, 'threshold'))
        thresholds = varargin{i + 1};
    elseif (strcmpi(varargin{i}, 'structfile'))
        structFile = varargin{i + 1};
    elseif (strcmpi(varargin{i}, 'prefix'))
        prefix = varargin{i + 1};
    elseif (strcmpi(varargin{i}, 'convert_to_z'))
        convert_to_z = varargin{i + 1};
    elseif (strcmpi(varargin{i}, 'interMediatePlot'))
        interMediatePlot = varargin{i + 1};
    elseif (strcmpi(varargin{i}, 'colorbar'))
        colorbarPlot = varargin{i + 1};
    end
end

useGUI = 0;
if (~exist('file_names', 'var'))
    useGUI = 1;
    file_names = icatb_selectEntry('typeSelection', 'multiple', 'typeEntity', 'file', 'filter', '*comp*.img;*comp*.nii', 'title', 'Select components', 'filetype', 'image');
    drawnow;
end

if (isempty(file_names))
    error('Files are not selected');
end

file_names = icatb_rename_4d_file(file_names);

compGroupNames = '';

if (useGUI)
    
    figData.numICs = size(file_names, 1); %% Draw graphics
    figData.comp = [];
    figureTag = 'setup_comp_networks_gui';
    figHandle = findobj('tag', figureTag);
    if (~isempty(figHandle))
        delete(figHandle);
    end
    
    % Setup figure for GUI
    InputHandle = icatb_getGraphics('Group networks', 'displaygui', figureTag, 'on');
    set(InputHandle, 'menubar', 'none');
    set(InputHandle, 'userdata', figData);
    
    promptWidth = 0.52; controlWidth = 0.32; promptHeight = 0.05;
    xOffset = 0.02; yOffset = promptHeight; listboxHeight = 0.18; listboxWidth = controlWidth; yPos = 0.9;
    okWidth = 0.12; okHeight = promptHeight;
    
    promptPos = [xOffset, yPos - 0.5*yOffset, promptWidth, promptHeight];
    
    compNetworkNameData = getCompData(file_names);
    drawnow;
    
    
    
    listboxXOrigin = promptPos(1) + promptPos(3) + xOffset;
    listboxYOrigin = promptPos(2) - 0.5*listboxHeight;
    
    %%  Components listbox (Group components by name)
    %promptPos(2) = listboxYOrigin - yOffset - 0.5*listboxHeight;
    textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Add Components', 'tag', ...
        'prompt_components', 'fontsize', UI_FS - 1);
    icatb_wrapStaticText(textH);
    listboxYOrigin = promptPos(2) - 0.5*listboxHeight;
    listboxPos = [listboxXOrigin, listboxYOrigin, listboxWidth, listboxHeight];
    listCompH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'listbox', 'position', listboxPos, 'string', compGroupNames, 'tag', ...
        'comp', 'fontsize', UI_FS - 1, 'min', 0, 'max', 2, 'value', 1, 'callback', {@addCompNetwork, InputHandle}, 'userdata', compNetworkNameData);
    
    addButtonPos = [listboxPos(1) + listboxPos(3) + xOffset, listboxPos(2) + 0.5*listboxPos(4) + 0.5*promptHeight, promptHeight + 0.01, promptHeight - 0.01];
    removeButtonPos = [listboxPos(1) + listboxPos(3) + xOffset, listboxPos(2) + 0.5*listboxPos(4) - 0.5*promptHeight, promptHeight + 0.01, promptHeight - 0.01];
    icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', addButtonPos, 'string', '+', 'tag', 'add_comp_button', 'fontsize',...
        UI_FS - 1, 'callback', {@addCompNetwork, InputHandle});
    icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', removeButtonPos, 'string', '-', 'tag', 'remove_comp_button', 'fontsize',...
        UI_FS - 1, 'callback', {@removeCompNetwork, InputHandle});
    
    
    promptPos(2) = listboxYOrigin - yOffset - 0.5*promptHeight;
    promptPos(3) = promptWidth;
    textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Enter prefix', 'tag', 'prompt_prefix', 'fontsize', UI_FS - 1);
    icatb_wrapStaticText(textH);
    
    editPos = promptPos;
    editPos(1) = editPos(1) + editPos(3) + xOffset;
    editPos(3) = controlWidth;
    editH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'edit', 'position', editPos, 'string', 'IC', 'tag', 'prefix', 'fontsize', UI_FS - 1);
    
    promptPos(2) = promptPos(2) - yOffset - 0.5*promptHeight;
    promptPos(3) = promptWidth;
    textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Convert To Z-scores?', 'tag', 'prompt_z', 'fontsize', UI_FS - 1);
    icatb_wrapStaticText(textH);
    
    editPos = promptPos;
    editPos(1) = editPos(1) + editPos(3) + xOffset;
    editPos(3) = controlWidth;
    zOptions = char('Yes', 'No');
    chkInd = strmatch(lower(CONVERT_Z), lower(zOptions), 'exact');
    popupH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'popup', 'position', editPos, 'string', zOptions, 'tag', 'convert_to_z', 'fontsize', UI_FS - 1, ...
        'value', chkInd);
    
    
    promptPos(2) = promptPos(2) - yOffset - 0.5*promptHeight;
    promptPos(3) = promptWidth;
    textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Enter threshold', 'tag', 'prompt_threshold', 'fontsize', ...
        UI_FS - 1);
    icatb_wrapStaticText(textH);
    
    editPos = promptPos;
    editPos(1) = editPos(1) + editPos(3) + xOffset;
    editPos(3) = controlWidth;
    editH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'edit', 'position', editPos, 'string', THRESHOLD_VALUE, 'tag', 'threshold', 'fontsize', UI_FS - 1);
    
    
    promptPos(2) = promptPos(2) - yOffset - 0.5*promptHeight;
    promptPos(3) = promptWidth;
    textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Select Image Values', 'tag', 'prompt_display', ...
        'fontsize', UI_FS - 1);
    icatb_wrapStaticText(textH);
    
    imOptions = char('Positive and Negative', 'Positive', 'Absolute value', 'Negative');
    returnValue = strmatch(lower(IMAGE_VALUES), lower(imOptions), 'exact');
    
    editPos = promptPos;
    editPos(1) = editPos(1) + editPos(3) + xOffset;
    editPos(3) = controlWidth;
    popupH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'popup', 'position', editPos, 'string', imOptions, 'tag', 'image_values', 'fontsize', UI_FS - 1, ...
        'value', returnValue);
    
    anatWidth = 0.2;
    anatHeight = 0.05;
    anatPos = [0.25 - 0.5*okWidth, promptPos(2) - yOffset - 0.5*okHeight, anatWidth, anatHeight];
    anatPos(2) = anatPos(2) - 0.5*anatPos(4);
    icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', anatPos, 'string', 'Select Anatomical', 'tag', 'anat_button', 'fontsize',...
        UI_FS - 1, 'callback', {@selectAnatomical, InputHandle});
    
    
    %% Add cancel, save and run buttons
    okPos = [0.75 - 0.5*okWidth, promptPos(2) - yOffset - 0.5*okHeight, okWidth, okHeight];
    okPos(2) = okPos(2) - 0.5*okPos(4);
    icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', okPos, 'string', 'Ok', 'tag', 'run_button', 'fontsize',...
        UI_FS - 1, 'callback', {@runCallback, InputHandle});
    
    waitfor(InputHandle);
    
    
    % check application data
    if isappdata(0, 'inputDispData')
        % get the application data
        answers = getappdata(0, 'inputDispData');
        rmappdata(0, 'inputDispData'); % remove the application data
        network_names = answers.network_names;
        network_vals = answers.network_vals;
        thresholds = answers.threshold;
        image_values = answers.image_values;
        convert_to_z = answers.convert_to_z;
        prefix = answers.prefix;
        try
            structFile = answers.structFile;
        catch
        end
    end
    
    
end

drawnow;


if (~exist('network_names', 'var') || isempty(network_names))
    error('Network names doesn''t exist or is empty.');
end

if (~exist('network_vals', 'var') || isempty(network_vals))
    error('Network values doesn''t exist or is empty.');
end

if (length(network_names) ~= length(network_vals))
    error('Network names must match the network vals');
end

if ((size(thresholds, 2) == 2) && (size(thresholds, 1) == 1))
    thresholds = repmat(thresholds, length(network_names), 1);
end

if (length(thresholds) == 1)
    thresholds = thresholds.*ones(length(network_names), 1);
end

maxCol = [];
for nR = 1:length(network_names)
    maxCol = max([length(network_vals{nR}), maxCol]);
end

% num_rows = ceil(sqrt(length(network_names)));
% num_cols = ceil(length(network_names)/num_rows);

if (interMediatePlot)
    graphicsH = repmat(struct('H', []), 1, length(network_names) + 1);
else
    graphicsH(1).H = [];
end
sz = get(0, 'ScreenSize');

allStrs = cell(length(network_names), 1);
Ndat = allStrs;

figPosition = [25, 25, sz(3) - 50, sz(4) - 50];

for N = 1:length(network_names)
    
    if (interMediatePlot)
        graphicsH(N).H = icatb_getGraphics('Networks', 'graphics', 'group_networks', 'on');
        set(graphicsH(N).H , 'resize', 'on');
        set(graphicsH(N).H , 'position', figPosition);
    end
    
    %     sz = get(0, 'screensize');
    %     sz(1) = 50;
    %     sz(2) = 50;
    %     sz(3) = sz(3) - 2*sz(1);
    %     sz(4) = sz(4) - 2*sz(2);
    %     set(graphicsH(N).H, 'position', sz);
    
    comps = network_vals{N};
    
    %subInds = getInds(length(comps), maxCol);
    tmpDat = [];
    tmpStrs = cell(length(comps), 1);
    for nComp = 1:length(comps)
        
        %sh = subplot(1, maxCol, subInds(nComp));
        if (interMediatePlot)
            sh = subplot(1, maxCol, nComp);
            set(sh, 'color', BG_COLOR);
            set(sh, 'Xcolor', FONT_COLOR);
            set(sh, 'Ycolor', FONT_COLOR);
        end
        
        currentFile = deblank(file_names(comps(nComp), :));
        hD = icatb_orth_views(currentFile, 'structfile', structFile, 'image_values', image_values, 'convert_to_zscores', convert_to_z, 'labels', ...
            [prefix, ' ', num2str(comps(nComp))], 'get_interp_data', 1, 'set_to_max_voxel', 1, 'threshold', thresholds(N, :));
        [sliceXY, sliceXZ, sliceYZ] = returnSlices(hD.data, hD.maxVoxelPos);
        sliceYZ = rot90(sliceYZ);
        sliceXZ = rot90(sliceXZ);
        sliceXY = rot90(sliceXY);
        tmp = stackData({sliceYZ, sliceXZ, sliceXY}, 101);
        
        realCoords = (hD.maxVoxelPos - hD.voxelOrigin).*hD.VOX;
        
        if (interMediatePlot)
            imagesc(tmp, [1, 200]);
            colormap(hD.cmap);
            axis(sh, 'image');
            set(sh, 'Ytick', []);
            set(sh, 'Xtick', []);
            if (nComp == 1)
                text(-0.1, 0.5, 0, network_names{N}, 'rotation', 90, 'parent', sh, 'units', 'normalized');
            end
        end
        
        tmpStr = ['IC',  num2str(comps(nComp)), ' (', num2str(realCoords(1)), ',', num2str(realCoords(2)), ',', num2str(realCoords(3)), ')'];
        tmpStrs{nComp} = tmpStr;
        %tmpStrs{nComp, 1} = ['IC',  num2str(comps(nComp))];
        %tmpStrs{nComp, 2} = ['(', num2str(realCoords(1)), ',', num2str(realCoords(2)), ',', num2str(realCoords(3)), ')'];
        str = charC (['IC',  num2str(comps(nComp))], ['(', num2str(realCoords(1)), ',', num2str(realCoords(2)), ',', num2str(realCoords(3)), ')']);
        if (interMediatePlot)
            title(str, 'parent', sh, 'horizontalalignment', 'center', 'fontsize', 9);
            
            axis(sh, 'off');
        end
        
        tmpDat = [tmpDat, tmp];
        
    end
    
    Ndat{N} = tmpDat;
    allStrs{N} = char(tmpStrs);
    
end


graphicsH(end).H = icatb_getGraphics('Networks', 'graphics', 'group_networks', 'on');
set(graphicsH(end).H , 'resize', 'on');
set(graphicsH(end).H , 'position', figPosition);

lengths = cellfun(@length, network_vals, 'UniformOutput', false);
lengths = cellfun(@length, network_vals, 'UniformOutput', false);
[maxLength, maxRInds] = max([lengths{:}]);

nrows = ceil(sqrt(length(network_names)));
ncols = ceil(length(network_names) / nrows);
colormap(hD.cmap);
dd=cellfun(@size, Ndat, 'UniformOutput', false);
dd = cat(1, dd{:});
dims = max(dd);
midSize = ceil(dims/2);
count = 0;
yPos = 0.97;
yOffset = 0.08;
xOffset = 0.062;
axesWidth = ((yPos - ncols*xOffset) / ncols);
axesHeight = ((yPos - nrows*yOffset) / nrows);
axesWidth = min([axesWidth, axesHeight]);
axesHeight = axesWidth;

for nr = 1:nrows
    
    XOrigin = 0.05;
    YOrigin = yPos - nr*axesHeight - nr*yOffset;
    for nc = 1:ncols
        count = count + 1;
        if (count > length(Ndat))
            break;
        end
        tmp = 101*ones(dims);
        tmp2 = Ndat{count};
        %tmp(1:size(tmp2, 1), 1:size(tmp2, 2)) = tmp2;
        midTmp = ceil(size(tmp2)/2);
        org = midSize - midTmp;
        tmp(org(1) + 1: org(1) + size(tmp2, 1), org(2) + 1: org(2) + size(tmp2, 2)) = tmp2;
        %tmp(1:size(tmp2, 1), 1:size(tmp2, 2)) = tmp2;
        sh = axes('units', 'normalized', 'position', [XOrigin, YOrigin, axesWidth, axesHeight]);
        %sh = subplot(nrows, ncols, count);
        %set(sh, 'units', 'normalized');
        %axesPos = get(sh, 'position');
        imagesc(tmp, [1, 200]);
        axis(sh, 'image');
        currentStr = allStrs{count};
        %for nF = 1:size(currentStr); dd2{nF} = cat(2, currentStr{nF, 1}, currentStr{nF, 2});end;
        %         offset = 0.1;
        %         compWidth = (1 - offset*maxLength) /maxLength;
        %         xPos = -0.001;
        %         for nComps = 1:size(currentStr, 1)
        %             text(xPos, 1.2, char(currentStr{nComps, 1}, currentStr{nComps, 2}), 'units', 'normalized', 'color', FONT_COLOR, 'fontsize', 6, 'HorizontalAlignment', 'center');
        %             xPos = xPos + offset + compWidth;
        %         end
        title(currentStr, 'fontsize', 7.5, 'HorizontalAlignment', 'center');
        ylabel(network_names{count}, 'fontsize', 11, 'color', FONT_COLOR);
        set(sh, 'color', BG_COLOR);
        set(sh, 'Xcolor', BG_COLOR);
        set(sh, 'Ycolor', BG_COLOR);
        set(sh, 'XTick', []);
        set(sh, 'YTick', []);
        set(sh, 'box', 'off');
        XOrigin = XOrigin + axesWidth + xOffset;
    end
end


if (colorbarPlot && (size(thresholds, 2) == 2))
    ch = colorbar('peer', sh);
    set(ch, 'units', 'normalized');
    colorbarPos = get(ch, 'position');
    colorbarPos(1) = colorbarPos(1) + colorbarPos(3) + 0.03;
    set(ch, 'position', colorbarPos);
    %ch = colorbar(pos);
    set(ch, 'yLim', [1, 100]);
    set(ch, 'xTick', []);
    set(ch, 'yTick', [1, 100]);
    set(ch, 'yTicklabel', cellstr(num2str([hD.minICAIM; hD.maxICAIM], '%0.1f')));
    set(ch, 'XColor', FONT_COLOR, 'YColor', FONT_COLOR);
    set(ch, 'tag', 'colorbar');
    set(ch, 'units', 'normalized');
end

if (interMediatePlot)
    icatb_plotNextPreviousExitButtons(graphicsH);
end





function [sliceXY, sliceXZ, sliceYZ] = returnSlices(data, voxelcoords)
%% Slices (XY, XZ, YZ)
%

sliceXY = reshape(data(:, :, voxelcoords(end)), size(data, 1), size(data, 2));
sliceXZ = reshape(data(:, voxelcoords(2), :), size(data, 1),size(data, 3));
sliceYZ = reshape(data(voxelcoords(1), :, :), size(data, 2), size(data, 3));

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

[m1, n1] = size(slices{1});
[m2, n2] = size(slices{2});
[m3, n3] = size(slices{3});

maxSizeX = max([m1, m2, m3]);
maxSizeY = max([n1, n2, n3]);

%data = minVal*ones(maxSizeX, maxSizeY);
data = minVal*ones(m1 + m2 + m3, maxSizeY);

e = 0;
for nS = 1:length(slices)
    tmp = slices{nS};
    s = e + 1;
    e = e + size(tmp, 1);
    inda = ceil(maxSizeY/2) - ceil(size(tmp, 2)/2) + 1;
    indb = inda + size(tmp, 2) - 1;
    data(s:e, inda:indb) = tmp;
end

function c = charC(a, b)
%% Character array. Strings are centered.
%

len = max([length(a), length(b)]);

c = repmat(' ', 2, len);

s = ceil(len/2 - length(a)/2) + 1;
c(1, s:s+length(a) - 1) = a;
s = ceil(len/2 - length(b)/2) + 1;
c(2, s:s + length(b) - 1) = b;


function compNetworkNameData = getCompData(file_names)

structFile = deblank(file_names(1, :));

%% Get colormap associated with the image values
structData2 =  icatb_spm_read_vols(icatb_spm_vol(structFile));
structData2(isfinite(structData2) == 0) = 0;
structDIM = [size(structData2, 1), size(structData2, 2), 1];

for nC = 1:size(file_names, 1)
    
    tmp = icatb_spm_read_vols(icatb_spm_vol(file_names(nC, :)));
    tmp(isfinite(tmp)==0) = 0;
    
    tmp(tmp ~= 0) = detrend(tmp(tmp ~= 0), 0) ./ std(tmp(tmp ~= 0));
    tmp(abs(tmp) < 1.0) = 0;
    
    if (nC == 1)
        compData = zeros(size(tmp, 1), size(tmp, 2), size(file_names, 1));
    end
    
    [dd, inds] = max(tmp(:));
    
    [x, y, z] = ind2sub(size(tmp), inds);
    
    [tmp, maxICAIM, minICAIM, minInterval, maxInterval] = icatb_overlayImages(reshape(tmp(:, :, z), [1, size(tmp, 1), size(tmp, 2), 1]), structData2(:, :, z), structDIM, structDIM, 1);
    
    compData(:, :, nC) = reshape(tmp, structDIM);
    
end


clim = [minInterval, 2*maxInterval];
cmap = icatb_getColormap(1, 1, 1);

compNetworkNameData.clim = clim;
compNetworkNameData.cmap = cmap;
compNetworkNameData.compData = compData;


function addCompNetwork(hObject, event_data, figH)
%% Add Component network
%

icatb_defaults;
global UI_FS;

figureTag = 'add_comp_dfnc';
compFigHandle = findobj('tag', figureTag);
if (~isempty(compFigHandle))
    delete(compFigHandle);
end

figData = get(figH, 'userdata');

listH = findobj(figH, 'tag', 'comp');

compVals = [];
networkName = '';
if (listH == hObject)
    if (~strcmpi(get(figH, 'selectionType'), 'open'))
        return;
    end
    val = get(listH, 'value');
    try
        networkName = figData.comp(val).name;
        compVals = figData.comp(val).value;
    catch
    end
end

compStr = num2str((1:figData.numICs)');

compFigHandle = icatb_getGraphics('Select Component Networks', 'normal', figureTag);
set(compFigHandle, 'menubar', 'none');

promptWidth = 0.52; controlWidth = 0.32; promptHeight = 0.05;
xOffset = 0.02; yOffset = promptHeight; listboxHeight = 0.6; yPos = 0.9;
okWidth = 0.12; okHeight = promptHeight;

%% Features text and listbox
promptPos = [xOffset, yPos - 0.5*yOffset, promptWidth, promptHeight];
textH = icatb_uicontrol('parent', compFigHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Enter Network Name', 'fontsize', UI_FS - 1);
icatb_wrapStaticText(textH);

promptPos = get(textH, 'position');

editPos = promptPos;
editPos(1) = editPos(1) + editPos(3) + xOffset;
editPos(3) = controlWidth;
icatb_uicontrol('parent', compFigHandle, 'units', 'normalized', 'style', 'edit', 'position', editPos, 'string', networkName, 'tag', 'comp_network_name', 'fontsize', UI_FS - 1);

%% Right Listbox
listbox2Wdith = 0.1;
promptPos(2) = promptPos(2) - yOffset - 0.5*promptHeight;
promptPos(1) = (1 - xOffset - 2*listbox2Wdith);
promptPos(3) = 2*listbox2Wdith;
textH = icatb_uicontrol('parent', compFigHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Components', 'fontsize', UI_FS - 1);
icatb_wrapStaticText(textH);
promptPos = get(textH, 'position');
listboxYOrigin = promptPos(2) - 0.5*yOffset - listboxHeight;
listboxXOrigin = promptPos(1) + 0.5*listbox2Wdith;
listboxPos = [listboxXOrigin, listboxYOrigin, listbox2Wdith, listboxHeight];
compListH = icatb_uicontrol('parent', compFigHandle, 'units', 'normalized', 'style', 'listbox', 'position', listboxPos, 'string', compStr, 'tag', 'components', 'fontsize', UI_FS - 1, ...
    'min', 0, 'max', 2, 'value', compVals);

%% Show components
showWidth = 0.08; showHeight = 0.04;
showButtonPos = [listboxXOrigin + 0.5*listbox2Wdith - 0.5*showWidth, listboxYOrigin - yOffset - 0.5*showHeight, showWidth, showHeight];
showH = icatb_uicontrol('parent', compFigHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', showButtonPos, 'string', 'Show', 'fontsize', UI_FS - 1, 'callback', ...
    {@drawComp, figH, compFigHandle});

%% Plot image on the left hand side
axesPos = [xOffset, listboxYOrigin, listboxHeight, listboxHeight];
axesH = axes('parent', compFigHandle, 'units', 'normalized', 'position', axesPos, 'tag', 'axes_display_comp');

promptPos = axesPos;

%% Add cancel and run buttons
okPos = [0.5 - 0.5*okWidth, promptPos(2) - yOffset - 0.5*okHeight, okWidth, okHeight];
okPos(2) = okPos(2) - 0.5*okPos(4);
icatb_uicontrol('parent', compFigHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', okPos, 'string', 'Done', 'tag', 'done_button', 'fontsize', UI_FS - 1, 'callback', ...
    {@setCompCallback, compFigHandle, figH});

%% Draw components on the left hand side
drawComp(compListH, [], figH, compFigHandle);


function removeCompNetwork(hObject, event_data, figH)
%% Remove Component network
%

figData = get(figH, 'userdata');
listH = findobj(figH, 'tag', 'comp');
val = get(listH, 'value');

if (~isempty(val))
    check = icatb_questionDialog('title', 'Remove Component Network', 'textbody', 'Do you want to remove the component network from the list?');
    if (~check)
        return;
    end
end

try
    strs = cellstr(char(figData.comp.name));
    figData.comp(val) = [];
    strs(val) = [];
    set(listH, 'value', 1);
    set(listH, 'string', strs);
    set(figH, 'userdata', figData);
catch
end

function drawComp(hObject, event_data, figH, compFigHandle)
%% Draw component

icatb_defaults;
global UI_FONTNAME;
global FONT_COLOR;

fontSizeText = 8;
set(compFigHandle, 'pointer', 'watch');

listH = findobj(figH, 'tag', 'comp');
compNetworkNameData = get(listH, 'userdata');

axesH = get(compFigHandle, 'currentaxes');

clim = compNetworkNameData.clim;
cmap = compNetworkNameData.cmap;
compData = compNetworkNameData.compData;

sel_comp = get(findobj(compFigHandle, 'tag', 'components'), 'value');

if (~isempty(sel_comp))
    DIM = [size(compData, 1), size(compData, 2), length(sel_comp)];
    [im, numImagesX, numImagesY, textToPlot] = icatb_returnMontage(compData(:, :, sel_comp), [], DIM, [1, 1, 1], sel_comp);
    image(im, 'parent', axesH, 'CDataMapping', 'scaled');
    set(axesH, 'clim', clim); % set the axis positions to the specified
    axis(axesH, 'off');
    axis(axesH, 'image');
    colormap(cmap);
    textCount = 0;
    dim = size(im);
    yPos = 1 + dim(1) / numImagesY;
    for nTextRows = 1:numImagesY
        xPos = 1;
        for nTextCols = 1:numImagesX
            textCount = textCount + 1;
            if textCount <= DIM(3)
                text(xPos, yPos, num2str(round(textToPlot(textCount))), 'color', FONT_COLOR,  ...
                    'fontsize', fontSizeText, 'HorizontalAlignment', 'left', 'verticalalignment', 'bottom', ...
                    'FontName', UI_FONTNAME, 'parent', axesH);
            end
            xPos = xPos + (dim(2) / numImagesX);
        end
        % end for cols
        yPos = yPos + (dim(1) / numImagesY); % update the y position
    end
else
    cla(axesH);
end

set(compFigHandle, 'pointer', 'arrow');

function runCallback(hObject, event_data, handles)
%% Select values
%

figData = get(handles, 'userdata');

%% Prefix
figData.prefix = get(findobj(handles, 'tag', 'prefix'), 'string');

%% Convert to z-scores
convertToZH = findobj(handles, 'tag', 'convert_to_z');
strs = cellstr(get(convertToZH, 'string'));
val = get(convertToZH, 'value');
convert_to_z = lower(strs{val});
figData.convert_to_z = convert_to_z;

%% Threshold
threshH = findobj(handles, 'tag', 'threshold');
figData.threshold = str2num(get(threshH, 'string'));

%% Image values
imageH = findobj(handles, 'tag', 'image_values');
strs = cellstr(get(imageH, 'string'));
val = get(imageH, 'value');
image_values = lower(strs{val});
figData.image_values = image_values;

if (isempty(figData.comp))
    error('Please select components using + button');
end

network_names = cellstr(char(figData.comp.name));
network_values = cell(1, length(network_names));

for n = 1:length(network_values)
    network_values{n} = figData.comp(n).value;
end

%% Network names and values
figData.network_vals = network_values;
figData.network_names = network_names;

setappdata(0, 'inputDispData', figData);

drawnow;

delete(handles);


function setCompCallback(hObject, event_data, compFigH, handles)
%% Get fields from component

figData = get(handles, 'userdata');
networkNameH = findobj(compFigH, 'tag', 'comp_network_name');
networkName = deblank(get(networkNameH, 'string'));

try
    
    if (isempty(networkName))
        error('You must enter a component network name');
    end
    
    listH = findobj(compFigH, 'tag', 'components');
    comps = get(listH, 'value');
    
    if (isempty(comps))
        error('Components are not selected');
    end
    
    if (length(figData.comp) > 0)
        chk = strmatch(lower(networkName), lower(cellstr(char(figData.comp.name))), 'exact');
        if (~isempty(chk))
            ind = chk;
        end
    end
    
    if (~exist('ind', 'var'))
        ind = length(figData.comp) + 1;
    end
    
    %% Set user selected information in figure
    figData.comp(ind).name = networkName;
    figData.comp(ind).value =  comps(:)';
    set(handles, 'userdata', figData);
    compListH = findobj(handles, 'tag', 'comp');
    set(compListH, 'string', cellstr(char(figData.comp.name)));
    delete(compFigH);
    
catch
    icatb_errorDialog(lasterr, 'Component Selection');
end


function selectAnatomical(hObject, event_data, figH)
%% Anatomical callback
%

figData = get(figH, 'userdata');

startPath = fileparts(which('gift.m'));
startPath = fullfile(startPath, 'icatb_templates');

oldDir = pwd;

if (~exist(startPath, 'dir'))
    startPath = pwd;
end

% get the structural file
structFile = icatb_selectEntry('typeEntity', 'file', 'title', 'Select Structural File', 'filter', ...
    '*.img;*.nii', 'fileType', 'image', 'fileNumbers', 1, 'startpath', startPath);

drawnow;

cd(oldDir);

if (~isempty(structFile))
    figData.structFile = structFile;
    set(figH, 'userdata', figData);
end

function subInds = getInds(N, maxCol)

if (N == maxCol)
    
    subInds = 1:maxCol;
    
else
    
    midPoint = ceil(maxCol/2);
    if (mod(N, 2) == 0)
        nLoops = N/2;
    else
        nLoops = round(N/2) - 1;
    end
    if (N == 1)
        subInds = midPoint;
    else
        subInds = [];
        for n = 1:nLoops
            subInds = [subInds, midPoint - n];
            subInds = [subInds, midPoint + n];
        end
        if (length(subInds) < N)
            subInds = [subInds, midPoint];
        end
        subInds = unique(subInds);
    end
    
    
end


function cdata = crop(figH)
%% Crop figure
%

allAxesH = findobj(figH, 'type', 'axes');
set(allAxesH, 'units', 'pixels');
pos = get(allAxesH(1), 'position');

Offset = 20;
axesWidth = pos(3);
axesHeight = pos(4);
xOrigin = pos(1);
yOrigin = pos(2);
maxX = xOrigin;
for nA = 2:length(allAxesH)
    pos = get(allAxesH(nA), 'position');
    xOrigin = min([xOrigin, pos(1)]);
    maxX = max([maxX, pos(1)]);
end

aspectRatio = get(allAxesH(1), 'PlotBoxAspectRatio');
height2 = (aspectRatio(2)/aspectRatio(1))*axesWidth;


newYOrigin = ((axesHeight - height2)/2) + yOrigin-Offset;
newXOrigin = xOrigin - Offset;
height2 = height2 + Offset*3;
width2 = maxX + axesWidth+20 - newXOrigin;

rect = [newXOrigin, newYOrigin, width2, height2];

im = getframe(figH, rect);
cdata = im.cdata;



