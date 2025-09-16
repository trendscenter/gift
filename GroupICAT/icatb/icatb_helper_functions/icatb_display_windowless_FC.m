function icatb_display_windowless_FC(param_file)
%% Display windowless functional connectivity results
%
% Inputs:
% param_file - Windowless functional connectivity parameter file (*windowless_fc.mat)
%
%


icatb_defaults;
global UI_FS;
global THRESHOLD_VALUE;
global CONVERT_Z;
global IMAGE_VALUES;


if (~exist('param_file', 'var'))
    param_file = icatb_selectEntry('typeEntity', 'file', 'typeSelection', 'single', 'filter', '*windowless_fc.mat', 'title', 'Select a valid parameter file');
end


% Setup figure for GUI
InputHandle = icatb_getGraphics('Group networks', 'displaygui', 'Display Windowless FC', 'on');
set(InputHandle, 'menubar', 'none');

% Windowless FC param file
load(param_file);
if (~exist('wfcInfo', 'var'))
    error('Selected file is not a valid windowless FC parameter file');
end


fnc_colorbar_label = 'Corr';
load(wfcInfo.param_file);

file_names = fullfile(fileparts(wfcInfo.param_file), sesInfo.icaOutputFiles(1).ses(1).name);
file_names = icatb_rename_4d_file(file_names);

figData.numICs = size(file_names, 1); %% Draw graphics
figData.comp = [];
MenuOptions = {'fnc_colorbar_label', 'Corr', 'conn_threshold', [], 'display_type', 'slices', 'slice_plane', 'sagittal', 'imWidth', [], 'cmap', 'coldhot', 'CLIM', []};
figData.Options = MenuOptions;

set(InputHandle, 'userdata', figData);

figureTag = 'summary_comp_networks_gui';
figHandle = findobj('tag', figureTag);
if (~isempty(figHandle))
    delete(figHandle);
end


optionsMenuH = uimenu('parent', InputHandle, 'label', 'FNC/Connectogram Defaults', 'callback', {@options_callback, InputHandle});

promptWidth = 0.52; controlWidth = 0.32; promptHeight = 0.06;
xOffset = 0.02; yOffset = promptHeight; listboxHeight = 0.18; listboxWidth = controlWidth; yPos = 0.98;
okWidth = 0.12; okHeight = promptHeight;

promptPos = [xOffset, yPos - 0.5*yOffset, promptWidth, promptHeight];

compNetworkNameData = getCompData(file_names);
drawnow;


compGroupNames = '';

listboxXOrigin = promptPos(1) + promptPos(3) + xOffset;

%%  Components listbox (Group components by name)
%promptPos(2) = listboxYOrigin - yOffset - 0.5*listboxHeight;
promptPos(2) = promptPos(2) - 0.3*listboxHeight - yOffset;
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


promptPos(2) = promptPos(2) - 2.2*yOffset - 0.5*promptHeight;
promptPos(3) = promptWidth;
textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Enter display threshold (Spatial maps)', 'tag', ...
    'prompt_threshold', 'fontsize', UI_FS - 1);
icatb_wrapStaticText(textH);

editPos = promptPos;
editPos(1) = editPos(1) + editPos(3) + xOffset;
editPos(3) = controlWidth;
editH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'edit', 'position', editPos, 'string', THRESHOLD_VALUE, 'tag', 'threshold', 'fontsize', UI_FS - 1);


promptPos(2) = promptPos(2) - yOffset - 0.5*promptHeight;
promptPos(3) = promptWidth;
textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Convert to z-scores', 'tag', 'prompt_z', ...
    'fontsize', UI_FS - 1);
icatb_wrapStaticText(textH);

zOptions = char('Yes', 'No');
returnValue = strmatch(lower(CONVERT_Z), lower(zOptions), 'exact');

editPos = promptPos;
editPos(1) = editPos(1) + editPos(3) + xOffset;
editPos(3) = controlWidth;
popupH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'popup', 'position', editPos, 'string', zOptions, 'tag', 'convert_to_z', 'fontsize', UI_FS - 1, ...
    'value', returnValue);


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


promptPos(2) = promptPos(2) - yOffset - 0.5*promptHeight;
promptPos(3) = 0.6;
promptPos(1) = 0.5 - 0.5*promptPos(3);
textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'checkbox', 'position', promptPos, 'string', 'Display Connectogram', 'tag', 'display_connectogram', ...
    'fontsize', UI_FS - 1, 'value', 0);
icatb_wrapStaticText(textH);

%
% promptPos(2) = promptPos(2) - yOffset - 0.5*promptHeight;
% promptPos(3) = promptWidth;
% textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Enter FNC label (Connectogram and FNC matrix)', 'tag', 'prompt_fnc', ...
%     'fontsize', UI_FS - 1);
% icatb_wrapStaticText(textH);
%
%
% editPos = promptPos;
% editPos(1) = editPos(1) + editPos(3) + xOffset;
% editPos(3) = controlWidth;
% editH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'edit', 'position', editPos, 'string', 'Corr', 'tag', 'fnc_label', 'fontsize', UI_FS - 1);


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


GraphicsHandle = [];

% check application data
if isappdata(0, 'inputDispData')
    
    % get the application data
    answers = getappdata(0, 'inputDispData');
    rmappdata(0, 'inputDispData'); % remove the application data
    answers.file_names = file_names;
    %answers.fnc_matrix_file = fnc_matrix_file;
    network_names = answers.network_names;
    network_vals = answers.network_vals;
    thresholds = answers.threshold;
    image_values = answers.image_values;
    convert_to_z = answers.convert_to_z;
    prefix = answers.prefix;
    fnc_colorbar_label = answers.fnc_colorbar_label;
    try
        structFile = answers.structFile;
    catch
    end
    
    
    if (~exist('network_names', 'var') || isempty(network_names))
        error('Network names doesn''t exist or is empty.');
    end
    
    if (~exist('network_vals', 'var') || isempty(network_vals))
        error('Network values doesn''t exist or is empty.');
    end
    
    if (length(network_names) ~= length(network_vals))
        error('Network names must match the network vals');
    end
    
    comp_network_names = [network_names(:), network_vals(:)];
    
    
    %% FNC correlations: Correlations are visualized in a matrix plot
    display_info.FNC = cell(1, size(wfcInfo.ksvd.Dictionary, 2));
    display_info.comp_network_names = comp_network_names;
    
    for nDict = 1:size(wfcInfo.ksvd.Dictionary, 2)
        
        FNCM = wfcInfo.ksvd.Dictionary(:, nDict);
        
        fncComps = comp_network_names(:, 2);
        fncComps = [fncComps{:}];
        fncComps = fncComps(:);
        
        FNCM = FNCM(fncComps, 1)*FNCM(fncComps, 1)';
        
        
        if (~isempty(FNCM))
            if ((size(FNCM, 1) ~= length(fncComps)) && (size(FNCM, 1) ~= size(file_names, 1)))
                error('FNC matrix dimensions does not match the number of components entered');
            end
            
            if (size(FNCM, 1) == size(file_names, 1))
                % Truncate
                FNCM = FNCM(fncComps, fncComps);
            end
        end
        
        FNCM(isfinite(FNCM) == 0) = 0;
        CLIM = max(abs(FNCM(:)));
        
        
        display_info.FNC{nDict} = FNCM;
        
        network_values = cellfun(@length, comp_network_names(:,2), 'UniformOutput', false);
        network_values = [network_values{:}];
        gH = icatb_getGraphics('FNC Correlations', 'graphics', 'fnc_corrs', 'on');
        colormap(answers.cmap);
        drawnow;
        [FH,AH,CH,IH] = icatb_plot_FNC(FNCM, answers.CLIM, cellstr(num2str(fncComps)), (1:length(fncComps)), gH, ['State #', num2str(nDict), ' ', answers.fnc_colorbar_label], gca, ...
            network_values, comp_network_names(:,1));
        
        GraphicsHandle(end + 1).H = FH;
        
        %% Connectogram view - FNC correlations are shown using bezier curves and thumbnails of spatial maps are shown in a circle. Components within the same network are shown in the same color.
        if (answers.display_connectogram)
            FH = icatb_plot_connectogram([], comp_network_names, 'C', FNCM, 'threshold', thresholds(1), 'image_file_names', file_names, 'colorbar_label', ['State #', num2str(nDict), ' ', answers.fnc_colorbar_label], 'cmap', answers.cmap, ...
                'slice_plane', answers.slice_plane, 'conn_threshold', answers.conn_threshold, 'imwidth', answers.imWidth, 'display_type', answers.display_type, 'CLIM', answers.CLIM);
            GraphicsHandle(end + 1).H = FH;
        end
        
    end
    
else
    
    error('Figure window was quit');
    
end


wfcInfo.display_info = display_info;

save(param_file, 'wfcInfo');

icatb_plotNextPreviousExitButtons(GraphicsHandle);



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



function options_callback(hObject, event_data, handles)
%% Options callback
%

figData = get(handles, 'userdata');
Options = figData.Options;

%imOptions = char('Positive and Negative', 'Positive', 'Absolute value', 'Negative');
numParameters = 1;

% matchedInd = strmatch('image_values', lower(Options(1:2:end)), 'exact');
% image_values = Options{2*matchedInd};
% returnValue = strmatch(lower(image_values), lower(imOptions), 'exact');
% inputText(numParameters).promptString = 'Select Image Values';
% inputText(numParameters).answerString = imOptions;
% inputText(numParameters).uiType = 'popup';
% inputText(numParameters).dataType = 'string';
% inputText(numParameters).tag = 'image_values';
% inputText(numParameters).enable = 'on';
% inputText(numParameters).value = returnValue;
% inputText(numParameters).help = struct('title', 'Image Values', 'string', 'Option is provided to display positive and negative, positive, negative and absolute image values.');
%
% numParameters = numParameters + 1;
% matchedInd = strmatch('convert_to_z', lower(Options(1:2:end)), 'exact');
% convert_to_z = Options{2*matchedInd};
% zOptions = char('Yes', 'No');
% chkInd = strmatch(lower(convert_to_z), lower(zOptions), 'exact');
% inputText(numParameters).promptString = 'Convert Spatial Maps To Z-scores?';
% inputText(numParameters).answerString = zOptions;
% inputText(numParameters).uiType = 'popup';
% inputText(numParameters).dataType = 'string';
% inputText(numParameters).tag = 'convert_to_z';
% inputText(numParameters).enable = 'on';
% inputText(numParameters).value = chkInd;
% inputText(numParameters).help = struct('title', 'Threshold', 'string', 'If you have selected option yes, images are converted to z-scores.');
%
% numParameters = numParameters + 1;
% matchedInd = strmatch('threshold', lower(Options(1:2:end)), 'exact');
% threshold = Options{2*matchedInd};
% inputText(numParameters).promptString = 'Enter Display Threshold Value';
% inputText(numParameters).answerString = num2str(threshold);
% inputText(numParameters).uiType = 'edit';
% inputText(numParameters).dataType = 'numeric';
% inputText(numParameters).tag = 'threshold';
% inputText(numParameters).enable = 'on';
% inputText(numParameters).help = struct('title', 'Threshold', 'string', 'Threshold is applied on the spatial maps.');

%numParameters = numParameters + 1;
matchedInd = strmatch('conn_threshold', lower(Options(1:2:end)), 'exact');
inputText(numParameters).promptString = 'Enter FNC threshold';
inputText(numParameters).answerString = Options{2*matchedInd};
inputText(numParameters).uiType = 'edit';
inputText(numParameters).dataType = 'string';
inputText(numParameters).tag = 'conn_threshold';
inputText(numParameters).enable = 'on';
inputText(numParameters).help = struct('title', 'Threshold', 'string', ...
    'Enter functional connectivity threshold. Values less than the connectivity thresholded are not included in the plot.');

numParameters = numParameters + 1;
matchedInd = strmatch('fnc_colorbar_label', lower(Options(1:2:end)), 'exact');
inputText(numParameters).promptString = 'Enter FNC label (in Colorbar)';
inputText(numParameters).answerString =  Options{2*matchedInd};
inputText(numParameters).uiType = 'edit';
inputText(numParameters).dataType = 'string';
inputText(numParameters).tag = 'fnc_colorbar_label';
inputText(numParameters).enable = 'on';
inputText(numParameters).help = struct('title', 'FNC label', 'string', 'Enter FNC label.');

numParameters = numParameters + 1;
matchedInd = strmatch('display_type', lower(Options(1:2:end)), 'exact');
display_type_opts = {'Slices', 'Render'};
display_type = Options{2*matchedInd};
chkInd = strmatch(lower(display_type), lower(display_type_opts), 'exact');
inputText(numParameters).promptString = 'Select display type';
inputText(numParameters).answerString = display_type_opts;
inputText(numParameters).uiType = 'popup';
inputText(numParameters).dataType = 'string';
inputText(numParameters).tag = 'display_type';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = chkInd;
inputText(numParameters).help = struct('title', 'Display Type', 'string', 'You could select either surface or slices to be plotted in a circle...');

numParameters = numParameters + 1;
sliceOptions = {'Axial', 'Coronal', 'Sagittal'};
matchedInd = strmatch('slice_plane', lower(Options(1:2:end)), 'exact');
slice_plane = Options{2*matchedInd};
chkInd = strmatch(lower(slice_plane), lower(sliceOptions), 'exact');
inputText(numParameters).promptString = 'Select Anatomical Plane';
inputText(numParameters).answerString = sliceOptions;
inputText(numParameters).uiType = 'popup';
inputText(numParameters).dataType = 'string';
inputText(numParameters).tag = 'slice_plane';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = chkInd;
inputText(numParameters).help = struct('title', 'Slice Plane', 'string', ...
    'Thumbnails are plotted using the selected slice plane. If you have selected render option, rendered images are shown.');

numParameters = numParameters + 1;
matchedInd = strmatch('cmap', lower(Options(1:2:end)), 'exact');
inputText(numParameters).promptString = 'Enter colormap for correlations';
inputText(numParameters).answerString = Options{2*matchedInd};
inputText(numParameters).uiType = 'edit';
inputText(numParameters).dataType = 'string';
inputText(numParameters).tag = 'cmap';
inputText(numParameters).enable = 'on';
inputText(numParameters).help = struct('title', 'colormap', 'string', 'You could use any valid command to get the colormap like coldhot, jet, winter.');


numParameters = numParameters + 1;
matchedInd = strmatch('clim', lower(Options(1:2:end)), 'exact');
inputText(numParameters).promptString = 'Enter min and max of correlations.';
inputText(numParameters).answerString = num2str(Options{2*matchedInd});
inputText(numParameters).uiType = 'edit';
inputText(numParameters).dataType = 'numeric';
inputText(numParameters).tag = 'CLIM';
inputText(numParameters).enable = 'on';
inputText(numParameters).help = struct('title', 'Colorbar range', 'string', 'Enter min and max of correlations. You could leave it as empty if max and min are determined using data.');

numParameters = numParameters + 1;
matchedInd = strmatch('imwidth', lower(Options(1:2:end)), 'exact');
inputText(numParameters).promptString = 'Enter image width of the spatial maps';
inputText(numParameters).answerString = num2str(Options{2*matchedInd});
inputText(numParameters).uiType = 'edit';
inputText(numParameters).dataType = 'numeric';
inputText(numParameters).tag = 'imWidth';
inputText(numParameters).enable = 'on';
inputText(numParameters).help = struct('title', 'Image Width', 'string', ...
    'You could provide a specified value like 0.06. If nothing is specified, image width is automatically determined.');


answer = icatb_inputDialog('inputtext', inputText, 'Title', 'Select connectogram options', 'handle_visibility',  'on');

if (~isempty(answer))
    % ICA options with flags and the values corresponding to it
    Options = cell(1, 2*length(answer));
    
    for i = 1:length(answer)
        Options{2*i - 1} = inputText(i).tag;
        Options{2*i} = answer{i};
    end
    figData.Options = Options;
    set(handles, 'userdata', figData);
end

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


function runCallback(hObject, event_data, handles)
%% Select values
%

figData = get(handles, 'userdata');

load('icatb_colors.mat', 'coldhot', 'coldhot_sensitive');

try
    
    for n = 1:2:length(figData.Options)
        figData.(figData.Options{n}) = figData.Options{n+1};
    end
    
    
    cmap = coldhot;
    try
        cmap = eval(figData.cmap);
    catch
    end
    
    if (mod(size(cmap, 1), 64) == 0)
        numToSkip = size(cmap, 1)/64;
        cmap = cmap(1:numToSkip:end, :);
    end
    
    figData.cmap = cmap;
    
catch
    
end

%% Prefix
figData.prefix = get(findobj(handles, 'tag', 'prefix'), 'string');

%% Image values
czH = findobj(handles, 'tag', 'convert_to_z');
strs = cellstr(get(czH, 'string'));
val = get(czH, 'value');
convert_to_z = lower(strs{val});
figData.convert_to_z = convert_to_z;


%% Threshold
threshH = findobj(handles, 'tag', 'threshold');
figData.threshold = str2num(get(threshH, 'string'));

%% FNC label
fnc_colorbar_label = get(findobj(handles, 'tag', 'fnc_label'), 'string');
figData.fnc_colorbar_label = fnc_colorbar_label;

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


displayConnH = findobj(handles, 'tag', 'display_connectogram');
figData.display_connectogram = get(displayConnH, 'value');

setappdata(0, 'inputDispData', figData);

drawnow;

delete(handles);


function SetColors(H)

pos = get(0, 'defaultFigurePosition');
set(H, 'color', 'w');
set(H, 'position', pos);
set(H, 'resize', 'on');
set(H, 'visible', 'on');
titleH = get(findobj(H, 'type', 'axes'), 'title');
if (iscell(titleH))
    % titleH = cell2mat(titleH);
    for nT = 1:length(titleH)
        set(titleH{nT}, 'color', 'k');
    end
else
    set(titleH, 'color', 'k');
end
axesH = findobj(H, 'type', 'axes');
set(axesH, 'YColor', 'k', 'XColor', 'k');
set(axesH, 'fontname', 'times');
set(axesH, 'fontsize', 8);
textH = findobj(H, 'type', 'text');
set(textH, 'color', 'k');
set(textH, 'fontname', 'times');
set(textH, 'fontsize', 8);
try
    C = findobj(H, 'type', 'colorbar');
    set(C, 'color', 'k');
    set(get(C, 'label'), 'color', 'k');
catch
end

