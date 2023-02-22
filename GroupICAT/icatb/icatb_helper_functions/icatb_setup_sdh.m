function icatb_setup_sdh(param_file, comp_network_names, resultsDir)
%% Run spatial dynamics hierarchy. Data is reconstructed using timecourses and spatial maps from each network label and components associated with the network. Clustering algorithms like k-means or k-svd is run on the
% reconstructed data.
%
% Inputs:
% 1. param_file - ICA Parameter file (*ica*param*mat)
% 2. comp_network_names - Component network names like:
%    comp_network_names = {'BG', [10, 12];
%                          'AUD', [15, 16, 17]};
%
% methodInfo - data structure containing method specific details
%   a. name - K-means
%   b. opts - Method specific options like number of states, max
%   iterations, etc.
%
% At the end file name prefix*sdh*mat file is written out with the
% necessary info required for display.
%

%% Initialize variables
icatb_defaults;
global PARAMETER_INFO_MAT_FILE;

filterP = ['*', PARAMETER_INFO_MAT_FILE, '*.mat'];

isGUI = 0;
if (~exist('param_file', 'var') || isempty(param_file))
    isGUI = 1;
    param_file = icatb_selectEntry('typeEntity', 'file', 'title', 'Select ICA Parameter File', 'filter', filterP);
    drawnow;
    if (isempty(param_file))
        error('ICA/SDH parameter file not selected');
    end
    resultsDir = icatb_selectEntry('typeEntity', 'directory', 'title', 'Select Spatial dynamics results directory ...');
    drawnow;
end

if (exist(resultsDir, 'dir') ~= 7)
    resultsDir = pwd;
end

if (isempty(param_file))
    error('ICA parameter file not selected');
end

outDir = fileparts(param_file);
if (isempty(outDir))
    outDir = pwd;
end

load(param_file);

sesInfo.userInput.pwd = outDir;
sesInfo.outputDir = outDir;

if (isGUI)
    inputParams = openGUI(sesInfo);
    comp_network_names = inputParams.comp_network_names;
    feature_params = inputParams.feature_params;
end

%% Initialize sdhInfo
sdhInfo.comp_network_names = comp_network_names;
sdhInfo.outputDir = resultsDir;
sdhInfo.numOfSub = sesInfo.numOfSub;
sdhInfo.numOfSess = sesInfo.numOfSess;
sdhInfo.ica_parameter_file = param_file;
try
    sdhInfo.feature_params = feature_params;
catch
end

%% Save params
fname = fullfile(resultsDir, [sesInfo.userInput.prefix, '_sdh_info.mat']);
disp(['Saving parameters info in ', fname]);
save(fname, 'sdhInfo');
fprintf('\n\n');


drawnow;
icatb_run_sdh(sdhInfo);

%diary('off');

function inputParams = openGUI(sesInfo)
%% Open GUI


icatb_defaults;
global UI_FS;

%% Draw graphics
figureTag = 'setup_sdh_gui';
figHandle = findobj('tag', figureTag);
if (~isempty(figHandle))
    delete(figHandle);
end

figData = struct;

figData.numOfSub = sesInfo.numOfSub;
figData.numOfSess = sesInfo.numOfSess;

if (~isfield(figData, 'comp'))
    figData.comp = [];
end

% Setup figure for GUI
InputHandle = icatb_getGraphics('Spatial Dynamics Hierarchy Setup Analysis', 'normal', figureTag, 'on');
set(InputHandle, 'menubar', 'none');


defaultsmenu = uimenu('parent', InputHandle, 'label', 'Spatial Dynamics Hierarchy Defaults', 'tag', 'sdh_defaults');
set(defaultsmenu, 'callback', {@defaultsCallback, InputHandle});

file_names = fullfile(sesInfo.outputDir, sesInfo.icaOutputFiles(1).ses(1).name);
file_names = icatb_rename_4d_file(file_names);

figData.numICs = size(file_names, 1); %% Draw graphics
figData.comp = [];

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


promptPos(2) = listboxPos(2) - 1.5*yOffset;
textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Enter TR in seconds', 'tag', 'prompt_num_clusters', 'fontsize', UI_FS - 1, 'tooltipstring', ...
    'Enter number of clusters or states. If you want to estimate clusters, leave string as empty');
icatb_wrapStaticText(textH);

editPos = promptPos;
editPos(1) = editPos(1) + editPos(3) + xOffset;
editPos(3) = controlWidth;
editH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'edit', 'position', editPos, 'string', '2', 'tag', 'TR', 'fontsize', UI_FS - 1);


% promptPos(2) = promptPos(2) - 1.5*yOffset;
% textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Select method', 'tag', 'prompt_num_clusters', 'fontsize', UI_FS - 1);
% icatb_wrapStaticText(textH);
%
% popupPos = promptPos;
% popupPos(1) = popupPos(1) + popupPos(3) + xOffset;
% popupPos(3) = controlWidth;
% popupH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'popup', 'position', popupPos, 'string', {'K-means', 'KSVD'}, 'tag', 'method', 'fontsize', UI_FS - 1, 'userdata', 'callback', ...
%     {@optsCallback, InputHandle});


%% Add cancel, save and run buttons
promptPos(2) = promptPos(2) - 1.5*yOffset;
cancelPos = [0.25 - 0.5*okWidth, promptPos(2) - yOffset - 0.5*okHeight, okWidth, okHeight];
cancelPos(2) = cancelPos(2) - 0.5*cancelPos(4);
icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', cancelPos, 'string', 'Cancel', 'tag', 'cancel_button', 'fontsize',...
    UI_FS - 1, 'callback', 'delete(gcbf);');

okPos = [0.75 - 0.5*okWidth, promptPos(2) - yOffset - 0.5*okHeight, okWidth, okHeight];
okPos(2) = okPos(2) - 0.5*okPos(4);
icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', okPos, 'string', 'Run', 'tag', 'run_button', 'fontsize',...
    UI_FS - 1, 'callback', {@applyCallback, InputHandle});


set(InputHandle, 'userdata', figData);
set(InputHandle, 'visible', 'on');
drawnow;

try
    waitfor(InputHandle);
catch
end

appName = 'sdh_app_data';
if (isappdata(0, appName))
    inputParams = getappdata(0, appName);
    rmappdata(0, appName);
else
    error('Input params not seleted');
end

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

figureTag = 'add_comp_sdh';
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


function applyCallback(hObject, event_data, handles)
%% Apply callback
%

figData = get(handles, 'userdata');
if (~isfield(figData, 'feature_params'))
    defaultsCallback(findobj(handles, 'tag', 'sdh_defaults'), [], handles, 'off');
    figData = get(handles, 'userdata');
end

if (~isfield(figData, 'comp') || isempty(figData.comp))
    error('Please select component network names');
end

comps = [figData.comp.value];

if (length(comps) ~= length(unique(comps)))
    error('Please check if there are duplicate entries of the same component in different networks');
end

comp_network_names = cell(length(figData.comp), 2);
for n = 1:size(comp_network_names, 1)
    comp_network_names{n, 1} = figData.comp(n).name;
    comp_network_names{n, 2} = figData.comp(n).value(:)';
end

%figData.methodInfo.num_states = str2num(get(findobj(handles, 'tag', 'num_clusters'), 'string'));
sdhInfo.comp_network_names = comp_network_names;
sdhInfo.feature_params = figData.feature_params;
%sdhInfo.methodInfo = figData.methodInfo;

setappdata(0, 'sdh_app_data', sdhInfo);

delete(gcbf);


function optsCallback(hObject, event_data, handles)
%% Options callback
%


getStr = get(hObject, 'string');
getval = get(hObject, 'value');
method_name = getStr{getval};

dlg_title = 'Select K-means options';
distance_opts = {'City', 'sqEuclidean', 'Hamming', 'Correlation', 'Cosine'};

numParameters = 1;
inputText(numParameters).promptString = 'Enter maximum number of iterations';
inputText(numParameters).uiType = 'edit';
inputText(numParameters).answerString = '150';
inputText(numParameters).dataType = 'numeric';
inputText(numParameters).tag = 'max_iter';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = 1;


if (strcmpi(method_name, 'k-means'))
    
    numParameters = numParameters + 1;
    
    inputText(numParameters).promptString = 'Select distance method';
    inputText(numParameters).uiType = 'popup';
    inputText(numParameters).answerString = distance_opts;
    inputText(numParameters).dataType = 'string';
    inputText(numParameters).tag = 'kmeans_distance_method';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 1;
    
    numParameters = numParameters + 1;
    
    inputText(numParameters).promptString = 'Number of times to repeat the clustering';
    inputText(numParameters).uiType = 'edit';
    inputText(numParameters).answerString = '10';
    inputText(numParameters).dataType = 'numeric';
    inputText(numParameters).tag = 'kmeans_num_replicates';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 1;
    
end


answers = icatb_inputDialog('inputtext', inputText, 'Title', dlg_title, 'handle_visibility', 'on');

opts.max_iter=150;
opts.kmeans_distance_method = 'City';
opts.kmeans_num_replicates = 10;

if (~isempty(answers))
    for n = 1:length(answers)
        opts.(inputText(n).tag) = answers{n};
    end
end


figData = get(handles, 'userdata');
methodInfo.name = method_name;
methodInfo.opts = opts;

figData.methodInfo = methodInfo;

set(handles, 'userdata', figData);


function defaultsCallback(hObject, ed, handles, figVisibility)
%% Defaults callback

if (~exist('figVisibility', 'var'))
    figVisibility = 'on';
end

sdhInfo = get_info_controls(handles);

try
    covInfo = sdhInfo.feature_params.final.tc_covariates;
catch
end

if (~exist('covInfo', 'var'))
    covInfo.numOfDataSets = sdhInfo.numOfSub*sdhInfo.numOfSess;
end

default_params = icatb_sdh_options('covInfo', covInfo);

if (~isfield(sdhInfo, 'feature_params'))
    feature_params = default_params;
else
    feature_params = sdhInfo.feature_params;
end

tags = cellstr(char(feature_params(1).inputParameters(1).options.tag));
motion_inds = strmatch('tc_covariates', tags, 'exact');
if (isempty(motion_inds))
    feature_params(1).inputParameters(1).options(end+1) = feature_params(1).inputParameters(1).options(end);
    feature_params(1).inputParameters(1).options(end).callback = [];
    feature_params(1).inputParameters(1).options(end) = default_params(1).inputParameters(1).options(end);
end


for nF = 1:length(feature_params.inputParameters)
    for nO = 1:length(feature_params.inputParameters(nF).options)
        if (~isfield(feature_params.inputParameters(nF).options(nO), 'enable') || isempty(feature_params.inputParameters(nF).options(nO).enable))
            feature_params.inputParameters(nF).options(nO).enable = 'on';
            feature_params.defaults(nF).options(nO).enable = 'on';
        end
    end
end

% tags = cellstr(char(feature_params(1).inputParameters(2).options.tag));
% check_window_type = strmatch('window_type', tags, 'exact');
% check_window_alpha = strmatch('window_alpha', tags, 'exact');
% if (~isempty(check_window_type))
%     if (feature_params(1).inputParameters(2).options(check_window_type).value == 1)
%         feature_params(1).inputParameters(2).options(check_window_alpha).answerString = '3';
%     end
%     feature_params(1).inputParameters(2).options(check_window_type) = [];
% end

out = icatb_OptionsWindow(feature_params.inputParameters, feature_params.defaults, figVisibility, ...
    'title', 'Spatial Dynamics Hierarchy Options');

drawnow;

feature_params.inputParameters = out.inputParameters;
feature_params.final = out.results;
sdhInfo.feature_params = feature_params;
set(handles, 'userdata', sdhInfo);


function sdhInfo = get_info_controls(handles)
%% Get information from user controls

sdhInfo = get(handles, 'userdata');
TRH = findobj(handles, 'tag', 'TR');
TR = str2num(get(TRH, 'string'));

sdhInfo.TR = TR;

set(handles, 'userdata', sdhInfo);


function tc = regress_cov(tc, file_name, scansToInclude)
%% Regress covariates from timecourses
%

file_name = deblank(file_name);

X = icatb_load_ascii_or_mat(file_name);

if (~exist('scansToInclude', 'var') || isempty(scansToInclude))
    scansToInclude = (1:size(X, 1));
end

scansToInclude(scansToInclude > size(X, 1)) = [];

if (isempty(scansToInclude))
    error(['Please check file numbers specified for file ', file_name]);
end

X = icatb_zscore(X);

X = X(scansToInclude, :);

if (size(X, 1) ~= size(tc, 1))
    error(['Please check the timepoints in file ', file_name]);
end

betas = pinv(X)*tc;

% Remove variance associated with the covariates
tc = tc - X*betas;