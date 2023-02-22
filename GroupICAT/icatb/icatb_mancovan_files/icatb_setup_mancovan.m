function icatb_setup_mancovan(param_file)
%% Setup mancovan
%
% Inputs:
% 1. param_file - ICA Parameter file or Mancovan parameter file
%

icatb_defaults;
global UI_FS;
global PARAMETER_INFO_MAT_FILE;

%% Select Mancovan file
if (~exist('param_file', 'var'))
    param_file = icatb_selectEntry('title', 'Select ICA/Mancovan Parameter File', 'typeEntity', 'file', 'typeSelection', 'single', 'filter', ['*', PARAMETER_INFO_MAT_FILE, '.mat;*mancovan.mat']);
    drawnow;
    if (isempty(param_file))
        error('ICA parameter file is not selected');
    end
end

drawnow;

[inDir, paramF, extn] = fileparts(param_file);
if (isempty(inDir))
    inDir = pwd;
end

param_file = fullfile(inDir, [paramF, extn]);
load(param_file);

if (~exist('mancovanInfo', 'var'))
    error(['Selected file ', param_file, ' is not a valid mancovan file. Please use create design matrix to setup mancovan design']);
end

if (~isfield(mancovanInfo, 'X'))
    mancovanInfo = icatb_mancovan_full_design(mancovanInfo);
end

outputDir = inDir;
mancovanInfo.userInput.outputDir = outputDir;
load(mancovanInfo.userInput.ica_param_file);
if (~exist(deblank(mancovanInfo.userInput.compFiles(1, :)), 'file'))
    unzipFiles(sesInfo, sesInfo.icaOutputFiles(1).ses(1).name, fileparts(mancovanInfo.userInput.ica_param_file));
end

cd(outputDir);

msg = 'Opening Setup Mancovan GUI ...';

disp(msg);

msgH = helpdlg(msg, 'Setup Mancovan');

drawnow;
structFile = icatb_rename_4d_file(deblank(mancovanInfo.userInput.compFiles(1, :)));
structFile = deblank(structFile(1, :));
if (~exist(icatb_parseExtn(structFile), 'file'))
    structFile = deblank(mancovanInfo.userInput.compFiles(1, :));
end

%% Get colormap associated with the image values
structData2 =  icatb_spm_read_vols(icatb_spm_vol(structFile));
structData2(isfinite(structData2) == 0) = 0;
structDIM = [size(structData2, 1), size(structData2, 2), 1];

compNetworkNameData = getCompData(mancovanInfo.userInput.compFiles);

if (mancovanInfo.userInput.numOfSub < 2)
    error('Cannot do stats for one subject');
end

selFeaturesVal = [];
compGroupNames = '';
covNames = '';
allFeatures = {'Spatial Maps', 'Timecourses Spectra', 'FNC Correlations', 'FNC Correlations (lag)'};

if (~isfield(mancovanInfo.userInput, 'features'))
    mancovanInfo.userInput.features = '';
end

selFeatures = mancovanInfo.userInput.features;

[dd, selFeaturesVal] = intersect(lower(cellstr(allFeatures)), lower(cellstr(selFeatures)));

selFeaturesVal = sort(selFeaturesVal);

if (~isfield(mancovanInfo.userInput, 'comp'))
    mancovanInfo.userInput.comp = [];
end

try
    compGroupNames = (cellstr(char(mancovanInfo.userInput.comp.name)));
catch
end

if (~isfield(mancovanInfo.userInput, 'doEstimation') || isempty(mancovanInfo.userInput.doEstimation))
    mancovanInfo.userInput.doEstimation = 0;
end

if (~isfield(mancovanInfo.userInput, 'numOfPCs'))
    mancovanInfo.userInput.numOfPCs = [];
end

if (~isfield(mancovanInfo.userInput, 'p_threshold'))
    mancovanInfo.userInput.p_threshold = 0.01;
end

if (~isfield(mancovanInfo.userInput, 'TR'))
    if (isfield(sesInfo, 'TR'))
        mancovanInfo.userInput.TR = sesInfo.TR;
    else
        mancovanInfo.userInput.TR = 1;
    end
end

clear sesInfo;

%% Draw graphics
figureTag = 'setup_mancovan_gui';
figHandle = findobj('tag', figureTag);
if (~isempty(figHandle))
    delete(figHandle);
end

% Setup figure for GUI
InputHandle = icatb_getGraphics('Mancovan Setup Analysis', 'normal', figureTag, 'off');
set(InputHandle, 'menubar', 'none');
set(InputHandle, 'userdata', mancovanInfo);

defaultsmenu = uimenu('parent', InputHandle, 'label', 'Mancovan-Defaults', 'tag', 'mancovan_defaults');
set(defaultsmenu, 'callback', {@defaultsCallback, InputHandle});

promptWidth = 0.52; controlWidth = 0.32; promptHeight = 0.05;
xOffset = 0.02; yOffset = promptHeight; listboxHeight = 0.18; listboxWidth = controlWidth; yPos = 0.9;
okWidth = 0.12; okHeight = promptHeight;

%% Features text and listbox
promptPos = [xOffset, yPos - 0.5*yOffset, promptWidth, promptHeight];
textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Select Features', 'tag', ...
    'prompt_features', 'fontsize', UI_FS - 1);
icatb_wrapStaticText(textH);

promptPos = get(textH, 'position');

listboxXOrigin = promptPos(1) + promptPos(3) + xOffset;
listboxYOrigin = promptPos(2) - 0.5*listboxHeight;
listboxPos = [listboxXOrigin, listboxYOrigin, listboxWidth, listboxHeight];
icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'listbox', 'position', listboxPos, 'string', ...
    allFeatures, 'tag', 'features', 'fontsize', UI_FS - 1, 'min', 0, 'max', 2, 'value', selFeaturesVal);


%%  Components listbox (Group components by name)
promptPos(2) = listboxYOrigin - yOffset - 0.5*listboxHeight;
textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Add Components', 'tag', ...
    'prompt_components', 'fontsize', UI_FS - 1);
icatb_wrapStaticText(textH);
listboxYOrigin = promptPos(2) - 0.5*listboxHeight;
listboxPos = [listboxXOrigin, listboxYOrigin, listboxWidth, listboxHeight];
listCompH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'listbox', 'position', listboxPos, 'string', compGroupNames, 'tag', ...
    'comp', 'fontsize', UI_FS - 1, 'min', 0, 'max', 1, 'value', 1, 'callback', {@addCompNetwork, InputHandle}, 'userdata', compNetworkNameData);

addButtonPos = [listboxPos(1) + listboxPos(3) + xOffset, listboxPos(2) + 0.5*listboxPos(4) + 0.5*promptHeight, promptHeight + 0.01, promptHeight - 0.01];
removeButtonPos = [listboxPos(1) + listboxPos(3) + xOffset, listboxPos(2) + 0.5*listboxPos(4) - 0.5*promptHeight, promptHeight + 0.01, promptHeight - 0.01];
icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', addButtonPos, 'string', '+', 'tag', 'add_comp_button', 'fontsize',...
    UI_FS - 1, 'callback', {@addCompNetwork, InputHandle});
icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', removeButtonPos, 'string', '-', 'tag', 'remove_comp_button', 'fontsize',...
    UI_FS - 1, 'callback', {@removeCompNetwork, InputHandle});


promptPos(2) = listboxYOrigin - yOffset - 0.5*promptHeight;
promptPos(3) = promptWidth;
textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Enter P-value Significance Threshold', 'tag', 'prompt_p_threshold', 'fontsize', UI_FS - 1);
icatb_wrapStaticText(textH);

editPos = promptPos;
editPos(1) = editPos(1) + editPos(3) + xOffset;
editPos(3) = controlWidth;
editH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'edit', 'position', editPos, 'string', num2str(mancovanInfo.userInput.p_threshold), 'tag', 'p_threshold', 'fontsize', UI_FS - 1);


promptPos(2) = promptPos(2) - yOffset - 0.5*promptHeight;
promptPos(3) = promptWidth;
textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Enter TR in seconds', 'tag', 'prompt_p_tr', 'fontsize', UI_FS - 1);
icatb_wrapStaticText(textH);

editPos = promptPos;
editPos(1) = editPos(1) + editPos(3) + xOffset;
editPos(3) = controlWidth;
editH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'edit', 'position', editPos, 'string', num2str(mancovanInfo.userInput.TR), 'tag', 'tr', 'fontsize', UI_FS - 1);


%% Add an option to do MDL estimation or give the user the option to enter
% PCS for each feature.
promptPos(2) = promptPos(2) - yOffset - 0.5*promptHeight;
promptPos(3) = 1 - 2*xOffset;
icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'checkbox', 'position', promptPos, 'string', 'Autoselect No. Of Components For Each Feature Using MDL', 'tag', ...
    'checkbox_mdl', 'fontsize', UI_FS - 1, 'value', mancovanInfo.userInput.doEstimation, 'callback', {@autoMDLCallback, InputHandle});

compEnable = 'on';
compVisible = 'on';
if (mancovanInfo.userInput.doEstimation == 1)
    compEnable = 'inactive';
    compVisible = 'off';
end

designCriteria = 'mancova';

try
    designCriteria = mancovanInfo.userInput.designCriteria;
catch
end

if (~strcmpi(designCriteria, 'mancova'))
    compEnable = 'inactive';
    compVisible = 'off';
end

promptPos(2) = promptPos(2) - yOffset - 0.5*promptHeight;
promptPos(3) = promptWidth;
textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Enter No. Of Components For Each Feature In a Vector', 'tag', ...
    'prompt_comps', 'fontsize', UI_FS - 1, 'enable', compEnable, 'visible', compVisible);
icatb_wrapStaticText(textH);
editPos = promptPos;
promptPos = get(textH, 'position');
promptPos(2)= promptPos(2) + (editPos(4) - promptPos(4));
set(textH, 'position', promptPos);

editPos(2) = promptPos(2) + 0.5*promptPos(4) - 0.5*editPos(4);
editPos(1) = promptPos(1) + promptPos(3) + xOffset;
editPos(3) = controlWidth;
editH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'edit', 'position', editPos, 'string', num2str(mancovanInfo.userInput.numOfPCs), 'tag', 'answer_comps', 'fontsize', UI_FS - 1, ...
    'enable', compEnable, 'visible', compVisible);

%% Add cancel, save and run buttons
cancelPos = [0.25 - 0.5*okWidth, promptPos(2) - yOffset - 0.5*okHeight, okWidth, okHeight];
cancelPos(2) = cancelPos(2) - 0.5*cancelPos(4);
icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', cancelPos, 'string', 'Cancel', 'tag', 'cancel_button', 'fontsize',...
    UI_FS - 1, 'callback', 'delete(gcbf);');

okPos = [0.75 - 0.5*okWidth, promptPos(2) - yOffset - 0.5*okHeight, okWidth, okHeight];
okPos(2) = okPos(2) - 0.5*okPos(4);
icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', okPos, 'string', 'Run', 'tag', 'run_button', 'fontsize',...
    UI_FS - 1, 'callback', {@runCallback, InputHandle});

savePos = [0.5 - 0.5*okWidth, promptPos(2) - yOffset - 0.5*okHeight, okWidth, okHeight];
savePos(2) = savePos(2) - 0.5*savePos(4);
icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', savePos, 'string', 'Save', 'tag', 'save_button', 'fontsize',...
    UI_FS - 1, 'callback', {@saveCallback, InputHandle});

try
    delete(msgH);
catch
end

set(InputHandle, 'visible', 'on');
drawnow;


function addCompNetwork(hObject, event_data, figH)
%% Add Component network
%

icatb_defaults;
global UI_FS;

figureTag = 'add_comp_mancovan';
compFigHandle = findobj('tag', figureTag);
if (~isempty(compFigHandle))
    delete(compFigHandle);
end

mancovanInfo = get(figH, 'userdata');

listH = findobj(figH, 'tag', 'comp');

compVals = [];
networkName = '';
if (listH == hObject)
    if (~strcmpi(get(figH, 'selectionType'), 'open'))
        return;
    end
    val = get(listH, 'value');
    try
        networkName = mancovanInfo.userInput.comp(val).name;
        compVals = mancovanInfo.userInput.comp(val).value;
    catch
    end
end

compStr = num2str((1:mancovanInfo.userInput.numICs)');

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
compListH = icatb_uicontrol('parent', compFigHandle, 'units', 'normalized', 'style', 'listbox', 'position', listboxPos, 'string', compStr, 'tag', 'components', 'fontsize', UI_FS - 1, 'min', 0, 'max', 2, 'value', compVals);

%% Show components
showWidth = 0.08; showHeight = 0.04;
showButtonPos = [listboxXOrigin + 0.5*listbox2Wdith - 0.5*showWidth, listboxYOrigin - yOffset - 0.5*showHeight, showWidth, showHeight];
showH = icatb_uicontrol('parent', compFigHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', showButtonPos, 'string', 'Show', 'fontsize', UI_FS - 1, 'callback', {@drawComp, figH, compFigHandle});

%% Plot image on the left hand side
axesPos = [xOffset, listboxYOrigin, listboxHeight, listboxHeight];
axesH = axes('parent', compFigHandle, 'units', 'normalized', 'position', axesPos, 'tag', 'axes_display_comp');

promptPos = axesPos;

%% Add cancel and run buttons
okPos = [0.5 - 0.5*okWidth, promptPos(2) - yOffset - 0.5*okHeight, okWidth, okHeight];
okPos(2) = okPos(2) - 0.5*okPos(4);
icatb_uicontrol('parent', compFigHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', okPos, 'string', 'Done', 'tag', 'done_button', 'fontsize', UI_FS - 1, 'callback', {@setCompCallback, ...
    compFigHandle, figH});

%% Draw components on the left hand side
drawComp(compListH, [], figH, compFigHandle);

function setCompCallback(hObject, event_data, compFigH, handles)
%% Get fields from component

mancovanInfo = get(handles, 'userdata');
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
    
    if (length(mancovanInfo.userInput.comp) > 0)
        chk = strmatch(lower(networkName), lower(cellstr(char(mancovanInfo.userInput.comp.name))), 'exact');
        if (~isempty(chk))
            ind = chk;
        end
    end
    
    if (~exist('ind', 'var'))
        ind = length(mancovanInfo.userInput.comp) + 1;
    end
    
    %% Set user selected information in figure
    mancovanInfo.userInput.comp(ind).name = networkName;
    mancovanInfo.userInput.comp(ind).value =  comps(:)';
    set(handles, 'userdata', mancovanInfo);
    compListH = findobj(handles, 'tag', 'comp');
    set(compListH, 'string', cellstr(char(mancovanInfo.userInput.comp.name)));
    delete(compFigH);
    
catch
    icatb_errorDialog(lasterr, 'Component Selection');
end


function removeCompNetwork(hObject, event_data, figH)
%% Remove Component network
%

mancovanInfo = get(figH, 'userdata');
listH = findobj(figH, 'tag', 'comp');
val = get(listH, 'value');

if (~isempty(val))
    check = icatb_questionDialog('title', 'Remove Component Network', 'textbody', 'Do you want to remove the component network from the list?');
    if (~check)
        return;
    end
end

try
    strs = cellstr(char(mancovanInfo.userInput.comp.name));
    mancovanInfo.userInput.comp(val) = [];
    strs(val) = [];
    set(listH, 'value', 1);
    set(listH, 'string', strs);
    set(figH, 'userdata', mancovanInfo);
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

function saveCallback(hObject, event_data, handles)
%% Save the mancovanInfo

mancovanInfo = get_info_controls(handles);
fileN = fullfile(mancovanInfo.userInput.outputDir, [mancovanInfo.userInput.prefix, '.mat']);
icatb_save(fileN, 'mancovanInfo');
wH = icatb_dialogBox('title', 'File Saved', 'textBody', ['Mancovan information is saved in file ', fileN], 'textType', 'large');
waitfor(wH);
delete(handles);

function runCallback(hObject, event_data, handles)
%% Run mancovan

mancovanInfo = get_info_controls(handles);

if (~isfield(mancovanInfo.userInput, 'feature_params'))
    defaultsMenuH = findobj(handles, 'tag', 'mancovan_defaults');
    defaultsCallback(defaultsMenuH, [], handles, 'off');
end

mancovanInfo = get(handles, 'userdata');

try
    
    if (~isfield(mancovanInfo.userInput, 'features') || isempty(mancovanInfo.userInput.features))
        error('Please select the features of interest');
    end
    
    
    if (~isfield(mancovanInfo.userInput, 'comp') || isempty(mancovanInfo.userInput.comp))
        error('Please select component network names');
    end
    
    comps = [mancovanInfo.userInput.comp.value];
    
    if (length(comps) ~= length(unique(comps)))
        error('Please check if there are duplicate entries of the same component in different networks');
    end
    
    if (~mancovanInfo.userInput.doEstimation)
        if (~isfield(mancovanInfo.userInput, 'numOfPCs'))
            error('Enter no. of PCs for each feature');
        end
        
        if (~isnumeric(mancovanInfo.userInput.numOfPCs))
            mancovanInfo.userInput.numOfPCs = str2num(mancovanInfo.userInput.numOfPCs);
        end
        
        if (length(mancovanInfo.userInput.numOfPCs) ~= length(mancovanInfo.userInput.features))
            error(['No. Of PCs must be a vector of length equal to the no. of features selected (', num2str(length(mancovanInfo.userInput.features)), ')']);
        end
        
        if (any(mancovanInfo.userInput.numOfPCs > mancovanInfo.userInput.numOfSub))
            error(['One/more PCs exceed the no. of subjects (', num2str(mancovanInfo.userInput.numOfSub), ')']);
        end
        
    end
    
    %% Run callback
    fileN = fullfile(mancovanInfo.userInput.outputDir, [mancovanInfo.userInput.prefix, '.mat']);
    icatb_save(fileN, 'mancovanInfo');
    delete(handles);
    drawnow;
    
    % Compute features
    icatb_run_mancovan(mancovanInfo, 2);
    
catch
    icatb_errorDialog(lasterr, 'Run Mancovan Error');
    rethrow(lasterror);
end

drawnow;



function mancovanInfo = get_info_controls(handles)
%% Get information from user controls

mancovanInfo = get(handles, 'userdata');
featuresH = findobj(handles, 'tag', 'features');
featuresVal = get(featuresH, 'value');
opts = cellstr(lower(get(featuresH, 'string')));
mancovanInfo.userInput.features = opts(featuresVal);
mancovanInfo.userInput.doEstimation = get(findobj(handles, 'tag', 'checkbox_mdl'), 'value');
try
    mancovanInfo.userInput.numOfPCs = str2num(get(findobj(handles, 'tag', 'answer_comps'), 'string'));
catch
end

pH = findobj(handles, 'tag', 'p_threshold');
TRH = findobj(handles, 'tag', 'tr');
p_threshold = str2num(get(pH, 'string'));
TR = str2num(get(TRH, 'string'));
mancovanInfo.userInput.p_threshold = p_threshold;
mancovanInfo.userInput.TR = TR;

set(handles, 'userdata', mancovanInfo);

function autoMDLCallback(hObject, event_data, handles)

compEnable = 'on';
compVisible = 'on';

if (get(hObject, 'value') == 1)
    compEnable = 'inactive';
    compVisible = 'off';
end

promptH = findobj(handles, 'tag', 'prompt_comps');
editH = findobj(handles, 'tag', 'answer_comps');

set(promptH, 'enable', compEnable);
set(editH, 'enable', compEnable);
set(promptH, 'visible', compVisible);
set(editH, 'visible', compVisible);

if (strcmpi(compVisible, 'on'))
    c = get(editH, 'backgroundcolor');
    set(editH, 'backgroundcolor', [0, 0, 0]);
    set(editH, 'backgroundcolor', c);
end


function unzipFiles(sesInfo, compFiles, inDir)
%% Unzip files

zipContents.zipFiles = {};
zipContents.files_in_zip(1).name = {};
if isfield(sesInfo, 'zipContents')
    zipContents = sesInfo.zipContents;
end

zipFile = icatb_getViewingSet_zip(compFiles, [], 'real', zipContents);
if (~isempty(zipFile))
    icatb_unzip(fullfile(inDir, zipFile), inDir);
end


function defaultsCallback(hObject, ed, handles, figVisibility)
%% Defaults callback

if (~exist('figVisibility', 'var'))
    figVisibility = 'on';
end


mancovanInfo = get_info_controls(handles);

try
    covInfo = mancovanInfo.userInput.feature_params.final.tc_covariates;
catch
end

if (~exist('covInfo', 'var'))
    covInfo.numOfDataSets = mancovanInfo.userInput.numOfSub*mancovanInfo.userInput.numOfSess;
end

mancovanInfo = get_info_controls(handles);
default_params = icatb_mancovan_feature_options('tr', mancovanInfo.userInput.TR, 'mask_dims', mancovanInfo.userInput.HInfo(1).dim(1:3), 'covInfo', covInfo);

if (~isfield(mancovanInfo.userInput, 'feature_params'))
    feature_params = default_params;
else
    feature_params = mancovanInfo.userInput.feature_params;
end


tags = cellstr(char(feature_params.inputParameters(3).options.tag));
motion_inds = strmatch('fnc_tc_covariates', tags, 'exact');
if (isempty(motion_inds))
    feature_params.inputParameters(3).options(end+1) = feature_params.inputParameters(3).options(end);
    feature_params.inputParameters(3).options(end).callback = [];
    feature_params.inputParameters(3).options(end) = default_params.inputParameters(3).options(end);
end

for nF = 1:length(feature_params.inputParameters)
    for nO = 1:length(feature_params.inputParameters(nF).options)
        if (~isfield(feature_params.inputParameters(nF).options(nO), 'enable'))
            feature_params.inputParameters(nF).options(nO).enable = 'on';
            feature_params.defaults(nF).options(nO).enable = 'on';
        end
    end
end

%% Include stats for selecting default mask
if (length(feature_params.inputParameters(1).options) == 2)
    feature_params.inputParameters(1).options(3) = default_params.inputParameters(1).options(3);
    feature_params.inputParameters(1).options(4) = default_params.inputParameters(1).options(4);
    feature_params.defaults(1).options(3) = default_params.defaults(1).options(3);
    feature_params.defaults(1).options(4) = default_params.defaults(1).options(4);
end

%% Set stats for selecting default mask to off for user specified mask
if (feature_params.inputParameters(1).options(1).value == 2)
    feature_params.inputParameters(1).options(3).enable = 'off';
    feature_params.defaults(1).options(3).enable = 'on';
end

tags = cellstr(char(feature_params.inputParameters(3).options.tag));

chk = strmatch('fnc_tc_filter', tags, 'exact');
if (strcmpi(feature_params.inputParameters(3).options(chk).uiType, 'popup') || strcmpi(feature_params.inputParameters(3).options(chk).uiType, 'popupmenu'))
    
    feature_params.inputParameters(3).options(chk).uiType = 'edit';
    feature_params.inputParameters(3).options(chk).answerType = 'numeric';
    
    tmp_opts = cellstr(feature_params.inputParameters(3).options(chk).answerString);
    feature_params.inputParameters(3).options(chk).value = 1;
    
    if (strcmpi(tmp_opts{feature_params.inputParameters(3).options(chk).value}, 'yes'))
        feature_params.inputParameters(3).options(chk).answerString = '0.15';
    else
        feature_params.inputParameters(3).options(chk).answerString = '0';
    end
    
end

out = icatb_OptionsWindow(feature_params.inputParameters, feature_params.defaults, figVisibility, 'title', 'Feature Options');

drawnow;

feature_params.inputParameters = out.inputParameters;
feature_params.final = out.results;
mancovanInfo.userInput.feature_params = feature_params;
set(handles, 'userdata', mancovanInfo);


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