function icatb_setup_dfnc(param_file)
%% Setup dFNC
%
% Inputs:
% 1. param_file - ICA Parameter file or dFNC parameter file
%

icatb_defaults;
global UI_FS;
global PARAMETER_INFO_MAT_FILE;

%% Select dFNC file
if (~exist('param_file', 'var'))
    param_file = icatb_selectEntry('title', 'Select ICA/dFNC Parameter File', 'typeEntity', 'file', 'typeSelection', 'single', 'filter', ['*', PARAMETER_INFO_MAT_FILE, '.mat;*dfnc.mat']);
    drawnow;
    if (isempty(param_file))
        error('ICA/dFNC parameter file is not selected');
    end
end

drawnow;

[inDir, paramF, extn] = fileparts(param_file);
if (isempty(inDir))
    inDir = pwd;
end

param_file = fullfile(inDir, [paramF, extn]);
load(param_file);


if (exist('sesInfo', 'var'))
    sesInfo.userInput.pwd = inDir;
    sesInfo.outputDir = inDir;
    outputDir = icatb_selectEntry('typeEntity', 'directory', 'title', 'Select output directory to place dFNC results ...');
    drawnow;
    if (isempty(outputDir))
        error('Output directory is not selected');
    end
    dfncInfo.userInput.ica_param_file = param_file;
    dfncInfo.userInput.outputDir = outputDir;
    dfncInfo.userInput.numOfSub = sesInfo.numOfSub;
    dfncInfo.userInput.numOfSess = sesInfo.numOfSess;
    dfncInfo.userInput.prefix = [sesInfo.userInput.prefix, '_dfnc'];
    compFiles = sesInfo.icaOutputFiles(1).ses(1).name;
    dfncInfo.userInput.compFiles = icatb_rename_4d_file(icatb_fullFile('directory', inDir, 'files', compFiles));
    dfncInfo.userInput.numICs = sesInfo.numComp;
    dfncInfo.userInput.HInfo = sesInfo.HInfo.V(1);
elseif (exist('dfncInfo', 'var'))
    load (dfncInfo.userInput.ica_param_file);
    dfncInfo.userInput.numOfSess = sesInfo.numOfSess;
    clear sesInfo;
    outputDir = inDir;
    dfncInfo.userInput.outputDir = outputDir;
    inDir = fileparts(dfncInfo.userInput.ica_param_file);
else
    error('Selected file is neither ICA parameter file nor dFNC parameter file');
end



drawnow;

if (exist(outputDir, 'dir') ~= 7)
    mkdir(outputDir);
end

cd(outputDir);

dfncInfo.userInput.outputDir = outputDir;
load(dfncInfo.userInput.ica_param_file);

if (~exist(deblank(dfncInfo.userInput.compFiles(1, :)), 'file'))
    unzipFiles(sesInfo, sesInfo.icaOutputFiles(1).ses(1).name, fileparts(dfncInfo.userInput.ica_param_file));
end

if (~exist(icatb_parseExtn(deblank(dfncInfo.userInput.compFiles(1, :))), 'file'))
    subjectICAFiles = icatb_parseOutputFiles('icaOutputFiles', sesInfo.icaOutputFiles, 'numOfSub', sesInfo.numOfSub, 'numOfSess', sesInfo.numOfSess);
    compFiles = subjectICAFiles(1).ses(1).name;
    dfncInfo.userInput.compFiles = icatb_rename_4d_file(icatb_fullFile('directory', inDir, 'files', compFiles));
    if (~exist(deblank(dfncInfo.userInput.compFiles(1, :)), 'file'))
        unzipFiles(sesInfo, compFiles, fileparts(dfncInfo.userInput.ica_param_file));
    end
end


msg = 'Opening Setup dFNC GUI ...';

disp(msg);

msgH = helpdlg(msg, 'Setup dFNC');

drawnow;

compNetworkNameData = getCompData(dfncInfo.userInput.compFiles);

compGroupNames = '';
covNames = '';

if (~isfield(dfncInfo.userInput, 'comp'))
    dfncInfo.userInput.comp = [];
end

try
    compGroupNames = (cellstr(char(dfncInfo.userInput.comp.name)));
catch
end


if (~isfield(dfncInfo.userInput, 'TR'))
    if (isfield(sesInfo, 'TR'))
        dfncInfo.userInput.TR = sesInfo.TR;
    else
        dfncInfo.userInput.TR = 1;
    end
end

clear sesInfo;

%% Draw graphics
figureTag = 'setup_dfnc_gui';
figHandle = findobj('tag', figureTag);
if (~isempty(figHandle))
    delete(figHandle);
end

% Setup figure for GUI
InputHandle = icatb_getGraphics('dFNC Setup Analysis', 'normal', figureTag, 'off');
set(InputHandle, 'menubar', 'none');
set(InputHandle, 'userdata', dfncInfo);

defaultsmenu = uimenu('parent', InputHandle, 'label', 'dFNC-Defaults', 'tag', 'dfnc_defaults');
set(defaultsmenu, 'callback', {@defaultsCallback, InputHandle});

promptWidth = 0.52; controlWidth = 0.32; promptHeight = 0.05;
xOffset = 0.02; yOffset = promptHeight; listboxHeight = 0.18; listboxWidth = controlWidth; yPos = 0.9;
okWidth = 0.12; okHeight = promptHeight;

promptPos = [xOffset, yPos - 0.5*yOffset, promptWidth, promptHeight];

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
    'comp', 'fontsize', UI_FS - 1, 'min', 0, 'max', 1, 'value', 1, 'callback', {@addCompNetwork, InputHandle}, 'userdata', compNetworkNameData);

addButtonPos = [listboxPos(1) + listboxPos(3) + xOffset, listboxPos(2) + 0.5*listboxPos(4) + 0.5*promptHeight, promptHeight + 0.01, promptHeight - 0.01];
removeButtonPos = [listboxPos(1) + listboxPos(3) + xOffset, listboxPos(2) + 0.5*listboxPos(4) - 0.5*promptHeight, promptHeight + 0.01, promptHeight - 0.01];
icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', addButtonPos, 'string', '+', 'tag', 'add_comp_button', 'fontsize',...
    UI_FS - 1, 'callback', {@addCompNetwork, InputHandle});
icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', removeButtonPos, 'string', '-', 'tag', 'remove_comp_button', 'fontsize',...
    UI_FS - 1, 'callback', {@removeCompNetwork, InputHandle});


promptPos(2) = listboxYOrigin - yOffset - 0.5*promptHeight;
promptPos(3) = promptWidth;
textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Enter TR in seconds', 'tag', 'prompt_p_tr', 'fontsize', UI_FS - 1);
icatb_wrapStaticText(textH);

editPos = promptPos;
editPos(1) = editPos(1) + editPos(3) + xOffset;
editPos(3) = controlWidth;
editH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'edit', 'position', editPos, 'string', num2str(dfncInfo.userInput.TR), 'tag', 'tr', 'fontsize', UI_FS - 1);


%% Add an option to do MDL estimation or give the user the option to enter
% PCS for each feature.
% promptPos(2) = promptPos(2) - yOffset - 0.5*promptHeight;
% promptPos(3) = 1 - 2*xOffset;
% icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'checkbox', 'position', promptPos, 'string', 'Autoselect No. Of Components For Each Feature Using MDL', 'tag', ...
%     'checkbox_mdl', 'fontsize', UI_FS - 1, 'value', dfncInfo.userInput.doEstimation, 'callback', {@autoMDLCallback, InputHandle});

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

figureTag = 'add_comp_dfnc';
compFigHandle = findobj('tag', figureTag);
if (~isempty(compFigHandle))
    delete(compFigHandle);
end

dfncInfo = get(figH, 'userdata');

listH = findobj(figH, 'tag', 'comp');

compVals = [];
networkName = '';
if (listH == hObject)
    if (~strcmpi(get(figH, 'selectionType'), 'open'))
        return;
    end
    val = get(listH, 'value');
    try
        networkName = dfncInfo.userInput.comp(val).name;
        compVals = dfncInfo.userInput.comp(val).value;
    catch
    end
end

compStr = num2str((1:dfncInfo.userInput.numICs)');

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

dfncInfo = get(handles, 'userdata');
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
    
    if (length(dfncInfo.userInput.comp) > 0)
        chk = strmatch(lower(networkName), lower(cellstr(char(dfncInfo.userInput.comp.name))), 'exact');
        if (~isempty(chk))
            ind = chk;
        end
    end
    
    if (~exist('ind', 'var'))
        ind = length(dfncInfo.userInput.comp) + 1;
    end
    
    %% Set user selected information in figure
    dfncInfo.userInput.comp(ind).name = networkName;
    dfncInfo.userInput.comp(ind).value =  comps(:)';
    set(handles, 'userdata', dfncInfo);
    compListH = findobj(handles, 'tag', 'comp');
    set(compListH, 'string', cellstr(char(dfncInfo.userInput.comp.name)));
    delete(compFigH);
    
catch
    icatb_errorDialog(lasterr, 'Component Selection');
end


function removeCompNetwork(hObject, event_data, figH)
%% Remove Component network
%

dfncInfo = get(figH, 'userdata');
listH = findobj(figH, 'tag', 'comp');
val = get(listH, 'value');

if (~isempty(val))
    check = icatb_questionDialog('title', 'Remove Component Network', 'textbody', 'Do you want to remove the component network from the list?');
    if (~check)
        return;
    end
end

try
    strs = cellstr(char(dfncInfo.userInput.comp.name));
    dfncInfo.userInput.comp(val) = [];
    strs(val) = [];
    set(listH, 'value', 1);
    set(listH, 'string', strs);
    set(figH, 'userdata', dfncInfo);
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
%% Save the dfncInfo

dfncInfo = get_info_controls(handles);
fileN = fullfile(dfncInfo.userInput.outputDir, [dfncInfo.userInput.prefix, '.mat']);
icatb_save(fileN, 'dfncInfo');
wH = icatb_dialogBox('title', 'File Saved', 'textBody', ['dFNC information is saved in file ', fileN], 'textType', 'large');
waitfor(wH);
delete(handles);

function runCallback(hObject, event_data, handles)
%% Run mancovan

dfncInfo = get_info_controls(handles);

if (~isfield(dfncInfo.userInput, 'feature_params'))
    defaultsMenuH = findobj(handles, 'tag', 'dfnc_defaults');
    defaultsCallback(defaultsMenuH, [], handles, 'off');
end

dfncInfo = get(handles, 'userdata');

try
    
    
    if (~isfield(dfncInfo.userInput, 'comp') || isempty(dfncInfo.userInput.comp))
        error('Please select component network names');
    end
    
    comps = [dfncInfo.userInput.comp.value];
    
    if (length(comps) ~= length(unique(comps)))
        error('Please check if there are duplicate entries of the same component in different networks');
    end
    
    %% Run callback
    fileN = fullfile(dfncInfo.userInput.outputDir, [dfncInfo.userInput.prefix, '.mat']);
    icatb_save(fileN, 'dfncInfo');
    delete(handles);
    drawnow;
    
    icatb_run_dfnc(dfncInfo);
    
catch
    icatb_errorDialog(lasterr, 'Run dFNC Error');
    rethrow(lasterror);
end

drawnow;



function dfncInfo = get_info_controls(handles)
%% Get information from user controls

dfncInfo = get(handles, 'userdata');
TRH = findobj(handles, 'tag', 'tr');
TR = str2num(get(TRH, 'string'));
dfncInfo.userInput.TR = TR;

set(handles, 'userdata', dfncInfo);


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

dfncInfo = get_info_controls(handles);

try
    covInfo = dfncInfo.userInput.feature_params.final.tc_covariates;
catch
end

if (~exist('covInfo', 'var'))
    covInfo.numOfDataSets = dfncInfo.userInput.numOfSub*dfncInfo.userInput.numOfSess;
end

default_params = icatb_dfnc_options('covInfo', covInfo);

if (~isfield(dfncInfo.userInput, 'feature_params'))
    feature_params = default_params;
else
    feature_params = dfncInfo.userInput.feature_params;
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

tags = cellstr(char(feature_params(1).inputParameters(2).options.tag));
check_window_type = strmatch('window_type', tags, 'exact');
check_window_alpha = strmatch('window_alpha', tags, 'exact');
if (~isempty(check_window_type))
    if (feature_params(1).inputParameters(2).options(check_window_type).value == 1)
        feature_params(1).inputParameters(2).options(check_window_alpha).answerString = '3';
    end
    feature_params(1).inputParameters(2).options(check_window_type) = [];
end

out = icatb_OptionsWindow(feature_params.inputParameters, feature_params.defaults, figVisibility, ...
    'title', 'dFNC Options');

drawnow;

feature_params.inputParameters = out.inputParameters;
feature_params.final = out.results;
dfncInfo.userInput.feature_params = feature_params;
set(handles, 'userdata', dfncInfo);


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
