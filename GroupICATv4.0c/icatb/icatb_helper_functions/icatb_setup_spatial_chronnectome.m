function icatb_setup_spatial_chronnectome(param_file)
%% Setup spatial chronnectome
%
% Inputs:
% 1. param_file - ICA Parameter file or spatial chronnectome parameter file
%
%

icatb_defaults;
global UI_FS;
global PARAMETER_INFO_MAT_FILE;

%% Select spatial chronnectome file
if (~exist('param_file', 'var'))
    param_file = icatb_selectEntry('title', 'Select ICA/spatial chronnectome Parameter File', 'typeEntity', 'file', 'typeSelection', 'single', 'filter', ['*', PARAMETER_INFO_MAT_FILE, '.mat;*schronn.mat']);
    drawnow;
    if (isempty(param_file))
        error('ICA/spatial chronnectome parameter file is not selected');
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
    outputDir = icatb_selectEntry('typeEntity', 'directory', 'title', 'Select output directory to place spatial chronnectome results ...');
    drawnow;
    if (isempty(outputDir))
        error('Output directory is not selected');
    end
    schronnInfo.userInput.ica_param_file = param_file;
    schronnInfo.userInput.outputDir = outputDir;
    schronnInfo.userInput.numOfSub = sesInfo.numOfSub;
    schronnInfo.userInput.numOfSess = sesInfo.numOfSess;
    schronnInfo.userInput.prefix = [sesInfo.userInput.prefix, '_schronn'];
    compFiles = sesInfo.icaOutputFiles(1).ses(1).name;
    schronnInfo.userInput.compFiles = icatb_rename_4d_file(icatb_fullFile('directory', inDir, 'files', compFiles));
    schronnInfo.userInput.numICs = sesInfo.numComp;
    schronnInfo.userInput.HInfo = sesInfo.HInfo.V(1);
elseif (exist('schronnInfo', 'var'))
    load (schronnInfo.userInput.ica_param_file);
    schronnInfo.userInput.numOfSess = sesInfo.numOfSess;
    clear sesInfo;
    outputDir = inDir;
    schronnInfo.userInput.outputDir = outputDir;
    inDir = fileparts(schronnInfo.userInput.ica_param_file);
else
    error('Selected file is neither ICA parameter file nor spatial chronnectome parameter file');
end



drawnow;

if (exist(outputDir, 'dir') ~= 7)
    mkdir(outputDir);
end

cd(outputDir);

schronnInfo.userInput.outputDir = outputDir;
load(schronnInfo.userInput.ica_param_file);

if (~exist(deblank(schronnInfo.userInput.compFiles(1, :)), 'file'))
    unzipFiles(sesInfo, sesInfo.icaOutputFiles(1).ses(1).name, fileparts(schronnInfo.userInput.ica_param_file));
end

if (~exist(icatb_parseExtn(deblank(schronnInfo.userInput.compFiles(1, :))), 'file'))
    subjectICAFiles = icatb_parseOutputFiles('icaOutputFiles', sesInfo.icaOutputFiles, 'numOfSub', sesInfo.numOfSub, 'numOfSess', sesInfo.numOfSess);
    compFiles = subjectICAFiles(1).ses(1).name;
    schronnInfo.userInput.compFiles = icatb_rename_4d_file(icatb_fullFile('directory', inDir, 'files', compFiles));
    if (~exist(deblank(schronnInfo.userInput.compFiles(1, :)), 'file'))
        unzipFiles(sesInfo, compFiles, fileparts(schronnInfo.userInput.ica_param_file));
    end
end


msg = 'Opening Setup Spatial Chronnectome GUI ...';

disp(msg);

msgH = helpdlg(msg, 'Setup Spatial Chronnectome');

drawnow;


covNames = '';

if (~isfield(schronnInfo.userInput, 'comp'))
    schronnInfo.userInput.comp = [];
end



if (~isfield(schronnInfo.userInput, 'TR'))
    if (isfield(sesInfo, 'TR'))
        schronnInfo.userInput.TR = sesInfo.TR;
    else
        schronnInfo.userInput.TR = 1;
    end
end


compStr = cellstr(num2str((1:sesInfo.numComp)'));

clear sesInfo;

%% Draw graphics
figureTag = 'setup_schronn_gui';
figHandle = findobj('tag', figureTag);
if (~isempty(figHandle))
    delete(figHandle);
end

% Setup figure for GUI
InputHandle = icatb_getGraphics('Spatial Chronnectome Setup Analysis', 'normal', figureTag, 'off');
set(InputHandle, 'menubar', 'none');
set(InputHandle, 'userdata', schronnInfo);

defaultsmenu = uimenu('parent', InputHandle, 'label', 'SpatialChronnectome-Defaults', 'tag', 'spatial_chronn_defaults');
set(defaultsmenu, 'callback', {@defaultsCallback, InputHandle});

promptWidth = 0.52; controlWidth = 0.32; promptHeight = 0.05;
xOffset = 0.02; yOffset = promptHeight; listboxHeight = 0.18; listboxWidth = controlWidth; yPos = 0.9;
okWidth = 0.12; okHeight = promptHeight;

promptPos = [xOffset, yPos - 0.5*yOffset, promptWidth, promptHeight];

listboxXOrigin = promptPos(1) + promptPos(3) + xOffset;
listboxYOrigin = promptPos(2) - 0.5*listboxHeight;

%%  Components listbox (Group components by name)
%promptPos(2) = listboxYOrigin - yOffset - 0.5*listboxHeight;
textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Select Components', 'tag', ...
    'prompt_components', 'fontsize', UI_FS - 1);
icatb_wrapStaticText(textH);
listboxYOrigin = promptPos(2) - 0.5*listboxHeight;
listboxPos = [listboxXOrigin, listboxYOrigin, listboxWidth, listboxHeight];
listCompH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'listbox', 'position', listboxPos, 'string', compStr, 'tag', ...
    'comp', 'fontsize', UI_FS - 1, 'min', 0, 'max', 2, 'value', schronnInfo.userInput.comp);

promptPos(2) = listboxYOrigin - yOffset - 0.5*promptHeight;
promptPos(3) = promptWidth;
textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Enter TR in seconds', 'tag', 'prompt_p_tr', 'fontsize', UI_FS - 1);
icatb_wrapStaticText(textH);

editPos = promptPos;
editPos(1) = editPos(1) + editPos(3) + xOffset;
editPos(3) = controlWidth;
editH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'edit', 'position', editPos, 'string', num2str(schronnInfo.userInput.TR), 'tag', 'tr', 'fontsize', UI_FS - 1);


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




function saveCallback(hObject, event_data, handles)
%% Save the schronnInfo

schronnInfo = get_info_controls(handles);
fileN = fullfile(schronnInfo.userInput.outputDir, [schronnInfo.userInput.prefix, '.mat']);
icatb_save(fileN, 'schronnInfo');
wH = icatb_dialogBox('title', 'File Saved', 'textBody', ['Spatial Chronnectome information is saved in file ', fileN], 'textType', 'large');
waitfor(wH);
delete(handles);

function runCallback(hObject, event_data, handles)
%% Run mancovan

schronnInfo = get_info_controls(handles);

if (~isfield(schronnInfo.userInput, 'feature_params'))
    defaultsMenuH = findobj(handles, 'tag', 'spatial_chronn_defaults');
    defaultsCallback(defaultsMenuH, [], handles, 'off');
end

schronnInfo = get(handles, 'userdata');

listCompH = findobj(handles, 'tag', 'comp');
compvals = get(listCompH, 'value');
compStr = get(listCompH, 'string');
schronnInfo.userInput.comp = (str2num(char(compStr(compvals))));


try
    
    
    if (~isfield(schronnInfo.userInput, 'comp') || isempty(schronnInfo.userInput.comp))
        error('Please select components');
    end
    
    
    %% Run callback
    fileN = fullfile(schronnInfo.userInput.outputDir, [schronnInfo.userInput.prefix, '.mat']);
    icatb_save(fileN, 'schronnInfo');
    delete(handles);
    drawnow;
    
    icatb_run_spatial_chronnectome(schronnInfo);
    
catch
    icatb_errorDialog(lasterr, 'Run dFNC Error');
    rethrow(lasterror);
end

drawnow;



function schronnInfo = get_info_controls(handles)
%% Get information from user controls

schronnInfo = get(handles, 'userdata');
TRH = findobj(handles, 'tag', 'tr');
TR = str2num(get(TRH, 'string'));
schronnInfo.userInput.TR = TR;

set(handles, 'userdata', schronnInfo);


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

schronnInfo = get_info_controls(handles);

try
    covInfo = schronnInfo.userInput.feature_params.final.tc_covariates;
catch
end

if (~exist('covInfo', 'var'))
    covInfo.numOfDataSets = schronnInfo.userInput.numOfSub*schronnInfo.userInput.numOfSess;
end

default_params = icatb_spatial_chronnectome_options('covInfo', covInfo);

if (~isfield(schronnInfo.userInput, 'feature_params'))
    feature_params = default_params;
else
    feature_params = schronnInfo.userInput.feature_params;
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
    'title', 'Spatial Chronnectome Defaults');

drawnow;

feature_params.inputParameters = out.inputParameters;
feature_params.final = out.results;
schronnInfo.userInput.feature_params = feature_params;
set(handles, 'userdata', schronnInfo);
