function icatb_setup_spatial_dfnc(param_file)
%% Setup spatial dfnc
%

icatb_defaults;
global UI_FS;
global PARAMETER_INFO_MAT_FILE;


%% Select sdFNC file
if (~exist('param_file', 'var'))
    param_file = icatb_selectEntry('title', 'Select ICA/spatial dFNC Parameter File', 'typeEntity', 'file', 'typeSelection', 'single', 'filter', ['*', PARAMETER_INFO_MAT_FILE, '.mat;*_sdfnc.mat']);
    drawnow;
    if (isempty(param_file))
        error('ICA/sdFNC parameter file is not selected');
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
    outputDir = icatb_selectEntry('typeEntity', 'directory', 'title', 'Select output directory to place spatial dFNC results ...');
    drawnow;
    if (isempty(outputDir))
        error('Output directory is not selected');
    end
    
    if (any(sesInfo.userInput.diffTimePoints ~= sesInfo.userInput.diffTimePoints(1)))
        error('No. of timepoints must be the same between subjects');
    end
    
    sdfncInfo.userInput.files = sesInfo.userInput.files;
    sdfncInfo.userInput.outputDir = outputDir;
    sdfncInfo.userInput.numOfSub = sesInfo.userInput.numOfSub;
    sdfncInfo.userInput.numOfSess = sesInfo.userInput.numOfSess;
    sdfncInfo.userInput.prefix = sesInfo.userInput.prefix;
    sdfncInfo.userInput.HInfo = sesInfo.userInput.HInfo.V(1);
    sdfncInfo.userInput.mask_ind = sesInfo.userInput.mask_ind;
    sdfncInfo.userInput.timepoints = sesInfo.userInput.diffTimePoints(1);
elseif (exist('sdfncInfo', 'var'))
    outputDir = inDir;
    sdfncInfo.userInput.outputDir = outputDir;
else
    error('Selected file is neither ICA parameter file nor spatial dFNC parameter file');
end


drawnow;

if (exist(outputDir, 'dir') ~= 7)
    mkdir(outputDir);
end

cd(outputDir);

outputDir = pwd;
prefix = '';
numComp = 8;
groupStr = '';
window_size = ceil(sdfncInfo.userInput.timepoints/4);

try
    prefix = sdfncInfo.userInput.prefix;
    groupStr =  cellstr(char(sdfncInfo.userInput.group.name));
    numComp = sdfncInfo.userInput.numComp;
    window_size = sdfncInfo.userInput.window_size;
catch
end


figureTag = 'spatial_dfnc_gui';
sdfncInfo.userInput.outputDir = outputDir;

figHandle = findobj('tag', figureTag);
if (~isempty(figHandle))
    delete(figHandle);
end


InputHandle = icatb_getGraphics('Spatial dFNC Setup Analysis', 'normal', figureTag, 'off');
set(InputHandle, 'menubar', 'none');
set(InputHandle, 'userdata', sdfncInfo);

yPos = 0.95; yOffset = 0.05; xOffset  = 0.02;
promptHeight = 0.05;promptWidth = 0.52;
okWidth = 0.12; okHeight = 0.05;
listboxHeight = 0.32;
listboxWidth = 0.32;

promptPos(1) = xOffset;
promptPos(2) = yPos - yOffset;
promptPos(3) = promptWidth;
promptPos(4) = promptHeight;
textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Enter output prefix', 'tag', 'prompt_prefix', ...
    'fontsize', UI_FS - 1);
icatb_wrapStaticText(textH);

editPos = promptPos;
editPos(1) = editPos(1) + editPos(3) + xOffset;
editPos(3) = 0.2;
editH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'edit', 'position', editPos, 'string', prefix, 'tag', 'prefix', 'fontsize', UI_FS - 1);

listboxYOrigin = promptPos(2) - 0.5*listboxHeight;
listboxXOrigin = promptPos(1) + promptPos(3) + xOffset;

%%  Groups listbox
promptPos(2) = listboxYOrigin - 1.2*yOffset;
%promptPos(3) = 0.32;
textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Add Groups', 'tag', ...
    'prompt_groups', 'fontsize', UI_FS - 1);
icatb_wrapStaticText(textH);
listboxYOrigin = promptPos(2) - 0.5*listboxHeight;
listboxPos = [listboxXOrigin, listboxYOrigin, listboxWidth, listboxHeight];
listCompH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'listbox', 'position', listboxPos, 'string', groupStr, 'tag', ...
    'group', 'fontsize', UI_FS - 1, 'min', 0, 'max', 1, 'value', 1, 'callback', {@addGroups, InputHandle});

addButtonPos = [listboxPos(1) + listboxPos(3) + xOffset, listboxPos(2) + 0.5*listboxPos(4) + 0.5*promptHeight, promptHeight + 0.01, promptHeight - 0.01];
removeButtonPos = [listboxPos(1) + listboxPos(3) + xOffset, listboxPos(2) + 0.5*listboxPos(4) - 0.5*promptHeight, promptHeight + 0.01, promptHeight - 0.01];
icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', addButtonPos, 'string', '+', 'tag', 'add_group_button', 'fontsize',...
    UI_FS - 1, 'callback', {@addGroups, InputHandle});
icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', removeButtonPos, 'string', '-', 'tag', 'remove_group_button', 'fontsize',...
    UI_FS - 1, 'callback', {@removeGroups, InputHandle});

promptPos(2) = listboxYOrigin - 1.2*yOffset - 0.5*promptHeight;
promptPos(3) = promptWidth;
promptPos(4) = promptHeight;
textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Enter window size (scans)', 'tag', 'prompt_window', ...
    'fontsize', UI_FS - 1);

editPos = promptPos;
editPos(1) = editPos(1) + editPos(3) + xOffset;
editPos(3) = 0.2;
editH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'edit', 'position', editPos, 'string', num2str(window_size), 'tag', 'window_size', 'fontsize', UI_FS - 1);


promptPos(2) = promptPos(2) - 1.5*yOffset - 0.5*promptHeight;
promptPos(3) = promptWidth;
promptPos(4) = promptHeight;
textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Number of IVA components to extract from the data', 'tag', ...
    'prompt_components', 'fontsize', UI_FS - 1);

editPos = promptPos;
editPos(1) = editPos(1) + editPos(3) + xOffset;
editPos(3) = 0.2;
editH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'edit', 'position', editPos, 'string', num2str(numComp), 'tag', 'numComp', 'fontsize', UI_FS - 1);

icatb_wrapStaticText(textH);
promptPos = get(textH, 'position');


promptPos(2) = promptPos(2) - 1.5*yOffset - 0.5*promptHeight;
promptPos(3) = promptWidth;
promptPos(4) = promptHeight;
textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Number Of Times IVA is run', 'tag', ...
    'prompt_iva_runs', 'fontsize', UI_FS - 1);

editPos = promptPos;
editPos(1) = editPos(1) + editPos(3) + xOffset;
editPos(3) = 0.2;
editH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'edit', 'position', editPos, 'string', num2str(5), 'tag', 'numRuns', 'fontsize', UI_FS - 1);

icatb_wrapStaticText(textH);
promptPos = get(textH, 'position');

okPos = [0.5 - 0.5*okWidth, promptPos(2) - 1.2*yOffset - 0.5*okHeight, okWidth, okHeight];
okPos(2) = okPos(2) - 0.5*okPos(4);
icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', okPos, 'string', 'Run', 'tag', 'run_button', 'fontsize',...
    UI_FS - 1, 'callback', {@runCallback, InputHandle});

if (~isfield(sdfncInfo.userInput, 'group'))
    sdfncInfo.userInput.group = [];
end

set(InputHandle, 'userdata', sdfncInfo);
set(InputHandle, 'visible', 'on');

function addGroups(hObject, event_data, figH)
%% Add groups

sdfncInfo = get(figH, 'userdata');
listH = findobj(figH, 'tag', 'group');
val = get(listH, 'value');

icatb_defaults;
global UI_FS;

figureTag = 'add_groups_sdfnc';
compFigHandle = findobj('tag', figureTag);
if (~isempty(compFigHandle))
    delete(compFigHandle);
end

groupVal = [];
groupName = '';
if (listH == hObject)
    if (~strcmpi(get(figH, 'selectionType'), 'open'))
        return;
    end
    val = get(listH, 'value');
    try
        groupName = sdfncInfo.userInput.group(val).name;
        groupVal = sdfncInfo.userInput.group(val).val;
    catch
    end
end

subjectString = cellstr([repmat('Subject ', sdfncInfo.userInput.numOfSub, 1), num2str((1:sdfncInfo.userInput.numOfSub)')]);

%[groupName, groupVal] = icatb_select_groups_gui(subjectString, groupName, 'select_subjects', groupVal);
[groupName, groupVal] = icatb_select_groups_gui(subjectString, 'Group', 'select_subjects', groupName, groupVal);

try
    
    if (isempty(groupName))
        error('Group name is not selected');
    end
    
    if (isempty(groupVal))
        error('Subjects are not selected');
    end
    
    if (length(sdfncInfo.userInput.group) > 0)
        chk = strmatch(lower(groupName), lower(cellstr(char(sdfncInfo.userInput.group.name))), 'exact');
        if (~isempty(chk))
            ind = chk;
        end
    end
    
    if (~exist('ind', 'var'))
        ind = length(sdfncInfo.userInput.group) + 1;
    end
    
    %% Set user selected information in figure
    sdfncInfo.userInput.group(ind).name = groupName;
    sdfncInfo.userInput.group(ind).val =  groupVal;
    set(figH, 'userdata', sdfncInfo);
    groupListH = findobj(figH, 'tag', 'group');
    set(groupListH, 'string', cellstr(char(sdfncInfo.userInput.group.name)));
    
catch
    icatb_errorDialog(lasterr, 'Group Selection');
end



function removeGroups(hObject, event_data, figH)
%% Remove groups
%

sdfncInfo = get(figH, 'userdata');
listH = findobj(figH, 'tag', 'group');
val = get(listH, 'value');
strs = cellstr(get(listH, 'string'));

if (~isempty(strs))
    check = icatb_questionDialog('title', 'Remove groups', 'textbody', ['Do you want to remove the group ', strs{val}, ' from the list?']);
    if (~check)
        return;
    end
end

try
    strs = cellstr(char(sdfncInfo.userInput.group.name));
    sdfncInfo.userInput.group(val) = [];
    strs(val) = [];
    set(listH, 'value', 1);
    set(listH, 'string', strs);
    set(figH, 'userdata', sdfncInfo);
catch
end


function runCallback(hObject, event_data, handles)
%% Run callback
%

sdfncInfo = get(handles, 'userdata');
prefix = get(findobj(handles, 'tag', 'prefix'), 'string');
window_size = str2num(get(findobj(handles, 'tag', 'window_size'), 'string'));
numComp = str2num(get(findobj(handles, 'tag', 'numComp'), 'string'));
numRuns = str2num(get(findobj(handles, 'tag', 'numRuns'), 'string'));

if (window_size < 0)
    error('Window size is less than zero');
end
if (isempty(sdfncInfo.userInput.group))
    error('Group/Groups are not selected.');
end

sdfncInfo.userInput.prefix = prefix;
sdfncInfo.userInput.window_size = window_size;
sdfncInfo.userInput.numComp = numComp;
sdfncInfo.userInput.num_iva_runs = numRuns;

IVA_Options = icatb_icaOptions([1, 1], 'iva-gl', 'on');
if (isempty(IVA_Options))
    error('IVA options are not selected');
end
sdfncInfo.userInput.IVA_Options = IVA_Options;

drawnow;

outputFile = fullfile(sdfncInfo.userInput.outputDir, [sdfncInfo.userInput.prefix, '_sdfnc.mat']);
icatb_save(outputFile, 'sdfncInfo');
disp(['Spatial dFNC information is saved in file ', outputFile]);
fprintf('\n');
delete(handles);

drawnow;

icatb_run_spatial_dfnc(sdfncInfo);