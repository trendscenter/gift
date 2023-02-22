function icatb_setup_mancovan_design(param_file)
%% Create mancovan design
%
% Inputs:
% 1. param_file - ICA Parameter file or Mancovan parameter file
%

icatb_defaults;
global UI_FS;
global PARAMETER_INFO_MAT_FILE;

if (~exist('param_file', 'var'))
    param_file = icatb_selectEntry('title', 'Select ICA/Mancovan Parameter File', 'typeEntity', 'file', 'typeSelection', 'single', 'filter', ['*', PARAMETER_INFO_MAT_FILE, '.mat;*mancovan.mat']);
    drawnow;
    if (isempty(param_file))
        error('ICA parameter file is not selected');
    end
end


[inDir, paramF, extn] = fileparts(param_file);
if (isempty(inDir))
    inDir = pwd;
end

param_file = fullfile(inDir, [paramF, extn]);


load(param_file);

if (exist('sesInfo', 'var'))
    outputDir = icatb_selectEntry('typeEntity', 'directory', 'title', 'Select output directory to place mancovan results ...');
    drawnow;
    if (isempty(outputDir))
        error('Output directory is not selected');
    end
    mancovanInfo.userInput.ica_param_file = param_file;
    mancovanInfo.userInput.outputDir = outputDir;
    mancovanInfo.userInput.numOfSub = sesInfo.numOfSub;
    mancovanInfo.userInput.numOfSess = sesInfo.numOfSess;
    mancovanInfo.userInput.prefix = [sesInfo.userInput.prefix, '_mancovan'];
    compFiles = sesInfo.icaOutputFiles(1).ses(1).name;
    mancovanInfo.userInput.compFiles = icatb_rename_4d_file(icatb_fullFile('directory', inDir, 'files', compFiles));
    mancovanInfo.userInput.numICs = sesInfo.numComp;
    mancovanInfo.userInput.HInfo = sesInfo.HInfo.V(1);
elseif (exist('mancovanInfo', 'var'))
    outputDir = inDir;
    mancovanInfo.userInput.outputDir = outputDir;
else
    error('Selected file is neither ICA parameter file nor Mancovan parameter file');
end

drawnow;

if (exist(outputDir, 'dir') ~= 7)
    mkdir(outputDir);
end

cd(outputDir);

msg = 'Opening Setup Mancovan Design GUI ...';

disp(msg);

msgH = helpdlg(msg, 'Setup Mancovan Design');

drawnow;

if (mancovanInfo.userInput.numOfSub < 2)
    error('Cannot do stats for one subject');
end

if (~exist('sesInfo', 'var'))
    load(mancovanInfo.userInput.ica_param_file);
    mancovanInfo.userInput.numOfSess = sesInfo.numOfSess;
end

clear sesInfo;

covNames = '';

if (~isfield(mancovanInfo.userInput, 'cov'))
    mancovanInfo.userInput.cov = [];
end

try
    covNames = (cellstr(char(mancovanInfo.userInput.cov.name)));
catch
end

%% Draw graphics
figureTag = 'setup_mancovan_gui';
figHandle = findobj('tag', figureTag);
if (~isempty(figHandle))
    delete(figHandle);
end

% Setup figure for GUI
InputHandle = icatb_getGraphics('Setup Mancovan Design', 'normal', figureTag, 'off');
set(InputHandle, 'menubar', 'none');
set(InputHandle, 'userdata', mancovanInfo);

controlWidth = 0.52;
promptHeight = 0.05;
promptWidth = controlWidth;
listboxHeight = controlWidth; listboxWidth = controlWidth;
xOffset = 0.02; yOffset = promptHeight; yPos = 0.9;
okWidth = 0.12; okHeight = promptHeight;

dropDownWidth = 0.4;

%% Design dropdown box
promptPos = [0.25 - 0.5*dropDownWidth, yPos - 0.5*yOffset, 0.4, promptHeight];

textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Select design criteria ...', 'tag', ...
    'prompt_design', 'fontsize', UI_FS - 1);
designOptions = {'MANCOVA', 'One sample t-test', 'Two sample t-test', 'Paired T-test'};

designValue = 1;
try
    designValue = strmatch(lower(mancovanInfo.userInput.designCriteria), lower(designOptions), 'exact');
catch
end

if (isempty(designValue))
    designValue = 1;
end

icatb_wrapStaticText(textH);
promptPos = get(textH, 'position');

promptPos(1) = promptPos(1) + promptPos(3) + xOffset;
promptPos(3) = 0.3;
icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'popup', 'position', promptPos, 'string', designOptions, 'tag', ...
    'design_criteria', 'fontsize', UI_FS - 1, 'value', designValue, 'callback', {@designCriteriaCallback, InputHandle});

yPos = yPos - promptPos(4) - yOffset;

%% Features text and listbox
promptPos = [0.5 - 0.5*controlWidth, yPos - 0.5*yOffset, promptWidth, promptHeight];

textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Add Covariates For Mancova', 'tag', ...
    'prompt_covariates', 'fontsize', UI_FS - 1);
icatb_wrapStaticText(textH);
listboxXOrigin = promptPos(1);
listboxYOrigin = promptPos(2) - yOffset - listboxHeight;
listboxPos = [listboxXOrigin, listboxYOrigin, listboxWidth, listboxHeight];
icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'listbox', 'position', listboxPos, 'string', covNames, 'tag', ...
    'cov', 'fontsize', UI_FS - 1, 'min', 0, 'max', 1, 'value', 1,  'callback', {@addCov, InputHandle});

addButtonPos = [listboxPos(1) + listboxPos(3) + xOffset, listboxPos(2) + 0.5*listboxPos(4) + 0.5*promptHeight, promptHeight + 0.01, promptHeight - 0.01];
removeButtonPos = [listboxPos(1) + listboxPos(3) + xOffset, listboxPos(2) + 0.5*listboxPos(4) - 0.5*promptHeight, promptHeight + 0.01, promptHeight - 0.01];
icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', addButtonPos, 'string', '+', 'tag', 'add_cov_button', 'fontsize',...
    UI_FS - 1, 'callback', {@addCov, InputHandle});
icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', removeButtonPos, 'string', '-', 'tag', 'remove_cov_button', 'fontsize',...
    UI_FS - 1, 'callback', {@removeCov, InputHandle});


promptPos = listboxPos;

%% Add cancel, save and run buttons
cancelPos = [0.25 - 0.5*okWidth, promptPos(2) - yOffset - 0.5*okHeight, okWidth, okHeight];
cancelPos(2) = cancelPos(2) - 0.5*cancelPos(4);
icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', cancelPos, 'string', 'Cancel', 'tag', 'cancel_button', 'fontsize',...
    UI_FS - 1, 'callback', 'delete(gcbf);');

okPos = [0.75 - 0.5*okWidth, promptPos(2) - yOffset - 0.5*okHeight, okWidth, okHeight];
okPos(2) = okPos(2) - 0.5*okPos(4);
icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', okPos, 'string', 'Create', 'tag', 'create_button', 'fontsize',...
    UI_FS - 1, 'callback', {@create_design_matrix, InputHandle});

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


function addCov(hObject, event_data, figH)
%% Add Covariates
%

icatb_defaults;
global UI_FS;

figureTag = 'add_cov_mancovan';
covFigHandle = findobj('tag', figureTag);
if (~isempty(covFigHandle))
    delete(covFigHandle);
end

mancovanInfo = get(figH, 'userdata');

listH = findobj(figH, 'tag', 'cov');

covName = '';
transformationName = '';
covVals = '';
cov_type = 'continuous';
covTypes = {'Continuous', 'Categorical'};
covTypeVal = 1;

if (listH == hObject)
    if (~strcmpi(get(figH, 'selectionType'), 'open'))
        return;
    end
    val = get(listH, 'value');
    try
        covName = mancovanInfo.userInput.cov(val).name;
        covVals = mancovanInfo.userInput.cov(val).value;
        if (isnumeric(covVals))
            covVals = covVals(:);
            covVals = num2str(covVals);
        end
        covVals = cellstr(covVals);
        covVals = covVals(:)';
        transformationName = mancovanInfo.userInput.cov(val).transformation;
        cov_type = mancovanInfo.userInput.cov(val).type;
        covTypeVal = strmatch(lower(cov_type), lower(covTypes), 'exact');
    catch
    end
end

covFigHandle = icatb_getGraphics('Select Covariates', 'normal', figureTag);
set(covFigHandle, 'menubar', 'none');

promptWidth = 0.52; controlWidth = 0.32; promptHeight = 0.05;
xOffset = 0.02; yOffset = promptHeight; yPos = 0.9;
okWidth = 0.12; okHeight = promptHeight;

%% Covariate name and value
promptPos = [xOffset, yPos - 0.5*yOffset, promptWidth, promptHeight];
textH = icatb_uicontrol('parent', covFigHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Enter Covariate Name', 'fontsize', UI_FS - 1);
icatb_wrapStaticText(textH);

promptPos = get(textH, 'position');

editPos = promptPos;
editPos(1) = editPos(1) + editPos(3) + xOffset;
editPos(3) = controlWidth;
editH = icatb_uicontrol('parent', covFigHandle, 'units', 'normalized', 'style', 'edit', 'position', editPos, 'string', covName, 'tag', 'cov_name', 'fontsize', UI_FS - 1);

%% Type of covariate (Continuous or categorical)
promptPos = [xOffset, editPos(2) - 0.5*promptHeight - yOffset, promptWidth, promptHeight];
textH = icatb_uicontrol('parent', covFigHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Select Type Of Covariate', 'fontsize', UI_FS - 1);
icatb_wrapStaticText(textH);

promptPos = get(textH, 'position');

editPos = promptPos;
editPos(1) = editPos(1) + editPos(3) + xOffset;
editPos(3) = controlWidth;
covH = icatb_uicontrol('parent', covFigHandle, 'units', 'normalized', 'style', 'popup', 'position', editPos, 'string', {'Continuous', 'Categorical'}, 'tag', 'cov_type', 'fontsize', UI_FS - 1, 'value', covTypeVal);

promptPos(2) = promptPos(2) - 0.5*promptHeight - yOffset;
promptPos(1) = 0.5 - 0.5*promptWidth;
promptPos(3) = promptWidth;

textH = icatb_uicontrol('parent', covFigHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Add Covariate', 'fontsize', UI_FS - 1);
icatb_wrapStaticText(textH);
promptPos = get(textH, 'position');

editWidth = promptWidth;
editHeight = 0.3;
editPos = [0.5 - 0.5*editWidth, promptPos(2) - yOffset - editHeight, editWidth, editHeight];
editH = icatb_uicontrol('parent', covFigHandle, 'units', 'normalized', 'style', 'edit', 'position', editPos, 'string', covVals, 'fontsize', UI_FS - 1, 'tag', 'cov_value', 'min', 0, 'max', 2, 'callback', {@covValueCallback, figH});

cmenu = uicontextmenu;

set(editH, 'uicontextmenu', cmenu);

uimenu(cmenu, 'Label', 'Load File', 'callback', {@editContextMenuCallback, covFigHandle, figH});


%% Transformation name and value
promptPos(2) = editPos(2) - promptHeight - yOffset;
promptPos(1) = xOffset;
editPos = promptPos;
textH = icatb_uicontrol('parent', covFigHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Enter transformation function (like log, atanh). Leave it as empty if you don''t want to apply transformation.', 'fontsize', UI_FS - 1, ...
    'tag', 'prompt_cov_transformation');
icatb_wrapStaticText(textH);
promptPos = get(textH, 'position');
promptPos(2) = promptPos(2) + (editPos(4) - promptPos(4));
set(textH, 'position', promptPos);
editPos(1) = editPos(1) + editPos(3) + xOffset;
editPos(2) = promptPos(2) + 0.5*promptPos(4) - 0.5*editPos(4);
editPos(3) = controlWidth;
transformH = icatb_uicontrol('parent', covFigHandle, 'units', 'normalized', 'style', 'edit', 'position', editPos, 'string', transformationName, 'tag', 'cov_transformation', 'fontsize', UI_FS - 1);

okPos = [0.5 - 0.5*okWidth, yOffset + 0.5*okHeight, okWidth, okHeight];
icatb_uicontrol('parent', covFigHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', okPos, 'string', 'Done', 'tag', 'cov_done', 'fontsize', UI_FS - 1, 'callback', {@setCovCallback, covFigHandle, figH});


set(findobj(covFigHandle, 'tag', 'cov_type'), 'callback', {@covTypeCallback, covFigHandle, transformH, textH});

covTypeCallback(findobj(covFigHandle, 'tag', 'cov_type'), [], covFigHandle, transformH, textH);


function editContextMenuCallback(hObject, event_data, handles, figH)
%% Context menu callback

mancovanInfo = get(figH, 'userdata');

txtFile = icatb_selectEntry('title', 'Select covariate file' , 'filter', '*.txt;*.asc', 'typeEntity', 'file', 'typeSelection', 'single');
drawnow;
covTypeH = findobj(handles, 'tag', 'cov_type');
opts = cellstr(get(covTypeH, 'string'));
covVal = get(covTypeH, 'value');

try
    val = icatb_mancovan_load_covariates(txtFile, opts{covVal}, mancovanInfo.userInput.numOfSub);
    covValueH = findobj(handles, 'tag', 'cov_value');
    set(covValueH, 'string', val);
catch
    icatb_errorDialog(lasterr, 'Covariate Selection');
end


function setCovCallback(hObject, event_data, covFigH, handles)
%% Set covariate name, value and type

mancovanInfo = get(handles, 'userdata');

covNameH = findobj(covFigH, 'tag', 'cov_name');
covValueH = findobj(covFigH, 'tag', 'cov_value');
covTypeH = findobj(covFigH, 'tag', 'cov_type');
covTransformationH = findobj(covFigH, 'tag', 'cov_transformation');

% Covariate name, value and type
cov_name = get(covNameH, 'string');
cov_value = get(covValueH, 'string');
opts = cellstr(get(covTypeH, 'string'));
val = get(covTypeH, 'value');
covType = lower(opts{val});
cov_transformation = get(covTransformationH, 'string');

try
    if (isempty(cov_name))
        error('Covariate name is not entered');
    end
    
    if (isempty(cov_value))
        error('Covariate vector is not entered');
    end
    
    if (length(mancovanInfo.userInput.cov) > 0)
        chk = strmatch(lower(cov_name), lower(cellstr(char(mancovanInfo.userInput.cov.name))), 'exact');
        if (~isempty(chk))
            ind = chk;
        end
    end
    
    if (~exist('ind', 'var'))
        ind = length(mancovanInfo.userInput.cov) + 1;
    end
    
    mancovanInfo.userInput.cov(ind).name = cov_name;
    mancovanInfo.userInput.cov(ind).value = deblank(cov_value(:)');
    mancovanInfo.userInput.cov(ind).transformation = lower(cov_transformation);
    mancovanInfo.userInput.cov(ind).type = covType;
    
    set(handles, 'userdata', mancovanInfo);
    
    covListH = findobj(handles, 'tag', 'cov');
    set(covListH, 'string', cellstr(char(mancovanInfo.userInput.cov.name)));
    delete(covFigH);
    
catch
    icatb_errorDialog(lasterr, 'Covariate Selection');
end


function removeCov(hObject, event_data, figH)
%% Remove Covariate
%

mancovanInfo = get(figH, 'userdata');
listH = findobj(figH, 'tag', 'cov');
val = get(listH, 'value');

if (~isempty(val))
    check = icatb_questionDialog('title', 'Remove Covariate', 'textbody', 'Do you want to remove the covariate from the list?');
    if (~check)
        return;
    end
end

try
    strs = cellstr(char(mancovanInfo.userInput.cov.name));
    mancovanInfo.userInput.cov(val) = [];
    strs(val) = [];
    set(listH, 'value', 1);
    set(listH, 'string', strs);
    set(figH, 'userdata', mancovanInfo);
catch
end

function covValueCallback(hObject, event_data, figH)
%% Covariate value callback
%

mancovanInfo = get(figH, 'userdata');
val = cellstr(get(hObject, 'string'));
inds = icatb_good_cells(val);
val = val(inds);
if (length(val) == 1)
    val = textscan(val{1}, '%s', 'delimiter', '\t,', 'multipleDelimsAsOne', 1);
    val = val{1};
end
eval('val2 = val;');
%eval(['val2 = [', [val{:}], '];']);

if (length(val2) ~= numel(val2))
    error('Covariate must be entered in a  vector');
end

if (length(val2) ~= mancovanInfo.userInput.numOfSub)
    error(['Covariate vector length must equal the no. of subjects (', num2str(mancovanInfo.userInput.numOfSub), ')']);
end

% if (isnumeric(val2))
%     val2 = val2(:);
%     val2 = num2str(val2);
% end

val2 = strtrim(cellstr(val2));

set(hObject, 'string', val2);

function saveCallback(hObject, event_data, handles)
%% Save the mancovanInfo

mancovanInfo = get(handles, 'userdata');
fileN = fullfile(mancovanInfo.userInput.outputDir, [mancovanInfo.userInput.prefix, '.mat']);
icatb_save(fileN, 'mancovanInfo');
wH = icatb_dialogBox('title', 'File Saved', 'textBody', ['Mancovan information is saved in file ', fileN], 'textType', 'large');
waitfor(wH);
delete(handles);


function covTypeCallback(hObject, ed, handles, transformH, promptH)
%% Covariates Type callback
%

str = cellstr(get(hObject, 'string'));
str = str{get(hObject, 'value')};
set(promptH, 'visible', 'on');
set(promptH, 'enable', 'on');
set(transformH, 'visible', 'on');
set(transformH, 'enable', 'on');
if (strcmpi(str, 'categorical'))
    set(promptH, 'visible', 'off');
    set(promptH, 'enable', 'off');
    set(transformH, 'visible', 'off');
    set(transformH, 'enable', 'off');
end


function create_design_matrix(hObject, event_data, handles)
%% Make design matrix
%

icatb_defaults;
global FONT_COLOR;
global AXES_COLOR;
global BG_COLOR;

set(handles, 'pointer', 'watch');

mancovanInfo = get(handles, 'userdata');

desCriteria = 'mancova';

try
    desCriteria = mancovanInfo.userInput.designCriteria;
catch
end



try
    %% full design matrix
    mancovanInfo = icatb_mancovan_full_design(mancovanInfo);
catch
    set(handles, 'pointer', 'arrow');
    rethrow(lasterror);
end

if (isfield(mancovanInfo, 'time'))
    mancovanInfo = rmfield(mancovanInfo, 'time');
end

%% Model time covariates
if (strcmpi(desCriteria, 'mancova'))
    timeMsgStr = 'Do you want to model time? At most paired samples are only modeled. Also, covariates are assumed to be the same between timepoints.';
    subjectStr = cell(length(mancovanInfo.good_sub_inds), 1);
    count = 0;
    for nSub = mancovanInfo.good_sub_inds
        count = count + 1;
        subjectStr{count} = ['Subject ', num2str(nSub)];
    end
    timeCovariates = icatb_questionDialog('title', 'Model Time?', 'textbody', timeMsgStr);
    if (timeCovariates)
        listMsgStr = 'Select subjects for time 1';
        [groupVal1, value] = icatb_listdlg('promptstring', listMsgStr, 'liststring', subjectStr, 'title_fig', listMsgStr);
        if (isempty(groupVal1))
            error('Time 1 subjects are not selected');
        end
        
        listMsgStr = 'Select subjects for time 2';
        [groupVal2, value] = icatb_listdlg('promptstring', listMsgStr, 'liststring', subjectStr, 'title_fig', listMsgStr);
        if (isempty(groupVal2))
            error('Time 2 subjects are not selected');
        end
        if (length(groupVal1) ~= length(groupVal2))
            error('Please select equal number of subjects across time');
        end
        mancovanInfo.time.subjects = {groupVal1, groupVal2};
        mancovanInfo.X = mancovanInfo.X(groupVal1, :);
    end
end


%% Save file
fileN = fullfile(mancovanInfo.outputDir, [mancovanInfo.prefix, '.mat']);
icatb_save(fileN, 'mancovanInfo');

set(handles, 'pointer', 'arrow');

if (strcmpi(desCriteria, 'mancova'))
    
    disp(['Design information is stored in field X in file ', fileN]);
    disp('Setup features for mancovan using the same file.');
    fprintf('\n');
    delete(handles);
    
    drawnow;
    
    InputH = figure('name', 'Mancovan Design Matrix', 'color', BG_COLOR);
    set(InputH, 'resize', 'on');
    sh = subplot(1, 1, 1);
    set(sh, 'color', AXES_COLOR);
    imagesc((1:size(mancovanInfo.X, 2)), (1:size(mancovanInfo.X, 1)), mancovanInfo.X);
    colormap(gray);
    axis(sh, 'square');
    xlabel('Covariates', 'parent', sh);
    ylabel('Subjects', 'parent', sh);
    all_regress = [{'Constant'}, mancovanInfo.regressors];
    set(sh, 'TickDir', 'out');
    set(sh, 'XTiCk', 1:length(all_regress));
    set(sh, 'XTickLabel', all_regress);
    set(sh, 'YColor', FONT_COLOR, 'XColor', FONT_COLOR);
    
    
    
    
else
    
    disp(['Design information is stored in field ttestOpts in file ', fileN]);
    disp('Setup features for ttest using the same file.');
    fprintf('\n');
    delete(handles);
    drawnow;
    
end

function designCriteriaCallback(hObject, event_data, handles)
%% Design criteria callback
%

covH = findobj(handles, 'tag', 'cov');
addH = findobj(handles, 'tag', 'add_cov_button');
removeH = findobj(handles, 'tag', 'remove_cov_button');

val = get(hObject, 'value');
strs = cellstr(get(hObject, 'string'));

mancovanInfo = get(handles, 'userdata');

desCriteria = strs{val};

mancovanInfo.userInput.designCriteria = desCriteria;

if (~strcmpi(strs{val}, 'mancova'))
    set(covH, 'enable', 'off');
    set(addH, 'enable', 'off');
    set(removeH, 'enable', 'off');
else
    set(covH, 'enable', 'on');
    set(addH, 'enable', 'on');
    set(removeH, 'enable', 'on');
    set(handles, 'userdata', mancovanInfo);
    return;
end


% avgRuns = 1;
% if (mancovanInfo.userInput.numOfSess ~= 1)
%     if (strcmpi(desCriteria, 'paired t-test'))
%         avgRuns = icatb_questionDialog('title', 'Average Runs?', 'textbody', 'Do You Want To Average Runs?');
%     end
% end

avgRuns = 1;
if (mancovanInfo.userInput.numOfSess ~= 1)
    if (strcmpi(desCriteria, 'paired t-test'))
        avgRuns = 0;
    end
end

if (~avgRuns)
    subjectStr = cell( mancovanInfo.userInput.numOfSub* mancovanInfo.userInput.numOfSess, 1);
    count = 0;
    for nSub = 1: mancovanInfo.userInput.numOfSub
        for nSess = 1: mancovanInfo.userInput.numOfSess
            count = count + 1;
            subjectStr{count} = ['Subject ', num2str(nSub), ' Session ', num2str(nSess)];
        end
    end
else
    subjectStr = cell( mancovanInfo.userInput.numOfSub, 1);
    count = 0;
    for nSub = 1:mancovanInfo.userInput.numOfSub
        count = count + 1;
        subjectStr{count} = ['Subject ', num2str(nSub)];
    end
end

ttestOpts.des = desCriteria;
%ttestOpts.avgRuns = avgRuns;

if (strcmpi(desCriteria, 'one sample t-test'))
    groupName = '';
    groupVal = [];
    try
        groupName = mancovanInfo.userInput.ttestOpts.t.name{1};
    catch
    end
    try
        groupVal = mancovanInfo.userInput.ttestOpts.t.val{1};
    catch
    end
    [groupName, groupVal] = icatb_select_groups_gui(subjectStr, 'Group', 'selGrp', groupName, groupVal);
    ttestOpts.t.name = {groupName};
    ttestOpts.t.val = {groupVal};
    
elseif (strcmpi(strs{val}, 'two sample t-test'))
    
    group1Name = '';
    groupVal1 = [];
    try
        group1Name = mancovanInfo.userInput.ttestOpts.t.name{1};
    catch
    end
    try
        groupVal1 = mancovanInfo.userInput.ttestOpts.t.val{1};
    catch
    end
    
    group2Name = '';
    groupVal2 = [];
    try
        group2Name = mancovanInfo.userInput.ttestOpts.t.name{2};
    catch
    end
    try
        groupVal2 = mancovanInfo.userInput.ttestOpts.t.val{2};
    catch
    end
    
    [group1Name, groupVal1] = icatb_select_groups_gui(subjectStr, 'Group 1', 'selGrp1', group1Name, groupVal1);
    [group2Name, groupVal2] = icatb_select_groups_gui(subjectStr, 'Group 2', 'selGrp2', group2Name, groupVal2);
    
    ttestOpts.t.name = {group1Name, group2Name};
    ttestOpts.t.val = {groupVal1, groupVal2};
    
elseif (strcmpi(strs{val}, 'paired t-test'))
    
    group1Name = '';
    groupVal1 = [];
    try
        group1Name = mancovanInfo.userInput.ttestOpts.t.name{1};
    catch
    end
    try
        groupVal1 = mancovanInfo.userInput.ttestOpts.t.val{1};
    catch
    end
    
    group2Name = '';
    groupVal2 = [];
    try
        group2Name = mancovanInfo.userInput.ttestOpts.t.name{2};
    catch
    end
    try
        groupVal2 = mancovanInfo.userInput.ttestOpts.t.val{2};
    catch
    end
    
    [group1Name, groupVal1] = icatb_select_groups_gui(subjectStr, 'Condition 1', 'selGrp1', group1Name, groupVal1);
    [group2Name, groupVal2] = icatb_select_groups_gui(subjectStr, 'Condition 2', 'selGrp2', group2Name, groupVal2);
    
    if (length(groupVal1) ~= length(groupVal2))
        error('Please select same number of data-sets for both conditions');
    end
    
    ttestOpts.t.name = {group1Name, group2Name};
    ttestOpts.t.val = {groupVal1, groupVal2};
    
end

mancovanInfo.userInput.ttestOpts = ttestOpts;
set(handles, 'userdata', mancovanInfo);

