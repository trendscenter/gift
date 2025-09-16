function icatb_statistical_testing_TC(regressionParamFile, userInputFile)
%% Function to do statistical testing of time courses (beta weights)
%
% Inputs:
% 1. regressionParamFile - Regression parameter text file
% 2. userInputFile - Input file containing necessary parameters to do stats
%
%

icatb_defaults;
global PARAMETER_INFO_MAT_FILE;
global UI_FS;

%% Select regression parameters text file
regressFilterName = '*_temporal_regression';
if ~exist('regressionParamFile', 'var')
    regressionParamFile = icatb_selectEntry('typeEntity', 'file', 'typeSelection', 'single', 'filter', ...
        [regressFilterName, '.txt'], 'title', 'Select regression parameters text file ...');
end

drawnow;

if isempty(regressionParamFile)
    error('Regression parameters text file is not selected');
end

if ~exist('userInputFile', 'var')
    userInputFile = '';
end

%% Initialise vars
groupData = [];
contrasts = [];
designCriteria = 1;
selConditions = 1;
groupsStr = '';
selGroups = 1;
contrastStr = '';
selContrasts = 1;

%% Design criteria options
desCriteriaOptions = {'One Sample t-test', 'Two Sample t-test', 'One Way Anova (Groups)', 'One Way Anova (Regressors)', ...
    'Two Way Anova (Groups & Regressors)', 'Multiple Regression', 'Paired t-test'};

helpMsg = 'Parsing regression parameters text file ...';

%% Get information from inputfile
if ~isempty(userInputFile)
    
    %% Check file path
    userInpPath = fileparts(userInputFile);
    
    if isempty(userInpPath)
        userInpPath = fileparts(which(userInputFile));
    end
    
    if isempty(userInpPath)
        error('Error:InputFile', 'File %s doesn''t exist\n', userInputFile);
    end
    
    %% Variables from input file
    keywd = 'averageRuns';
    inputData = icatb_read_variables(userInputFile, keywd, 'scalar', 'integer');
    averageRuns = getfield(inputData, keywd);
    clear inputData;
    
    keywd = 'desCriteria';
    inputData = icatb_read_variables(userInputFile, keywd, 'scalar', 'integer');
    designCriteria = getfield(inputData, keywd);
    clear inputData;
    
    if (designCriteria > length(desCriteriaOptions))
        error('Error:DesCriteria', 'Value for desCriteria (%s) exceeds the number of design criteria options (%s)', ...
            num2str(designCriteria), num2str(length(desCriteriaOptions)));
    end
    
    if (designCriteria == 6)
        keywd = 'multi_regress_files';
        inputData = icatb_read_variables(userInputFile, keywd, 'vector', 'cell');
        multi_regress_files = getfield(inputData, keywd);
        clear inputData;
    end
    
    % Group information
    keywd = 'groupInfo';
    inputData = icatb_read_variables(userInputFile, keywd, 'vector', 'cell');
    groupInfo = getfield(inputData, keywd);
    clear inputData;
    
    % Selected groups
    keywd = 'selGroups';
    inputData = icatb_read_variables(userInputFile, keywd, 'vector', 'integer');
    selGroups = getfield(inputData, keywd);
    clear inputData;
    
    
    if max(selGroups) > size(groupInfo, 1)
        error('Error:SelGroups', 'Maximum number specified for selGroups (%s) exceeds the length of groupInfo (%s) variable', ...
            num2str(max(selGroups)), num2str(size(groupInfo, 1)));
    end
    
    % Selected conditions
    keywd = 'selConditions';
    inputData = icatb_read_variables(userInputFile, keywd, 'vector', 'integer');
    selConditions = getfield(inputData, keywd);
    clear inputData;
    
    contrastNames = '';
    contrastMatrix = [];
    
    try
        % Contrast names
        keywd = 'contrastNames';
        inputData = icatb_read_variables(userInputFile, keywd, 'vector', 'cell');
        contrastNames = getfield(inputData, keywd);
        clear inputData;
    catch
    end
    
    if exist('contrastNames', 'var') && ~isempty(contrastNames)
        % Contrast matrix
        keywd = 'contrastMatrix';
        inputData = icatb_read_variables(userInputFile, keywd, 'vector', 'numeric');
        contrastMatrix = getfield(inputData, keywd);
        clear inputData;
    end
    
    contrastStr = '';
    selContrasts = [];
    
    
    if exist('contrastNames', 'var')
        
        if length(contrastNames) ~= size(contrastMatrix, 1)
            error('Error:CheckCon', 'Length of contrastNames (%s) doesn''t match number of contrasts (%s)', num2str(length(contrastNames)), ...
                num2str(size(contrastMatrix, 1)));
        end
        
        contrastStr = str2mat(contrastNames);
        selContrasts = [1:size(contrastStr, 1)];
        
    end
    
    groupsStr = str2mat(groupInfo(:, 1));
    
    
else
    
    helpHandle = helpdlg(helpMsg, helpMsg);
    
end
% End for checking user input


disp(helpMsg);

%% Parse text file
try
    parameters = icatb_parseRegressParamTextFile(regressionParamFile);
    if exist('helpHandle', 'var')
        if ishandle(helpHandle)
            delete(helpHandle);
        end
    end
catch
    if exist('helpHandle', 'var')
        if ishandle(helpHandle)
            delete(helpHandle);
        end
    end
    icatb_displayErrorMsg;
end

disp('Done parsing regression parameters text file');

[outputDir, fName] = fileparts(regressionParamFile);

if isempty(outputDir)
    outputDir = pwd;
end

cd(outputDir);

prefix = strrep(fName, regressFilterName(2:end), '');

% Find parameter file
spmMatFlag = 'Unspecified';
paramFile = fullfile(outputDir, [prefix, PARAMETER_INFO_MAT_FILE, '.mat']);
if exist(paramFile, 'file')
    load(paramFile);
    if isfield(sesInfo.userInput, 'spmMatFlag')
        spmMatFlag = sesInfo.userInput.spmMatFlag;
    end
    clear sesInfo;
end

% Get data from parameters structure
viewingSet = parameters.viewingSet;
numOfSub = parameters.numOfSub;
numOfSess = parameters.numOfSess;
componentNum = parameters.componentNum;
regressInfo = parameters.regressInfo;
numConditions = parameters.numConditions;


if ~exist('averageRuns', 'var')
    
    averageRuns = 0;
    
    if (numOfSub > 1) && (numOfSess > 1)
        msgStr = 'Do you want to average across sessions?';
        [averageRuns] = icatb_questionDialog('title', 'Average across sessions', 'textbody', msgStr);
    end
    
end

%%%% Average runs or sessions %%%%
if averageRuns
    
    disp('Averaging beta weights across sessions ...');
    countDataSet = 0;
    
    % % Initialise regressInfo variable
    newRegressInfo = repmat(struct('cond', repmat(struct('name', '', 'values', []), numConditions, 1)), ...
        numOfSub, 1);
    
    % Loop over subjects
    for nSub = 1:numOfSub
        % Get the condition and regressor name
        condNames = str2mat(regressInfo(countDataSet + 1).cond.name);
        regressorNames = str2mat(regressInfo(countDataSet + 1).regressorName.name);
        % Loop over sessions
        for nSess = 1:numOfSess
            countDataSet = countDataSet + 1;
            if nSess == 1
                paramValues = zeros(length(componentNum), numConditions, numOfSess);
            end
            % Loop over conditions
            for nCond = 1:numConditions
                paramValues(:, nCond, nSess) = regressInfo(countDataSet).cond(nCond).values(:);
            end
            % End loop over conditions
        end
        % End loop over sessions
        
        % Assign the mean of runs to new regressInfo structure
        paramValues = icatb_nanmean(paramValues, 3);
        
        % Loop over conditions
        for nCond = 1:numConditions
            newRegressInfo(nSub).cond(nCond).values = paramValues(:, nCond);
            newRegressInfo(nSub).cond(nCond).name = deblank(condNames(nCond, :));
            newRegressInfo(nSub).regressorName(nCond).name = deblank(regressorNames(nCond, :));
        end
        % End loop over conditions
        
        clear paramValues;
        
    end
    % End loop over subjects
    
    clear regressInfo;
    regressInfo = newRegressInfo;
    clear newRegressInfo;
    
    viewingSet = 'Session';
    numOfSess = 1;
    
    disp('Done averaging beta weights across sessions');
    fprintf('\n');
end
%%%% End for averaging across runs or sessions %%%


parameters.regressInfo = regressInfo;

if exist('multi_regress_files', 'var')
    parameters.multiple_regression.multi_regress_files = multi_regress_files;
end

if ~isempty(userInputFile)
    
    if max(selConditions) > numConditions
        error('Error:selConditions', 'Maximum value for selConditions (%s) exceeds the number of conditions (%s)', ...
            num2str(max(selConditions)), num2str(numConditions));
    end
    
    groupData = repmat(struct('name', '', 'subVal', [], 'sessVal', []), size(groupInfo, 1), 1);
    % Group data structure
    for nG = 1:length(groupData)
        groupData(nG).name = groupInfo{nG, 1};
        if ~isnumeric(groupInfo{nG, 2})
            groupData(nG).subVal = str2num(groupInfo{nG, 2});
        else
            groupData(nG).subVal = groupInfo{nG, 2};
        end
        
        % Do error check for number of subjects
        if max(groupData(nG).subVal) > numOfSub
            error('Error:groupInfo', ['Check groupInfo variable as maximum of the subjects (%s) entered \n', ...
                'exceeds the number of subjects (%s) for group (%s)'], ...
                num2str(max(groupData(nG).subVal)), num2str(numOfSub), groupData(nG).name);
        end
        % End for doing error check on subjects
        
        if numOfSess == 1
            groupData(nG).sessVal = 1;
        else
            if ~isnumeric(groupInfo{nG, 3})
                groupData(nG).sessVal = str2num(groupInfo{nG, 3});
            else
                groupData(nG).sessVal = groupInfo{nG, 3};
            end
        end
        
        groupData(nG).subVal = [groupData(nG).subVal(:)]';
        groupData(nG).sessVal = [groupData(nG).sessVal(:)]';
        
        % Do error check for number of sessions
        if max(groupData(nG).sessVal) > numOfSess
            error('Error:groupInfo', ['Check groupInfo variable as maximum of the sessions (%s) entered \n', ...
                'exceeds the number of sessions (%s) for group (%s)'], ...
                num2str(max(groupData(nG).sessVal)), num2str(numOfSess), groupData(nG).name);
        end
        % End for doing error check on sessions
        
    end
    
    
    contrasts = repmat(struct('name', '', 'vector', []), size(contrastMatrix, 1), 1);
    % Contrast structure
    for nC = 1:size(contrastMatrix, 1)
        contrasts(nC).name = contrastNames{nC};
        contrasts(nC).vector = num2str(contrastMatrix(nC, :));
        if sum(contrastMatrix(nC, :)) ~= 0
            error('Error:ConVec', 'Sum of contrasts of %s is not equal to zero', contrasts(nC).name);
        end
    end
end

% Delete all handles with tag stat_testing_tc
deleteHandles('stat_testing_tc');


% Plot Figure
graphicsHandle = icatb_getGraphics('Statistical Testing of Time courses (Beta Weights)', 'displaygui', ...
    'stat_testing_tc', 'off');

set(graphicsHandle, 'menubar', 'none');


% Plot a figure with the following:
% a. Select criteria: One sample, Two sample, Anova (N-way)
% b. Enter no. of groups.
% c. Draw Groups Listbox
% d. Draw Subjects Listbox with an ed button below:
% e. Draw Conditions Listbox.

% Subject string
subjectStr = repmat(struct('name', ''), numOfSub, 1);
sessionStr = repmat(struct('name', ''), numOfSess, 1);

% Loop over subjects
for nSub = 1:numOfSub
    subjectStr(nSub).name = ['Subject ', num2str(nSub)];
end
% End loop over subjects

subjectString = str2mat(subjectStr.name);

% Loop over sessions
for nSess = 1:numOfSess
    sessionStr(nSess).name = ['Session ', num2str(nSess)];
end
% End loop over sessions

sessionString = str2mat(sessionStr.name);

% Offsets
xOffset = 0.04; yOffset = 0.04;
popupTextHeight = 0.05; popupTextWidth = 0.45; yPos = 0.97;
popupTextPos = [xOffset, yPos - yOffset - popupTextHeight, popupTextWidth, popupTextHeight];


%%%%%%% Select design %%%%%%%%%

% Plot Text
designTextH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'text', ...
    'position', popupTextPos, 'String', 'Select design criteria', 'fontsize', UI_FS - 1, ...
    'horizontalalignment', 'center');

[designTextH] = icatb_wrapStaticText(designTextH);

popupTextPos = get(designTextH, 'position');

popupWidth = 0.4;
popupPos = popupTextPos;
popupPos(1) = popupPos(1) + popupPos(3) + xOffset;
popupPos(3) = popupWidth;

% Plot popup
designPopupH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'popup', ...
    'position', popupPos, 'String', desCriteriaOptions, 'fontsize', UI_FS - 1, 'tag', 'design_criteria', 'callback', ...
    {@popupCallback, graphicsHandle}, 'value', designCriteria);

%%%%%%% End for selecting design %%%%%%%%

listboxWidth = 0.3;

%%%%%%%%%%%%%% Groups Text and Listbox %%%%%%%%%%%%%%%%%%%%

% Plot textbox
textBoxWidth = listboxWidth; textboxHeight = 0.04;
textboxPos = [xOffset, popupTextPos(2) - popupTextPos(4) - yOffset, textBoxWidth, textboxHeight];

% Plot groups text
groupsTextH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'text', ...
    'position', textboxPos, 'String', {'Groups'}, 'fontsize', UI_FS - 1);

[groupsTextH] = icatb_wrapStaticText(groupsTextH);

textboxPos = get(groupsTextH, 'position');

listboxPos(1) = textboxPos(1); listboxPos(3) = textboxPos(3);
listboxPos(4) = 0.25;
listboxPos(2) = textboxPos(2) - 0.5*textboxPos(4) - listboxPos(4);

% Plot groups listbox
groupsListH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'listbox', ...
    'position', listboxPos, 'String', groupsStr, 'value', selGroups, 'min', 0, 'max', 2, 'tag', ...
    'selGroups', 'fontsize', UI_FS - 1, 'callback', {@doubleClickGroups, graphicsHandle});

%%%%%%%%%%%%%%% End for plotting groups text and listbox %%%%%%%%%%%%


%  Draw three push buttons right of listbox
buttonWidth = 0.08; buttonHeight = 0.04;

addPos = [listboxPos(1) + listboxPos(3) + xOffset, listboxPos(2) + listboxPos(4) - 0.01 - buttonHeight, buttonWidth, buttonHeight];
removePos = addPos;
removePos(2) = removePos(2) - removePos(4) - 0.02;
donePos = removePos;
donePos(2) = donePos(2) - donePos(4) - 0.02;

% Plot Add button
addGroupsH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'pushbutton', ...
    'position', addPos, 'String', '+', 'tag',  'addGroups', 'fontsize', UI_FS - 1, 'callback', {@addGroupsCallback, ...
    graphicsHandle});

% Plot Add button
removeGroupsH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'pushbutton', ...
    'position', removePos, 'String', '-', 'tag',  'removeGroups', 'fontsize', UI_FS - 1, 'callback', ...
    {@removeGroupsCallback, graphicsHandle});

%%%%%%%%%%%%%% Condition Text and Listbox %%%%%%%%%%%%%%%%%%%%

% Plot textbox
textBoxWidth = listboxWidth; textboxHeight = 0.04;
textboxPos(1) = addPos(1) + addPos(3) + xOffset; %popupPos(1);

% Plot conditions text
condTextH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'text', ...
    'position', textboxPos, 'String', {'Regressors'}, 'fontsize', UI_FS - 1);

[condTextH] = icatb_wrapStaticText(condTextH);

textboxPos = get(condTextH, 'position');

listboxPos(1) = textboxPos(1); listboxPos(3) = textboxPos(3);
listboxPos(2) = textboxPos(2) - 0.5*textboxPos(4) - listboxPos(4);


regressorsListStr = cellstr(str2mat(regressInfo(1).regressorName.name));
if strcmpi(spmMatFlag, 'diff_sub_diff_sess') || strcmpi(spmMatFlag, 'same_sub_diff_sess')
    regressorsListStr = regexprep(regressorsListStr, 'Sn\(\d+\) ', '');
end
regressorsListStr = str2mat(regressorsListStr);

% Plot condition listbox
condListH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'listbox', ...
    'position', listboxPos, 'String', regressorsListStr, 'value', selConditions, 'min', 0, 'max', ...
    2, 'tag', 'selConditions', 'fontsize', UI_FS - 1);

%%%%%%%%%%%%%%% End for plotting conds text and listbox %%%%%%%%%%%%


%%%%%%%%%%%%%% Other Conditions Text and Listbox %%%%%%%%%%%%%%%%%%%%

% Plot textbox
textBoxWidth = listboxWidth; textboxHeight = 0.04;
groupTextPos = get(groupsTextH, 'position');
textboxPos = groupTextPos;
textboxPos(2) = listboxPos(2) - yOffset - textboxHeight;

% Plot other conditions text
otherCondTextH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'text', ...
    'position', textboxPos, 'String', {'Other Conditions'}, 'fontsize', UI_FS - 1, 'visible', 'off');

[otherCondTextH] = icatb_wrapStaticText(otherCondTextH);

textboxPos = get(otherCondTextH, 'position');

listboxPos(1) = textboxPos(1); listboxPos(3) = textboxPos(3);
listboxPos(2) = textboxPos(2) - 0.5*textboxPos(4) - listboxPos(4);

% Plot other conditions listbox
otherCondListH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'listbox', ...
    'position', listboxPos, 'String', '', 'value', 1, 'min', 0, 'max', 2, 'tag', ...
    'selOtherConditions', 'fontsize', UI_FS - 1, 'callback', {@doubleClickOtherCond, graphicsHandle}, 'visible', 'off');

%%%%%%%%%%%%%%% End for plotting other conds text and listbox %%%%%%%%%%%%


%  Draw three push buttons right of listbox

addPos = [listboxPos(1) + listboxPos(3) + xOffset, listboxPos(2) + listboxPos(4) - 0.01 - buttonHeight, buttonWidth, buttonHeight];
removePos = addPos;
removePos(2) = removePos(2) - removePos(4) - 0.02;
donePos = removePos;
donePos(2) = donePos(2) - donePos(4) - 0.02;

% Plot Add button
addOtherCondH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'pushbutton', ...
    'position', addPos, 'String', '+', 'tag',  'addOtherConditions', 'fontsize', UI_FS - 1, 'enable', 'on', ...
    'callback', {@addConditions, graphicsHandle}, 'visible', 'off');

% Plot Add button
removeOtherCondH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'pushbutton', ...
    'position', removePos, 'String', '-', 'tag',  'removeOtherConditions', 'fontsize', UI_FS - 1, ...
    'enable', 'on', 'callback', {@removeCondCallback, graphicsHandle}, 'visible', 'off');


%%%%%%%%%%%%%% Contrasts %%%%%%%%%%%%%%%%%%%%

% Plot textbox
textBoxWidth = listboxWidth; textboxHeight = 0.04;
textboxPos(1) = addPos(1) + addPos(3) + xOffset;

% Plot contrasts text
contrastsTextH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'text', ...
    'position', textboxPos, 'String', {'Anova Contrasts'}, 'fontsize', UI_FS - 1);

[contrastsTextH] = icatb_wrapStaticText(contrastsTextH);

textboxPos = get(contrastsTextH, 'position');

listboxPos(1) = textboxPos(1); listboxPos(3) = textboxPos(3);
listboxPos(2) = textboxPos(2) - 0.5*textboxPos(4) - listboxPos(4);

% Plot contrast listbox
contrastListH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'listbox', ...
    'position', listboxPos, 'String', contrastStr, 'value', selContrasts, 'min', 0, 'max', 2, 'tag', ...
    'selContrasts', 'fontsize', UI_FS - 1, 'callback', {@doubleClickContrast, graphicsHandle});

%%%%%%%%%%%%%%% End for plotting contrasts text and listbox %%%%%%%%%%%%


%  Draw three push buttons right of listbox

addPos = [listboxPos(1) + listboxPos(3) + xOffset, listboxPos(2) + listboxPos(4) - 0.01 - buttonHeight, buttonWidth, buttonHeight];
removePos = addPos;
removePos(2) = removePos(2) - removePos(4) - 0.02;
donePos = removePos;
donePos(2) = donePos(2) - donePos(4) - 0.02;

% Plot Add button
addContrastsH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'pushbutton', ...
    'position', addPos, 'String', '+', 'tag',  'addContrasts', 'fontsize', UI_FS - 1, 'enable', 'on', 'callback', ...
    {@addContrasts, graphicsHandle});

% Plot remove button
removeContrastsH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'pushbutton', ...
    'position', removePos, 'String', '-', 'tag',  'removeContrasts', 'fontsize', UI_FS - 1, ...
    'enable', 'on', 'callback', {@removeContrasts, graphicsHandle});

%%%%%%%%%% Plot Calculate Button %%%%%%%%%%
calculateWidth = 0.2; calculateHeight = 0.05;
%calculatePos = get(condListH, 'position');
calculatePos = [0.75 - 0.5*calculateWidth, listboxPos(2) - yOffset - calculateHeight, ...
    calculateWidth, calculateHeight];
calculateH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'pushbutton', ...
    'position', calculatePos, 'String', 'Calculate', 'tag',  'calculate', 'fontsize', UI_FS - 1, 'callback', ...
    {@calculateCallback, graphicsHandle});

% Parameters structure
parameters.subjectString = subjectString;
parameters.sessionString = sessionString;
parameters.output_prefix = prefix;
parameters.outputDir = outputDir;
parameters.numOfSub = numOfSub;
parameters.numOfSess = numOfSess;

% Set figure data
figureData.parameters = parameters;
figureData.groupData = groupData;
figureData.otherCond = [];
figureData.contrasts = contrasts;
figureData.userInputFile = userInputFile;

set(graphicsHandle, 'userdata', figureData);


% Execute callbacks
popupCallback(findobj(graphicsHandle, 'tag', 'design_criteria'), [], graphicsHandle);

if ~isempty(userInputFile)
    calculateCallback(calculateH, [], graphicsHandle);
    try
        delete(graphicsHandle);
    catch
    end
else
    
    % Help Menu
    helpMenu = uimenu('parent', graphicsHandle, 'label', 'GIFT-Help');
    htmlHelpMenu = uimenu(helpMenu, 'label', 'Stats On Beta Weights', 'callback', 'icatb_openHTMLHelpFile(''icatb_stats_on_beta_weights.htm'');');
    
    % Add parameters to figure data
    set(graphicsHandle, 'visible', 'on');
end



%%%%% Function callbacks %%%%%


%%%%%%%%%%%%% Callbacks Related to Groups %%%%%%%%%
function addGroupsCallback(hObject, event_data, handles, flagAdd)
% Add groups to the listbox

icatb_defaults;
global UI_FS;


if ~exist('flagAdd', 'var')
    flagAdd = 'new';
end


% Get figure data
figureData = get(handles, 'userdata');

% Parameters structure
parameters = figureData.parameters;

% Data-set string
subjectString = parameters.subjectString;
sessionString = parameters.sessionString;

% Selected groups
selGroupsH = findobj(handles, 'tag', 'selGroups');

getListValue = get(selGroupsH, 'val');

% Data of different groups
groupData = figureData.groupData;

if isempty(groupData)
    flagAdd = 'new';
end

if strcmpi(flagAdd, 'new')
    groupNum = length(groupData) + 1;
    subListVal = 1;
    sessListVal = 1;
    thisGroupName = ['Group ', num2str(groupNum)];
else
    groupNum = getListValue;
    subListVal = groupData(groupNum).subVal;
    sessListVal = groupData(groupNum).sessVal;
    thisGroupName = groupData(groupNum).name;
end

% Delete figures that have tag select_subjects_stats
deleteHandles('select_subjects_stats');

% Plot Figure
[graphicsHandle] = icatb_getGraphics('Select subjects', 'displaygui', 'select_subjects_stats', 'off');
set(graphicsHandle, 'menubar', 'none');

if ispc
    set(graphicsHandle, 'windowstyle', 'modal');
end

% Offsets
xOffset = 0.05; yOffset = 0.05;
editTextHeight = 0.05; editTextWidth = 0.45; yPos = 0.95;

editTextPos = [xOffset, yPos - yOffset - editTextHeight, editTextWidth, editTextHeight];

% Plot Text
editTextH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'text', ...
    'position', editTextPos, 'String', ['Enter name for group ', num2str(groupNum)], 'fontsize', UI_FS - 1, ...
    'horizontalalignment', 'center');

[editTextH] = icatb_wrapStaticText(editTextH);

editTextPos = get(editTextH, 'position');

%%% Edit Box
editWidth = 0.4;
editPos = editTextPos;
editPos(1) = editPos(1) + editPos(3) + xOffset;
editPos(3) = editWidth;

% Plot edit control
editH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'edit', ...
    'position', editPos, 'String', thisGroupName, 'fontsize', UI_FS - 1, 'tag', 'group_name');

% List text width and height
listTextWidth = 0.35; listTextHeight = 0.05;
listTextPos = [xOffset, editTextPos(2) - 2*yOffset - 0.5*listTextHeight, listTextWidth, ...
    listTextHeight];

% Plot listbox
listTextH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'text', ...
    'position', listTextPos, 'String', 'Select Subjects', 'fontsize', UI_FS - 1);


[listTextH] = icatb_wrapStaticText(listTextH);

listTextPos = get(listTextH, 'position');

% Listbox position
listWidth = listTextWidth; listHeight = 0.4;
listPos = listTextPos;
listPos(2) = listPos(2) - yOffset - listHeight;
listPos(3) = listWidth;
listPos(4) = listHeight;

% Plot listbox
listH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'listbox', ...
    'position', listPos, 'String', subjectString, 'fontsize', UI_FS - 1, 'tag', 'subject_listbox', 'min', 0, ...
    'max', 2, 'value', subListVal);

editWidth = 0.2;
editHeight = 0.05;
editPos = [listPos(1) + 0.5*listPos(3) - 0.5*editWidth, listPos(2) - yOffset - 0.5*editHeight, ...
    editWidth, editHeight];

% Plot edit control
editH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'edit', ...
    'position', editPos, 'String', num2str(get(listH, 'value')), 'fontsize', UI_FS - 1, 'tag', ...
    'subject_editbox');


set(listH, 'callback', {@listSubjectsCallback, graphicsHandle, editH});
set(editH, 'callback', {@editSubjectsCallback, graphicsHandle, listH})

% List text width and height
listTextWidth = 0.35; listTextHeight = 0.05;
listTextPos = [1 - xOffset - listTextWidth, editTextPos(2) - 2*yOffset - 0.5*listTextHeight, listTextWidth, ...
    listTextHeight];

% Plot listbox
listTextH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'text', ...
    'position', listTextPos, 'String', 'Select Sessions', 'fontsize', UI_FS - 1);


[listTextH] = icatb_wrapStaticText(listTextH);

listTextPos = get(listTextH, 'position');

% Listbox position
listWidth = listTextWidth; listHeight = 0.4;
listPos = listTextPos;
listPos(2) = listPos(2) - yOffset - listHeight;
listPos(3) = listWidth;
listPos(4) = listHeight;

% Plot listbox
listH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'listbox', ...
    'position', listPos, 'String', sessionString, 'fontsize', UI_FS - 1, 'tag', 'session_listbox', 'min', 0, ...
    'max', 2, 'value', sessListVal);

editWidth = 0.2;
editHeight = 0.05;
editPos = [listPos(1) + 0.5*listPos(3) - 0.5*editWidth, listPos(2) - yOffset - 0.5*editHeight, ...
    editWidth, editHeight];

% Plot edit control
editH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'edit', ...
    'position', editPos, 'String', num2str(get(listH, 'value')), 'fontsize', UI_FS - 1, 'tag', ...
    'session_editbox');


set(listH, 'callback', {@listSubjectsCallback, graphicsHandle, editH});
set(editH, 'callback', {@editSubjectsCallback, graphicsHandle, listH})

% Ok button
okWidth = 0.2; okHeight = 0.05;
okPos = [0.5 - 0.5*okWidth, yOffset + 0.5*okHeight, okWidth, okHeight];
okH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'pushbutton', ...
    'position', okPos, 'String', 'OK', 'fontsize', UI_FS - 1, 'tag', 'ok_select_subjects', 'callback', ...
    {@okSubjectsCallback, graphicsHandle, handles, flagAdd});

% Set graphics handle visibility on
set(graphicsHandle, 'visible', 'on');


function okSubjectsCallback(hObject, event_data, graphicsHandle, handles, flagAdd)
% Ok button for selecting subjects


figureData = get(handles, 'userdata');

groupData = figureData.groupData;

% Listbox handle
selGroupsH = findobj(handles, 'tag', 'selGroups');

% Subject listbox
subjectListboxH = findobj(graphicsHandle, 'tag', 'subject_listbox');
sessListboxH = findobj(graphicsHandle, 'tag', 'session_listbox');

groupEditH = findobj(graphicsHandle, 'tag', 'group_name');

% Get data-set numbers
subjectVal = get(subjectListboxH, 'value');

sessVal = get(sessListboxH, 'value');

% Get group name
groupName = deblank(get(groupEditH, 'string'));


if strcmpi(flagAdd, 'new')
    
    % Add data to groupData structure
    newIndex = length(groupData) + 1;
    groupData(newIndex).name = groupName;
    groupData(newIndex).subVal = subjectVal;
    groupData(newIndex).sessVal = sessVal;
    
else
    
    newIndex = get(selGroupsH, 'value');
    groupData(newIndex).name = groupName;
    groupData(newIndex).subVal = subjectVal;
    groupData(newIndex).sessVal = sessVal;
    
end

if ~isempty(groupData)
    checkIndex = strmatch(groupName, str2mat(groupData.name), 'exact');
    if length(checkIndex) > 1
        error('Error:GroupName', 'Group name %s already exists', groupName);
    end
end

% Add field groupData to figureData
figureData.groupData = groupData;

set(handles, 'userdata', figureData);
set(selGroupsH, 'string', str2mat(str2mat(groupData.name)));

delete(graphicsHandle);


function editSubjectsCallback(hObject, event_data, handles, subjectListH)
% Set the listbox value

editString = deblank(get(hObject, 'string'));

editVal = str2num(editString);

numItems = size(get(subjectListH, 'string'), 1);

%%%% Do Error checking %%%%
if isempty(editVal)
    error('Error:EditBox', 'Check the edit box string (%s) \nas it must generate a valid numeric value', editString);
end

if ~isempty(find(editVal == 0))
    error('Error:EditBox', 'Edit box string (%s) cannot accept a value of 0', editString);
end

if numItems < length(editVal)
    error('Error:EditBox', 'Number of items in edit box (%d) is larger than the \nnumber of items (%d) in listbox ', ...
        length(editVal), numItems);
end

if numItems < max(editVal)
    error('Error:EditBox', 'Maximum value in editbox string (%s) is larger than the \nnumber of items (%d) in listbox ', ...
        editString, numItems);
end
%%%% End for doing error checking %%%%

checkInteger = [];
try
    checkInteger = strread(num2str(editVal), '%d');
catch
end
if isempty(checkInteger)
    error('Error:EditBox', 'Editbox string (%s) does not contain integer items', editString);
end

set(subjectListH, 'value', editVal);


function listSubjectsCallback(hObject, event_data, handles, editH)
% Set the listbox value

editString = num2str(get(hObject, 'value'));

set(editH, 'string', editString);

function removeGroupsCallback(hObject, event_data, handles)
%% Remove groups callback
%

% Groups listbox
selGroupsH = findobj(handles, 'tag', 'selGroups');

figureData = get(handles, 'userdata');
groupData = figureData.groupData;

groupVal = get(selGroupsH, 'value');

if ~isempty(groupData)
    msgStr = 'Do you want to remove the selected groups?';
    [removeGroups] = icatb_questionDialog('title', 'Remove Groups', 'textbody', msgStr);
    if removeGroups
        groupData(groupVal) = [];
    else
        return;
    end
end

set(selGroupsH, 'value', 1);

if ~isempty(groupData)
    set(selGroupsH, 'string', str2mat(groupData.name));
else
    set(selGroupsH, 'string', '');
end

% Set figure Data
figureData.groupData = groupData;
set(handles, 'userdata', figureData);



function calculateCallback(hObject, event_data, handles)
%% Calculate callback
%

try
    
    set(handles, 'pointer', 'watch');
    
    %% Get information from figure data
    figureData = get(handles, 'userdata');
    groupData = figureData.groupData;
    otherCond = figureData.otherCond;
    userInputFile = figureData.userInputFile;
    parameters = figureData.parameters;
    regressInfo = parameters.regressInfo;
    outputDir = parameters.outputDir;
    output_prefix = parameters.output_prefix;
    numOfSub = parameters.numOfSub;
    numOfSess = parameters.numOfSess;
    
    %% Selected conditions
    selConditionsH = findobj(handles, 'tag', 'selConditions');
    regressorNames = get(selConditionsH, 'string');
    selectedConditions = get(selConditionsH, 'value');
    
    %% Contrasts
    myContrasts = figureData.contrasts;
    
    clear figureData;
    
    %% Other conditions
    %     otherCondListH = findobj(handles, 'tag', 'selOtherConditions');
    %     if isempty(otherCond)
    %         selectedOtherConditions = [];
    %     else
    %         selectedOtherConditions = get(otherCondListH, 'value');
    %     end
    
    %% Contrasts
    contrastListH = findobj(handles, 'tag', 'selContrasts');
    if isempty(myContrasts)
        selContrasts = [];
    else
        selContrasts = get(contrastListH, 'value');
    end
    
    %% Component numbers
    componentNum = parameters.componentNum;
    
    if isempty(groupData)
        error('Please select groups in order to stats on beta weights');
    end
    
    groupNames = str2mat(groupData.name);
    
    %% Selected groups
    selGroupsH = findobj(handles, 'tag', 'selGroups');
    selectedGroups = get(selGroupsH, 'value');
    
    %%%% Check Design Criteria %%%%%%%%
    designCriteriaH = findobj(handles, 'tag', 'design_criteria');
    designOptions = get(designCriteriaH, 'string');
    designCriteriaVal = get(designCriteriaH, 'value');
    designCriteria = designOptions{designCriteriaVal};
    %%%% End for checking design criteria %%%%%
    
    %% Initial error check
    [selectedGroups, designValue] = errChk(designCriteria, selectedGroups, selectedConditions);
    
    if (designValue == 6)
        %% Check regressors information for multiple regression
        if (~isfield(parameters, 'multiple_regression'))
            error('Please select the regressors if you want to use Multiple Regression');
        end
    end
    
    %% Form equation for one way anova (groups)
    [eq_regressors, eqStr] = formEqRegressors(designValue, selectedConditions, userInputFile, regressorNames);
    
    drawnow;
    
    titlePrint = designCriteria;
    
    designCriteria = lower(designCriteria);
    
    chkAnova = 0;
    if ((designValue == 3) || (designValue == 4) || (designValue == 5))
        chkAnova = 1;
    end
    
    disp(['Doing ', designCriteria, ' on beta weights']);
    
    outFileName = fullfile(outputDir, [output_prefix, '_stats_summary.txt']);
    
    %% Print title
    fid = fopen(outFileName, 'a+');
    
    if fid == -1
        error('Error:FileOpen', 'File %s can''t be opened', outFileName);
    end
    
    fprintf(fid, '%s\n', '.............................................');
    fprintf(fid, '%s\n', '       SUMMARY OF STATS ON BETA WEIGHTS      ');
    fprintf(fid, '%s\n', '.............................................');
    fprintf(fid, '\n');
    fprintf(fid, '\n');
    
    fclose(fid);
    %% End for printing title
    
    %% Loop over selected groups
    for nG = 1:length(selectedGroups)
        if nG == 1
            strPrint1 = ['Selected group/groups is/are: ', deblank(groupNames(selectedGroups(nG), :))];
        else
            strPrint1 = [strPrint1, ', ', deblank(groupNames(selectedGroups(nG), :))];
        end
    end
    %% End loop over selected groups
    
    if ~isempty(eq_regressors)
        strPrint2 = ['Selected regressor is: ', eqStr];
    else
        % Loop over number of regressors
        for nRegress = 1:length(selectedConditions)
            if nRegress == 1
                strPrint2 = ['Selected regressor/regressors is/are: ', deblank(regressorNames(selectedConditions(nRegress), :))];
            else
                strPrint2 = [strPrint2, ', ', deblank(regressorNames(selectedConditions(nRegress), :))];
            end
        end
        % End loop over number of regressors
    end
    
    additionalInfo = str2mat(strPrint1, strPrint2, '');
    
    disp(additionalInfo);
    
    %% Group data
    groupData = groupData(selectedGroups);
    
    %% Get groupID
    [subinGroups, groupID] = get_groupID(groupData, numOfSess);
    
    %% Get data and condID vector
    [data, condInd] = get_regress_data(regressInfo, subinGroups, selectedConditions, componentNum);
    
    clear regressInfo;
    
    if isempty(eq_regressors)
        % Replicate group ID over conditions
        groupID = repmat(groupID, length(selectedConditions), 1);
    else
        %% Use eq_regressors to do operation on conditions
        dat = zeros(length(groupID), length(componentNum), length(selectedConditions));
        for nC = 1:length(selectedConditions)
            inds = find(condInd == nC);
            dat(:, :, nC) = eq_regressors(nC).*data(inds, :);
        end
        clear data;
        data = sum(dat, 3);
        clear dat;
    end
    
    if (chkAnova)
        
        %% Initialise anova results data structure
        modelType = 'linear';
        if (designValue == 3)
            anovaParameters{1} = groupID;
            varNames = {'Groups'};
        elseif (designValue == 4)
            anovaParameters{1} = condInd;
            varNames = {'Regressors'};
        else
            anovaParameters{1} = groupID;
            anovaParameters{2} = condInd;
            modelType = 'interaction';
            varNames = {'Groups'; 'Regressors'};
        end
        
    end
    
    %% Check design criteria
    switch (designValue)
        
        case {1, 2}
            %% For one sample t-test and two sample t-test print component
            % number and p value
            
            %% Initialise p and t values
            pValue = zeros(1, length(componentNum));
            tValue = zeros(1, length(componentNum));
            
            if (designValue == 1)
                %% Loop over number of components
                for nComp = 1:length(componentNum)
                    [pValue(nComp), tValue(nComp)] = icatb_ttest(data(:, nComp));
                end
            else
                %% Loop over number of components
                for nComp = 1:length(componentNum)
                    [pValue(nComp), tValue(nComp)] = icatb_ttest2(data(1:length(subinGroups(1).val), nComp), ...
                        data(length(subinGroups(1).val) + 1 : end, nComp));
                end
            end
            
            numPara = 1;
            varStruct(numPara).tag = 'Component Number';
            varStruct(numPara).value = componentNum(:);
            
            numPara = numPara + 1;
            varStruct(numPara).tag = 'p-value';
            varStruct(numPara).value = pValue(:);
            
            numPara = numPara + 1;
            varStruct(numPara).tag = 'T-value';
            varStruct(numPara).value = tValue(:);
            
            icatb_printToFile(outFileName, varStruct, titlePrint, 'column_wise', 'append', additionalInfo);
            
        case {3, 4, 5}
            %% For anova print table and contrasts table
            %
            
            % Open file for printing
            fid = fopen(outFileName, 'a+');
            
            for nn = 1:size(titlePrint, 1)
                fprintf(fid, '%s\n', titlePrint(nn, :));
            end
            
            fprintf(fid, '\n');
            for nn = 1:size(additionalInfo, 1)
                fprintf(fid, '%s\n', additionalInfo(nn, :));
            end
            
            %% Anova Contrasts
            if (~isempty(selContrasts))
                
                countTerms = 0;
                for nA = 1:length(anovaParameters)
                    countTerms = countTerms + length(unique(anovaParameters{nA}));
                end
                
                if strcmpi(designCriteria, 'one way anova (regressors)')
                    checkConStr = 'conditions';
                elseif strcmpi(designCriteria, 'one way anova (groups)')
                    checkConStr = 'groups';
                else
                    checkConStr = 'conditions and groups';
                end
                
                %% Initialise contrast vector and names
                conVec = zeros(countTerms, length(selContrasts));
                contrastNames = cell(length(selContrasts), 1);
                
                %% Loop over contrasts
                for nCon = 1:length(selContrasts)
                    % Current contrast name and vector
                    contrastName = deblank(myContrasts(selContrasts(nCon)).name);
                    contrastVector = str2num(myContrasts(selContrasts(nCon)).vector);
                    if (length(contrastVector) ~= countTerms)
                        error('Error:Contrasts', ['Contrast vector (%s) doesn''t match the length of ', ...
                            'selected ', checkConStr, ' (%d)'], myContrasts(selContrasts(nCon)).vector, countTerms);
                    end
                    % Store contrast name and vector
                    contrastNames{nCon} = contrastName;
                    conVec(:, nCon) = contrastVector(:);
                end
                %% End loop over contrasts
                
            end
            
            sepStrs = repmat('.', 1, 30);
            
            %% Print summary
            for nComp = 1:length(componentNum)
                
                fprintf(fid, '%s\n', sepStrs);
                fprintf(fid, '%s\n', ['Component ', num2str(nComp), ': ']);
                fprintf(fid, '%s\n', sepStrs);
                
                %% Anova
                anovaR = icatb_anova(data(:, nComp), anovaParameters, 'model', modelType, 'var_names', varNames);
                tbl = icatb_anova_table(anovaR);
                
                icatb_print_table(tbl, fid, 'a+', 0);
                fprintf(fid, '\n');
                
                %% Contrasts
                if (~isempty(selContrasts))
                    %fprintf(fid, '%s\n', sepStrs);
                    [conResults, con_tbl] = icatb_anova_contrasts(anovaR, (1:length(anovaParameters)), conVec, contrastNames);
                    clear conResults;
                    icatb_print_table(con_tbl{1}, fid, 'a+', 0);
                    fprintf(fid, '\n');
                end
                clear con_tbl;
                
            end
            %% End for printing summary
            
            
            fprintf(fid, '\n');
            fprintf(fid, '\n');
            
            fclose(fid);
            
        case 6
            %% Regression
            % Do a multiple regression on regressors like age, test-scores
            % and report the R-square stat for all the regressors
            
            multi_regress_data = parameters.multiple_regression.data;
            multi_regress_names = parameters.multiple_regression.names;
            
            % Number of regressors
            numRegressIn = size(multi_regress_data, 2);
            % Slope names
            slope_names = strcat('Slope (', multi_regress_names, ')');
            slope_names = slope_names';
            % Partial correlations
            partial_corr_names = strcat('Partial Corr (', multi_regress_names, ')');
            partial_corr_names = partial_corr_names';
            
            sepStrs = repmat('.', 1, 30);
            
            % Open file for printing
            fid = fopen(outFileName, 'a+');
            
            for nn = 1:size(titlePrint, 1)
                fprintf(fid, '%s\n', titlePrint(nn, :));
            end
            
            fprintf(fid, '\n');
            for nn = 1:size(additionalInfo, 1)
                fprintf(fid, '%s\n', additionalInfo(nn, :));
            end
            
            
            %% Loop over components
            for nComps = 1:length(componentNum)
                
                fprintf(fid, '%s\n', sepStrs);
                fprintf(fid, '%s\n', ['Component ', num2str(nComps), ': ']);
                fprintf(fid, '%s\n', sepStrs);
                
                % Number of columns and rows of regress table
                numColsRegressTable = 2*(numRegressIn + 1);
                numRowsRegressTable = 1 + (length(selectedConditions)*length(selectedGroups));
                tbl = repmat({''}, numRowsRegressTable, numColsRegressTable);
                % Initial heading
                tbl(1, :) = {'Source', 'R-square', slope_names{:}, partial_corr_names{:}};
                countTbl = 1;
                endDataInd = 0;
                startDataInd = 1;
                % Loop over selected conditions
                for nConds = 1:length(selectedConditions)
                    % Current condition name
                    currentCondName = deblank(regressorNames(selectedConditions(nConds), :));
                    % Loop over selected groups
                    for nGroups = 1:length(selectedGroups)
                        countTbl = countTbl + 1;
                        % Current group name
                        currentGroupName = deblank(groupNames(selectedGroups(nGroups), :));
                        currentRegressTitle = [currentGroupName, ' ', currentCondName];
                        
                        % Current subjects in group
                        currentSubjectsInGroup = subinGroups(nGroups).val;
                        % Ending data indices
                        endDataInd = endDataInd + length(currentSubjectsInGroup);
                        dataInds = (startDataInd:endDataInd);
                        % Remove mean for data and model
                        tmpData = detrend(data(dataInds, nComps), 0);
                        tmpModel = detrend(multi_regress_data(currentSubjectsInGroup, :), 0);
                        
                        % Do multiple regression (Remove only mean )
                        [rSquare_stat, beta_weights, ModelIndices, otherIndices, linearRegress, removeTrend, tmpData, ...
                            subject_partial_corr, partialCorrSlopes] = icatb_multipleRegression(tmpModel, tmpData, numRegressIn, 1, ...
                            size(tmpData, 1), 0);
                        
                        % Truncate the outputs of regression
                        beta_weights = beta_weights(ModelIndices);
                        % Note that GIFT gives absolute of partial corr.
                        % When reporting in the table use the sign of the
                        % slope
                        subject_partial_corr = sign(partialCorrSlopes).*subject_partial_corr;
                        
                        betaStr = cellstr(num2str(beta_weights(:), '%0.4f'));
                        partialCorrStr = cellstr(num2str(subject_partial_corr(:), '%0.4f'));
                        
                        % Store contents to table
                        tbl(countTbl, :) = {currentRegressTitle, num2str(rSquare_stat, '%0.4f'), betaStr{:}, partialCorrStr{:}};
                        
                        % Starting data indices
                        startDataInd = endDataInd + 1;
                    end
                    % End loop over selected groups
                end
                % End loop over selected conditions
                % Print table to file
                icatb_print_table(tbl, fid, 'a+', 0);
                fprintf(fid, '\n');
                
            end
            % End for loop over components
            
            fprintf(fid, '\n');
            fprintf(fid, '\n');
            
            fclose(fid);
            
            
        case 7
            %% Paired t-test
            
            %% Initialise p and t values
            pValue = zeros(1, length(componentNum));
            tValue = zeros(1, length(componentNum));
            
            nLen = ceil(size(data, 1)/2);
            
            if (length(data(1:nLen, 1)) ~= length(data(nLen + 1:end, 1)))
                error('Vector length must be the same between conditions');
            end
            
            for nComp = 1:length(componentNum)
                tmp1 = data(1:nLen, nComp);
                tmp2 = data(nLen+1:end, nComp);
                [pValue(nComp), tValue(nComp)] = icatb_ttest(tmp1 - tmp2);
            end
            
            numPara = 1;
            varStruct(numPara).tag = 'Component Number';
            varStruct(numPara).value = componentNum(:);
            
            numPara = numPara + 1;
            varStruct(numPara).tag = 'p-value';
            varStruct(numPara).value = pValue(:);
            
            numPara = numPara + 1;
            varStruct(numPara).tag = 'T-value';
            varStruct(numPara).value = tValue(:);
            
            icatb_printToFile(outFileName, varStruct, titlePrint, 'column_wise', 'append', additionalInfo);
            
    end
    %% End for checking design criteria
    
    
    %% Compute mean and standard deviation
    moreInfo(1).str = sprintf('%s\n', '.............................................');
    moreInfo(length(moreInfo) + 1).str = sprintf('%s\n', '          MEAN AND STANDARD DEVIATION        ');
    moreInfo(length(moreInfo) + 1).str = sprintf('%s\n', '.............................................');
    
    % Loop over components
    for nComp = 1:length(componentNum)
        countN = 0;
        startInd = 1;
        endInd = 0;
        % Loop over selected groups
        for nGroup = 1:length(subinGroups)
            groupNum = subinGroups(nGroup).val;
            endInd = endInd + length(groupNum);
            if ~isempty(eq_regressors)
                % Handle the case where operation is done on multiple
                % regressors to treat as a single regressor
                tmpStr = deblank(groupNames(selectedGroups(nGroup), :));
                if nGroup == 1
                    md = zeros(length(subinGroups), 1);
                    sd = md;
                end
                md(nGroup) = mean(data(startInd:endInd, nComp));
                sd(nGroup) = std(data(startInd:endInd, nComp));
                if nGroup == 1
                    groupCondStr = tmpStr;
                else
                    groupCondStr = str2mat(groupCondStr, tmpStr);
                end
            else
                
                % Loop over conditions
                for nC = 1:length(selectedConditions)
                    countN = countN + 1;
                    if countN == 1
                        md = zeros(length(subinGroups)*length(selectedConditions), 1);
                        sd = md;
                    end
                    inds = find(condInd == nC);
                    inds = inds(startInd:endInd);
                    md(countN) = mean(data(inds, nComp));
                    sd(countN) = std(data(inds, nComp));
                    tmpStr = [deblank(groupNames(selectedGroups(nGroup), :)), ' ', deblank(regressorNames(selectedConditions(nC), :))];
                    if countN == 1
                        groupCondStr = tmpStr;
                    else
                        groupCondStr = str2mat(groupCondStr, tmpStr);
                    end
                end
                % End loop over conditions
            end
            startInd = endInd + 1;
        end
        % End loop over selected groups
        
        tmpStr = str2mat(['Component ', num2str(nComp), ':'], '', 'Name', '', groupCondStr, '');
        meanStr = str2mat('', '', 'Mean', '', num2str(md(:)), '');
        stdStr = str2mat('', '', 'Standard Deviation', '', num2str(sd(:)), '');
        
        if nComp == 1
            strToPrint = tmpStr;
            mStr = meanStr;
            sdStr = stdStr;
        else
            strToPrint = str2mat(strToPrint, tmpStr);
            mStr = str2mat(mStr, meanStr);
            sdStr = str2mat(sdStr, stdStr);
        end
        
        
    end
    % End loop over components
    
    
    %% Print mean and standard deviation
    numPara = 1;
    varStruct(numPara).tag = '';
    varStruct(numPara).value = strToPrint;
    
    numPara = numPara + 1;
    varStruct(numPara).tag = '';
    varStruct(numPara).value = mStr;
    
    numPara = numPara + 1;
    varStruct(numPara).tag = '';
    varStruct(numPara).value = sdStr;
    
    icatb_printToFile(outFileName, varStruct, '', 'column_wise', 'append', str2mat(moreInfo.str));
    
    disp(['Done ', designCriteria, ' on beta weights']);
    
    disp(['Results are saved in output file ', outFileName]);
    
    set(handles, 'pointer', 'arrow');
    
catch
    
    msg = lasterror;
    
    if (exist('fid', 'var'))
        try
            fclose(fid);
        catch
        end
    end
    
    set(handles, 'pointer', 'arrow');
    rethrow(msg);
    
end

function addConditions(hObject, event_data, handles, flagAdd)
% Add other conditions callback

% Figure data
figureData = get(handles, 'userdata');
groupData = figureData.groupData;

if ~exist('flagAdd', 'var')
    flagAdd = 'new';
end

if isempty(groupData)
    error('Enter groups before adding other conditions');
end

parameters = figureData.parameters;

% Other conditions listbox
otherCondListH = findobj(handles, 'tag', 'selOtherConditions');

otherCond = figureData.otherCond;

condListVal =  get(otherCondListH, 'value');

if isempty(otherCond)
    flagAdd = 'new';
end

if strcmpi(flagAdd, 'new')
    otherCondNum = length(otherCond) + 1;
    otherCondName = ['Condition ', num2str(otherCondNum)];
    otherCondVector = '';
else
    otherCondNum = condListVal;
    otherCondName = otherCond(condListVal).name;
    otherCondVector = otherCond(condListVal).vector;
end


numParameters = 1;
inputText(numParameters).promptString = ['Enter other condition ', num2str(otherCondNum)];
inputText(numParameters).uiType = 'edit';
inputText(numParameters).answerString = otherCondName;
inputText(numParameters).uiPos = [0.4, 0.05];
inputText(numParameters).dataType = 'string';
inputText(numParameters).tag = 'other_cond_name';
inputText(numParameters).enable = 'on';

numParameters = numParameters + 1;

inputText(numParameters).promptString = 'Enter vector ([group1 group2 group3])';
inputText(numParameters).uiType = 'edit';
inputText(numParameters).answerString = otherCondVector;
inputText(numParameters).uiPos = [0.4, 0.05];
inputText(numParameters).dataType = 'string';
inputText(numParameters).tag = 'other_cond_vector';
inputText(numParameters).enable = 'on';

% Open an input dialog box to get answers
answers = icatb_inputDialog('inputtext', inputText, 'Title', 'Enter condition name and vector');

if isempty(answers)
    return;
end

for nn = 1:length(answers)
    if isempty(answers{nn})
        error('Error:OtherCond', '%s field is empty', inputText(nn).promptString);
    end
end


% Other condition name and vector
other_cond_name = answers{1};
other_cond_vector = answers{2};

try
    str2num(other_cond_vector);
catch
    error('Error:OtherCond', 'Value entered for %s doesn''t generate a valid number', ...
        inputText(1).promptString);
end

% Set other condition name and vector
otherCond(otherCondNum).name = other_cond_name;
otherCond(otherCondNum).vector = other_cond_vector;

if ~isempty(otherCond)
    checkIndex = strmatch(other_cond_name, str2mat(otherCond.name), 'exact');
    if length(checkIndex) > 1
        error('Error:CondName', 'Condition name %s already exists', other_cond_name);
    end
end

figureData.otherCond = otherCond;

set(otherCondListH, 'string', str2mat(otherCond.name));

% Set userdata to the figure
set(handles, 'userdata', figureData);


function removeCondCallback(hObject, event_data, handles)
% Remove conditions callback


% Other conditions listbox
otherCondListH = findobj(handles, 'tag', 'selOtherConditions');

figureData = get(handles, 'userdata');
otherCond = figureData.otherCond;

otherCondVal = get(otherCondListH, 'value');

if ~isempty(otherCond)
    msgStr = 'Do you want to remove the selected conditions?';
    [removeConds] = icatb_questionDialog('title', 'Remove conditions', 'textbody', msgStr);
    if removeConds
        otherCond(otherCondVal) = [];
    else
        return;
    end
end

set(otherCondListH, 'value', 1);

if ~isempty(otherCond)
    set(otherCondListH, 'string', str2mat(otherCond.name));
else
    set(otherCondListH, 'string', '');
end

% Set figure Data
figureData.otherCond = otherCond;

set(handles, 'userdata', figureData);


function deleteHandles(tagName)
% Delete figures that have matching tag

allObj = findobj('tag', tagName);
for nn = 1:length(allObj)
    delete(allObj(nn));
end


function doubleClickGroups(hObject, event_data, handles)
% Detect double click of listbox

addGroupsH = findobj(handles, 'tag',  'addGroups');

if strcmpi(get(handles, 'SelectionType'), 'open')
    addGroupsCallback(addGroupsH, [], handles, 'reuse');
end


function doubleClickOtherCond(hObject, event_data, handles)
% Detect double click of listbox

addCondH = findobj(handles, 'tag',  'addOtherConditions');

if strcmpi(get(handles, 'SelectionType'), 'open')
    addConditions(addCondH, [], handles, 'reuse');
end


function addContrasts(hObject, event_data, handles, flagAdd)
% Add contrasts

% Figure data
figureData = get(handles, 'userdata');
groupData = figureData.groupData;

if ~exist('flagAdd', 'var')
    flagAdd = 'new';
end

if isempty(groupData)
    error('Enter groups before adding contrasts');
end

parameters = figureData.parameters;

%%%% Check Design Criteria %%%%%%%%
designCriteriaH = findobj(handles, 'tag', 'design_criteria');
designOptions = get(designCriteriaH, 'string');
designCriteriaVal = get(designCriteriaH, 'value');
designCriteria = designOptions{designCriteriaVal};
%%%% End for checking design criteria %%%%%

% Selected conditions
selConditionsH = findobj(handles, 'tag', 'selConditions');
regressorNames = get(selConditionsH, 'string');
selectedConditions = get(selConditionsH, 'value');

% Selected groups
selGroupsH = findobj(handles, 'tag', 'selGroups');
selectedGroups = get(selGroupsH, 'value');

numGroups = length(selectedGroups);
numCond = length(selectedConditions);

checkConStr = '';
if strcmpi(designCriteria, 'one way anova (groups)')
    checkConStr = ['of length (', num2str(numGroups), ') equal to number of selected groups '];
elseif strcmpi(designCriteria, 'one way anova (regressors)')
    checkConStr = ['of length (', num2str(numCond), ') equal to number of selected regressors'];
elseif strcmpi(designCriteria, 'two way anova (groups & regressors)')
    checkConStr = ['of length (', num2str(numGroups + numCond), ') equal to number of selected groups and regressors'];
end

% contrasts
contrastsH = findobj(handles, 'tag', 'selContrasts');

otherCond = figureData.otherCond;

% Contrasts data structure
contrasts = figureData.contrasts;

contrastListVal =  get(contrastsH, 'value');

if isempty(contrasts)
    flagAdd = 'new';
end

if strcmpi(flagAdd, 'new')
    contrastNum = length(contrasts) + 1;
    contrastName = ['Contrast ', num2str(contrastNum)];
    contrastVector = '';
else
    contrastNum = contrastListVal;
    contrastName = contrasts(contrastListVal).name;
    contrastVector = contrasts(contrastListVal).vector;
end

numParameters = 1;

inputText(numParameters).promptString = ['Enter contrast name ', num2str(contrastNum)];
inputText(numParameters).uiType = 'edit';
inputText(numParameters).answerString = contrastName;
inputText(numParameters).uiPos = [0.4, 0.05];
inputText(numParameters).dataType = 'string';
inputText(numParameters).tag = 'contrast_name';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = 1;

numParameters = numParameters + 1;
inputText(numParameters).promptString = ['Enter contrast vector ', checkConStr];
inputText(numParameters).uiType = 'edit';
inputText(numParameters).answerString = contrastVector;
inputText(numParameters).uiPos = [0.4, 0.05];
inputText(numParameters).dataType = 'string';
inputText(numParameters).tag = 'contrast_vector';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = 1;

% Open an input dialog box to get answers
answers = icatb_inputDialog('inputtext', inputText, 'Title', 'Enter contrast name and vector');

if isempty(answers)
    return;
end

for nn = 1:length(answers)
    if isempty(answers{nn})
        error('Error:Contrasts', '%s field is empty', inputText(nn).promptString);
    end
end


% Contrast name and vector
contrastName = answers{1};
contrastVector = answers{2};

try
    checkSum = (sum(str2num(contrastVector)) == 0);
    if ~checkSum
        error('Sum of contrasts must be zero');
    end
catch
    icatb_displayErrorMsg;
end

% Set other condition name and vector
contrasts(contrastNum).name = contrastName;
contrasts(contrastNum).vector = contrastVector;

if ~isempty(contrasts)
    checkIndex = strmatch(contrastName, str2mat(contrasts.name), 'exact');
    if length(checkIndex) > 1
        error('Error:Contrasts', 'Contrast name %s already exists', contrastName);
    end
end

figureData.contrasts = contrasts;

set(contrastsH, 'string', str2mat(contrasts.name));

% Set userdata to the figure
set(handles, 'userdata', figureData);


function removeContrasts(hObject, event_data, handles)
% Remove contrasts callback

contrastListH = findobj(handles, 'tag', 'selContrasts');

figureData = get(handles, 'userdata');
contrasts = figureData.contrasts;

contrastVal = get(contrastListH, 'value');

if ~isempty(contrasts)
    msgStr = 'Do you want to remove the selected contrasts?';
    [removeContrasts] = icatb_questionDialog('title', 'Remove contrasts', 'textbody', msgStr);
    if removeContrasts
        contrasts(contrastVal) = [];
    else
        return;
    end
end

set(contrastListH, 'value', 1);

if ~isempty(contrasts)
    set(contrastListH, 'string', str2mat(contrasts.name));
else
    set(contrastListH, 'string', '');
end

% Set figure Data
figureData.contrasts = contrasts;

set(handles, 'userdata', figureData);



function doubleClickContrast(hObject, event_data, handles)
% Detect double click of listbox

addContrastsH = findobj(handles, 'tag',  'addContrasts');

if strcmpi(get(handles, 'SelectionType'), 'open')
    addContrasts(addContrastsH, [], handles, 'reuse');
end


function popupCallback(hObject, event_data, handles)
% Popup callback

popupOptions = get(hObject, 'string');

popupVal = get(hObject, 'value');

selectedStr = popupOptions{popupVal};

RegressorListH = findobj(handles, 'tag', 'selConditions');
groupsListH = findobj(handles, 'tag', 'selGroups');
selContrastsH = findobj(handles, 'tag', 'selContrasts');
addContrastsH = findobj(handles, 'tag', 'addContrasts');
removeContrastsH = findobj(handles, 'tag', 'removeContrasts');

set(groupsListH, 'min', 0);
set(groupsListH, 'max', 2);

groupVal = get(groupsListH, 'value');
regressVal = get(RegressorListH, 'value');

if strcmpi(selectedStr, 'one sample t-test') || strcmpi(selectedStr, 'two sample t-test')
    %% One sample t-test and two sample t-test
    
    % Make single selection listbox
    set(RegressorListH, 'value', regressVal(1));
    set(RegressorListH, 'min', 0);
    set(RegressorListH, 'max', 1);
    
    if strcmpi(selectedStr, 'one sample t-test')
        set(groupsListH, 'value', groupVal(1));
        set(groupsListH, 'min', 0);
        set(groupsListH, 'max', 1);
    end
    
    set(selContrastsH, 'enable', 'off');
    set(addContrastsH, 'enable', 'off');
    set(removeContrastsH, 'enable', 'off');
    
elseif strcmpi(selectedStr, 'multiple regression')
    %% Regression
    set(selContrastsH, 'enable', 'off');
    set(addContrastsH, 'enable', 'off');
    set(removeContrastsH, 'enable', 'off');
    
    % Multiple selection listbox
    set(RegressorListH, 'min', 0);
    set(RegressorListH, 'max', 2);
    
    addMultiRegressors(handles);
    
else
    %% Anova
    set(selContrastsH, 'enable', 'on');
    set(addContrastsH, 'enable', 'on');
    set(removeContrastsH, 'enable', 'on');
    
    % Multiple selection listbox
    set(RegressorListH, 'min', 0);
    set(RegressorListH, 'max', 2);
    
end

function [subinGroups, groupID] = get_groupID(groupData, numOfSess)

subinGroups = repmat(struct('val', []), length(groupData), 1);
% Loop over groups
for nGroup = 1:length(groupData)
    
    subListVal = groupData(nGroup).subVal;
    sessListVal = groupData(nGroup).sessVal;
    subjectNumbers = zeros(1, length(subListVal)*length(sessListVal));
    
    countSet = 0;
    
    % Loop over number of subjects
    for nSub = subListVal
        % Loop over number of sessions
        for nSess = sessListVal
            countSet = countSet + 1;
            subjectNumbers(countSet) = (nSub - 1)*numOfSess + nSess;
        end
        % End loop over sessions
    end
    % End loop over subjects
    
    % Subject numbers
    tmpGroup = nGroup*ones(length(subjectNumbers), 1);
    if nGroup == 1
        groupID = tmpGroup;
    else
        groupID = [groupID; tmpGroup];
    end
    
    subinGroups(nGroup).val = subjectNumbers;
    
    clear tmpGroup;
    clear subjectNumbers;
    
end
% End loop over groups


function [data, condInd] = get_regress_data(regressInfo, subinGroups, selectedConditions, componentNum)
% Get data

allSubjects = [subinGroups.val];

% Loop over selected conditions
for nCond = 1:length(selectedConditions)
    tmp = zeros(length(allSubjects), length(componentNum));
    for nSub = 1:length(allSubjects)
        tmp(nSub, :) = (regressInfo(allSubjects(nSub)).cond(selectedConditions(nCond)).values(:))';
    end
    tmpCondInd = nCond*ones(size(tmp, 1), 1);
    if nCond == 1
        data = tmp;
        condInd = tmpCondInd;
    else
        data = [data; tmp];
        condInd = [condInd; tmpCondInd];
    end
    clear tmp;
end
% End loop over selected conditions



function  [selectedGroups, designVal] = errChk(designCriteria, selectedGroups, selectedConditions)
%% Initial error check on groups and conditions
%


switch (lower(designCriteria))
    case 'one sample t-test'
        %% One sample t-test
        
        selectedGroups = selectedGroups(1);
        designVal = 1;
        
    case 'two sample t-test'
        %% Two sample t-test
        
        if (length(selectedGroups) == 1)
            error('Need two groups to do two sample t-test');
        end
        
        if (length(selectedGroups) > 2)
            selectedGroups = selectedGroups(1:2);
        end
        
        designVal = 2;
        
    case 'one way anova (groups)'
        %% One way anova for groups
        
        if (length(selectedGroups) < 2)
            error('Need atleast two groups to do one way anova (groups)');
        end
        
        designVal = 3;
        
    case 'one way anova (regressors)'
        %% One way anova for regressors
        
        if (length(selectedConditions) < 2)
            error('Need atleast two regressors to do one way anova (regressors)');
        end
        
        designVal = 4;
        
    case  'two way anova (groups & regressors)'
        %% Two way anova for groups and regressors
        
        if (length(selectedGroups) < 2)
            error('Need atleast two groups to do two way anova (groups, regressors)');
        end
        if (length(selectedConditions) < 2)
            error('Need atleast two regressors to do two way anova (groups, regressors)');
        end
        
        designVal = 5;
        
    case 'multiple regression'
        %% Multiple regression
        
        designVal = 6;
        
    case 'paired t-test'
        %% Paired t-test
        
        designVal = 7;
        
        if ((length(selectedGroups) ~= 2) && (length(selectedConditions) ~= 2))
            error('Need two conditions to do paired t-test');
        end
        
    otherwise
        
        error('Unknown design criteria specified');
        
end

function [eq_regressors, eqStr] = formEqRegressors(designValue, selectedConditions, userInputFile, regressorNames)
%% Form equations for regressors

eq_regressors = [];
eqStr = '';

%% Check if it is one way anova (groups)
if ((designValue == 3) && (length(selectedConditions) > 1))
    
    if isempty(userInputFile)
        % Open an input dialog box for specifying equation for regressors if
        % one way anova (Groups) is selected
        % open input dialog box
        prompt = {['Enter equation for regressors of length (', num2str(length(selectedConditions)), ') or leave empty.']};
        dlg_title = 'Equation for regressors';
        num_lines = 1;
        def = {''};
        % equation for regressors
        answer = icatb_inputdlg2(prompt, dlg_title, num_lines, def);
    else
        % Read equation for regressors from input file
        keywd = 'eq_regressors';
        try
            inputData = icatb_read_variables(userInputFile, keywd, 'vector', 'numeric');
            answer{1} = num2str(getfield(inputData, keywd));
            clear inputData;
        catch
            answer{1} = '';
        end
    end
    
    if isempty(answer)
        answer = {''};
    end
    
    try
        eq_regressors = str2num(answer{1});
    catch
        error('Error:EqRegress', 'Please check the equation %s for regressors', answer{1});
    end
    
    if ~isempty(eq_regressors)
        if length(eq_regressors) ~= length(selectedConditions)
            error('Error:EqRegress', 'Length of equation for regressors (%s) doesn''t match \nlength of selected conditions (%s)', ...
                num2str(length(eq_regressors)), num2str(length(selectedConditions)));
        end
        
        fprintf('\n');
        for nC = 1:length(eq_regressors)
            
            tmpStr = [num2str(eq_regressors(nC)),  '*', deblank(regressorNames(selectedConditions(nC), :))];
            if nC == 1
                eqStr = tmpStr;
            else
                eqStr = [eqStr, ' + ', tmpStr];
            end
        end
        
        disp(['Equation ', eqStr, ' will be applied on data for doing one way anova (groups)']);
        
        fprintf('\n');
    end
    
end
%% End for checking one way anova (groups)



function addMultiRegressors(handles)
%% Add regressors for multiple regression criteria
%

%% Get information from figure data
figureData = get(handles, 'userdata');
parameters = figureData.parameters;
numOfSub = parameters.numOfSub;
numOfSess = parameters.numOfSess;
userInputFile = figureData.userInputFile;

if isempty(userInputFile)
    
    %% Select regressors
    selectRegress = 1;
    if (isfield(parameters, 'multiple_regression'))
        msgStr = 'Do you want to re-select the regressor files for Multiple Regression?';
        selectRegress = icatb_questionDialog('title', 'Select regressor files', 'textbody', msgStr);
    end
    
    if selectRegress
        %% Load GUI for getting the regressors
        regressorFiles = icatb_selectEntry('typeEntity', 'file', 'typeSelection', 'multiple', 'filter', '*.txt;*.dat;*.mat', 'title', ...
            'Select regressor file/files (Age, scores, etc) for using Multiple Regression ...');
    else
        return;
    end
    
else
    %% Use regressor files from input file
    regressorFiles = parameters.multiple_regression.multi_regress_files;
    
end

drawnow;

if (isempty(regressorFiles))
    error('Regressor files are not selected for using Multiple Regression');
end

regressorFiles = str2mat(regressorFiles);

%% Loop over files
for nF = 1:size(regressorFiles, 1)
    
    currentFile = deblank(regressorFiles(nF, :));
    
    [pathstr, fName] = fileparts(currentFile);
    
    data = load(currentFile);
    
    if (isstruct(data))
        %% Load as MAT file
        fN = fieldnames(data);
        for nField = 1:length(fN)
            temp = getfield(data, fN{nField});
            if isnumeric(temp)
                clear data;
                data = temp;
                break;
            end
        end
    end
    
    %% Error check on data
    if ~(isnumeric(data))
        error('Error:Data', 'Data is not numeric for file %s\n', currentFile);
    end
    
    size_data = size(data);
    if length(size_data) > 2
        error('Error:Data', 'Data must be in 2 dimensions for file %s\n', currentFile);
    end
    
    if (size_data(1) == 1)
        data = data';
    end
    
    clear size_data;
    
    if (size(data, 1) ~= numOfSub)
        error('Error:Data', 'No. of rows of data (%d) doesn''t match the number of subjects (%d) \n', size(data, 1), numOfSub);
    end
    
    if (size(data, 2) > 1)
        tmpNames = cellstr(strcat(fName, '_', num2str((1:size(data, 2))')));
    else
        tmpNames = {fName};
    end
    
    %% Replicate data over sessions
    if (numOfSess > 1)
        newData = zeros(size(data, 1)*numOfSess, size(data, 2));
        startInd = 1; endInd = 0;
        % Loop over subjects
        for nD = 1:size(data, 1)
            endInd = endInd + numOfSess;
            temp = data(nD, :);
            temp = repmat(temp, numOfSess, 1);
            newData(startInd:endInd, :) = temp;
            startInd = endInd + 1;
        end
        % End loop over subjects
        clear data;
        data = newData;
        clear newData;
    end
    %% End for replicating data over sessions
    
    if (nF == 1)
        regressData = data;
        regressNames = tmpNames;
    else
        regressData = [regressData, data];
        regressNames = [regressNames; tmpNames];
    end
    
end
%% End loop over files

%% Set figure data
parameters.multiple_regression.regressorFiles = regressorFiles;
parameters.multiple_regression.data = regressData;
parameters.multiple_regression.names = regressNames;

figureData.parameters = parameters;
set(handles, 'userdata', figureData);