function icatb_spatial_chronn_stats_gui(param_file)
%% Compute statistics on the transition matrices and subject states
%

global UI_FS;

if (~exist('param_file', 'var'))
    param_file = icatb_selectEntry('title', 'Select spatial chronnectome Parameter File', 'typeEntity', 'file', 'typeSelection', ...
        'single', 'filter', '*schronn.mat');
    drawnow;
end

load(param_file);

covNames = '';

try
    statsInfo = schronnInfo.statsInfo;
catch
    statsInfo.numOfSub = schronnInfo.userInput.numOfSub;
end

if (~isfield(statsInfo, 'cov'))
    statsInfo.cov = [];
end

try
    covNames = (cellstr(char(statsInfo.cov.name)));
catch
end


figureTag = 'stats_schronn';

% Setup figure for GUI
InputHandle = icatb_getGraphics('Setup Stats', 'normal', figureTag, 'on');
set(InputHandle, 'menubar', 'none');
set(InputHandle, 'userdata', statsInfo);


controlWidth = 0.4;
promptHeight = 0.05;
promptWidth = controlWidth;
listboxHeight = 0.2; listboxWidth = controlWidth;
xOffset = 0.02; yOffset = promptHeight; yPos = 0.98;
okWidth = 0.12; okHeight = promptHeight;

dropDownWidth = 0.4;

%% Design dropdown box
promptPos = [0.25 - 0.5*dropDownWidth, yPos - 0.5*yOffset, 0.4, promptHeight];
yPos = yPos - promptPos(4) - yOffset;

%% Features text and listbox
promptPos = [0.5 - 0.5*controlWidth, yPos - 0.5*yOffset, promptWidth, promptHeight];

textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Add Covariates', 'tag', ...
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

promptPos(1) = 0.05;

promptPos(2) = promptPos(2) - 1.5*yOffset;

textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', ...
    'Select p-threshold', 'tag', 'prompt_p_threshold', 'fontsize', UI_FS - 1);
icatb_wrapStaticText(textH);
promptPos = get(textH, 'position');

promptPos(1) = promptPos(1) + promptPos(3) + xOffset;
promptPos(3) = 0.3;
icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'edit', 'position', promptPos, 'string', '0.05', 'tag', ...
    'p_threshold', 'fontsize', UI_FS - 1);


promptPos(1) = 0.05;

promptPos(2) = promptPos(2) - 1.5*yOffset - promptPos(4);

textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', ...
    'Select threshold criteria', 'tag', 'prompt_threshold_criteria', 'fontsize', UI_FS - 1);
icatb_wrapStaticText(textH);
promptPos = get(textH, 'position');

promptPos(1) = promptPos(1) + listboxPos(3) + xOffset;
promptPos(3) = 0.3;
icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'popup', 'position', promptPos, 'string', {'fdr', 'MAFDR', 'none'}, ...
    'tag', 'threshdesc', 'fontsize', UI_FS - 1);

%% Add cancel, save and run buttons
okPos = [0.5 - 0.5*okWidth, promptPos(2) - yOffset - 0.5*okHeight, okWidth, okHeight];
okPos(2) = okPos(2) - 0.5*okPos(4);
icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', okPos, 'string', 'Done', 'tag', 'stats_button', 'fontsize',...
    UI_FS - 1, 'callback', {@okCallback, InputHandle});

appName = 'statsSchronnAppData';

waitfor(InputHandle);


if (isappdata(0, appName))
    
    statsInfo = getappdata(0, appName);
    rmappdata(0, appName);
    schronnInfo.postprocess.statsInfo = statsInfo;
    drawnow;
    icatb_spatial_chronn_stats(schronnInfo);
    
else
    
    error('StatsWindow was quit');
    
end



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

statsInfo = get(figH, 'userdata');

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
        covName = statsInfo.cov(val).name;
        covVals = statsInfo.cov(val).value;
        if (isnumeric(covVals))
            covVals = covVals(:);
            covVals = num2str(covVals);
        end
        covVals = cellstr(covVals);
        covVals = covVals(:)';
        transformationName = statsInfo.cov(val).transformation;
        cov_type = statsInfo.cov(val).type;
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

statsInfo = get(figH, 'userdata');

txtFile = icatb_selectEntry('title', 'Select covariate file' , 'filter', '*.txt;*.asc', 'typeEntity', 'file', 'typeSelection', 'single');
drawnow;
covTypeH = findobj(handles, 'tag', 'cov_type');
opts = cellstr(get(covTypeH, 'string'));
covVal = get(covTypeH, 'value');

try
    val = icatb_mancovan_load_covariates(txtFile, opts{covVal}, statsInfo.numOfSub);
    covValueH = findobj(handles, 'tag', 'cov_value');
    set(covValueH, 'string', val);
catch
    icatb_errorDialog(lasterr, 'Covariate Selection');
end


function setCovCallback(hObject, event_data, covFigH, handles)
%% Set covariate name, value and type

statsInfo = get(handles, 'userdata');

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
    
    if (length(statsInfo.cov) > 0)
        chk = strmatch(lower(cov_name), lower(cellstr(char(statsInfo.cov.name))), 'exact');
        if (~isempty(chk))
            ind = chk;
        end
    end
    
    if (~exist('ind', 'var'))
        ind = length(statsInfo.cov) + 1;
    end
    
    statsInfo.cov(ind).name = cov_name;
    statsInfo.cov(ind).value = deblank(cov_value(:)');
    statsInfo.cov(ind).transformation = lower(cov_transformation);
    statsInfo.cov(ind).type = covType;
    
    set(handles, 'userdata', statsInfo);
    
    covListH = findobj(handles, 'tag', 'cov');
    set(covListH, 'string', cellstr(char(statsInfo.cov.name)));
    delete(covFigH);
    
catch
    icatb_errorDialog(lasterr, 'Covariate Selection');
end


function removeCov(hObject, event_data, figH)
%% Remove Covariate
%

statsInfo = get(figH, 'userdata');
listH = findobj(figH, 'tag', 'cov');
val = get(listH, 'value');

if (~isempty(val))
    check = icatb_questionDialog('title', 'Remove Covariate', 'textbody', 'Do you want to remove the covariate from the list?');
    if (~check)
        return;
    end
end

try
    strs = cellstr(char(statsInfo.cov.name));
    statsInfo.cov(val) = [];
    strs(val) = [];
    set(listH, 'value', 1);
    set(listH, 'string', strs);
    set(figH, 'userdata', statsInfo);
catch
end

function covValueCallback(hObject, event_data, figH)
%% Covariate value callback
%

statsInfo = get(figH, 'userdata');
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

if (length(val2) ~= statsInfo.numOfSub)
    error(['Covariate vector length must equal the no. of subjects (', num2str(statsInfo.numOfSub), ')']);
end

% if (isnumeric(val2))
%     val2 = val2(:);
%     val2 = num2str(val2);
% end

val2 = strtrim(cellstr(val2));

set(hObject, 'string', val2);


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


function okCallback(hObject, ed, handles)

statsInfo = get(handles, 'userdata');

p_threshold = str2num(get(findobj(handles, 'tag', 'p_threshold'), 'string'));
threshval = get(findobj(handles, 'tag', 'threshdesc'), 'value');
threshstr = get(findobj(handles, 'tag', 'threshdesc'), 'string');
threshdesc = lower(deblank(threshstr{threshval}));

statsInfo.threshdesc = threshdesc;
statsInfo.p_threshold = p_threshold;

univariate_tests = {};

try
    all_covariates = cellstr(char(statsInfo.cov.name));
    if (length(all_covariates) >= 1)
        if (length(all_covariates) == 1)
            univariate_tests = {all_covariates{1}, {}};
        else
            tests = icatb_univariate_cov_sel(all_covariates, 'Do you want to specify nuisance covariates for each covariate?');
            if (isempty(tests))
                %tests = [];
                for nT = 1:length(all_covariates)
                    tests(end+1).name = all_covariates{nT};
                end
            end
            univariate_tests = cell(length(tests), 2);
            for nS = 1:length(tests)
                univariate_tests{nS, 1} = tests(nS).name;
                if (isfield(tests(nS), 'str') && ~isempty(tests(nS).str))
                    univariate_tests{nS, 2} = tests(nS).str(:)';
                end
            end
            
        end
    end
    
catch
    
end

statsInfo.univariate_tests = univariate_tests;

appName = 'statsSchronnAppData';
setappdata(0, appName, statsInfo);

%set(handles, 'userdata', 'statsInfo');

delete(handles);