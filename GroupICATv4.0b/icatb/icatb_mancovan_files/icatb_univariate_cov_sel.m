function regressOut = icatb_univariate_cov_sel(covariates)
%% Select covariates of interest and include covariates in the reduced model.

%% Load defaults
icatb_defaults;
global UI_FS;

regressOut = [];

%% Open figure
InputHandle = icatb_getGraphics('Univariate Tests', 'normal', 'sel_covariates_univariate', 'on');
set(InputHandle, 'menubar', 'none');

controlWidth = 0.25;
promptHeight = 0.05;
promptWidth = controlWidth;
listboxHeight = controlWidth; listboxWidth = controlWidth;
xOffset = 0.02; yOffset = promptHeight; yPos = 0.9;
okWidth = 0.12; okHeight = promptHeight;

%% Skip Multivariate?
promptPos = [0.1, yPos - 0.5*yOffset, 0.6, promptHeight];

textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Do You Want To Skip Multivariate tests?', 'tag', ...
    'prompt_design', 'fontsize', UI_FS - 1);
icatb_wrapStaticText(textH);
promptPos = get(textH, 'position');

promptPos(1) = promptPos(1) + promptPos(3) + xOffset;
promptPos(3) = 0.16;
multH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'popup', 'position', promptPos, 'string', {'No', 'Yes'}, 'tag', ...
    'chk_multivariate', 'fontsize', UI_FS - 1, 'value', 1, 'callback', {@multCallback, InputHandle});

yOffset = 0.12;

yPos = yPos - promptPos(4) - yOffset;

%% All covariates
promptPos = [xOffset, yPos - 0.5*yOffset, promptWidth, promptHeight];

textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'All covariates', 'tag', ...
    'prompt_all_covariates', 'fontsize', UI_FS - 1);
icatb_wrapStaticText(textH);
listboxXOrigin = promptPos(1);
listboxYOrigin = promptPos(2) - yOffset - listboxHeight;
listboxPos = [listboxXOrigin, listboxYOrigin, listboxWidth, listboxHeight];
icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'listbox', 'position', listboxPos, 'string', covariates, 'tag', ...
    'cov', 'fontsize', UI_FS - 1, 'min', 0, 'max', 1, 'value', 1);

addButtonPos = [listboxPos(1) + listboxPos(3) + xOffset, listboxPos(2) + 0.5*listboxPos(4) + 0.5*promptHeight, promptHeight + 0.01, promptHeight - 0.01];
icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', addButtonPos, 'string', '+', 'tag', 'add_cov_button', 'fontsize',...
    UI_FS - 1, 'callback', {@addCov, InputHandle});
removeButtonPos = [listboxPos(1) + listboxPos(3) + xOffset, listboxPos(2) + 0.5*listboxPos(4) - 0.5*promptHeight, promptHeight + 0.01, promptHeight - 0.01];
icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', removeButtonPos, 'string', '-', 'tag', 'remove_cov_button', 'fontsize',...
    UI_FS - 1, 'callback', {@removeCov, InputHandle});

%% Selected covariates
promptPos(1) = addButtonPos(1) + addButtonPos(3) + 2*xOffset;
textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Selected covariates', 'tag', ...
    'prompt_all_covariates', 'fontsize', UI_FS - 1);
icatb_wrapStaticText(textH);
listboxXOrigin = promptPos(1);
listboxYOrigin = promptPos(2) - yOffset - listboxHeight;
listboxPos = [listboxXOrigin, listboxYOrigin, listboxWidth, listboxHeight];
icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'listbox', 'position', listboxPos, 'string', '', 'tag', ...
    'sel_cov', 'fontsize', UI_FS - 1, 'min', 0, 'max', 1, 'value', 1, 'callback', {@selectedCov, InputHandle});

%% Reduced model
promptPos(1) = promptPos(1) + promptPos(3) + 2*xOffset;

textH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Reduced Model (Nuisance covariates)', 'tag', ...
    'prompt_reduced_covariates', 'fontsize', UI_FS - 1);
icatb_wrapStaticText(textH);
listboxXOrigin = promptPos(1);
listboxYOrigin = promptPos(2) - yOffset - listboxHeight;
listboxPos = [listboxXOrigin, listboxYOrigin, listboxWidth, listboxHeight];
icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'listbox', 'position', listboxPos, 'string', '', 'tag', ...
    'reduced_cov', 'fontsize', UI_FS - 1, 'min', 0, 'max', 2, 'value', [],  'callback', {@reducedCov, InputHandle});


%% Add ok button

okPos = [0.6 - 0.5*okWidth, 0.2, okWidth, okHeight];
okPos(2) = okPos(2) - 0.5*okPos(4);
icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', okPos, 'string', 'Done', 'tag', 'Select', 'fontsize',...
    UI_FS - 1, 'callback', {@okCovCallback, InputHandle});

helpPos = [0.4 - 0.5*okWidth, 0.2, okWidth, okHeight];
helpPos(2) = helpPos(2) - 0.5*helpPos(4);
icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', helpPos, 'string', '?', 'tag', 'help', 'fontsize',...
    UI_FS - 1, 'callback', {@helpCallback, InputHandle}, 'foregroundcolor', [1, 1, 0]);

multCallback(findobj(InputHandle, 'tag', 'chk_multivariate'), [], InputHandle);

waitfor(InputHandle);

appName = 'UnivariateAppdata';
if (isappdata(0, appName))
    regressOut = getappdata(0, appName);
    rmappdata(0, appName);
else
    error('Figure window was quit');
end


function multCallback(hObject, event_data, handles)
%% Skip multivariate option?

addH = findobj(handles, 'tag', 'add_cov_button');
removeH = findobj(handles, 'tag', 'remove_cov_button');
covH = findobj(handles, 'tag', 'cov');
selCovH = findobj(handles, 'tag', 'sel_cov');
reducedH = findobj(handles, 'tag', 'reduced_cov');

val = get(hObject, 'value');
strs = get(hObject, 'string');

set(addH, 'enable', 'on');
set(removeH, 'enable', 'on');
set(covH, 'enable', 'on');
set(reducedH, 'enable', 'on');
set(selCovH, 'enable', 'on');

if (strcmpi(strs{val}, 'no'))
    set(addH, 'enable', 'inactive');
    set(removeH, 'enable', 'inactive');
    set(covH, 'enable', 'inactive');
    set(reducedH, 'enable', 'inactive');
    set(selCovH, 'enable', 'inactive');
end


function addCov(hObject, event_data, handles)
%% Add covariates callback
%

hd = get(handles, 'userdata');
covH = findobj(handles, 'tag', 'cov');
selCovH = findobj(handles, 'tag', 'sel_cov');
allCovariates = get(covH, 'string');
val = get(covH, 'value');
hd(end + 1).name = allCovariates{val};

set(selCovH, 'value', 1);
set(selCovH, 'string', cellstr(char(hd.name)));
set(handles, 'userdata', hd);

selectedCov(selCovH, [], handles);


function removeCov(hObject, event_data, handles)
%% Remove covariates callback
%

hd = get(handles, 'userdata');
covH = findobj(handles, 'tag', 'cov');
selCovH = findobj(handles, 'tag', 'sel_cov');
reducedH = findobj(handles, 'tag', 'reduced_cov');
val = get(selCovH, 'value');
if (~isempty(val))
    hd(val) = [];
    set(selCovH, 'value', 1);
    set(selCovH, 'string', cellstr(char(hd.name)));
    set(reducedH, 'value', []);
    set(reducedH, 'string', '');
    set(handles, 'userdata', hd);
end

function selectedCov(hObject, event_data, handles)
%% Selected covariates callback
%

hd = get(handles, 'userdata');
covH = findobj(handles, 'tag', 'cov');
selCovH = findobj(handles, 'tag', 'sel_cov');
reducedH = findobj(handles, 'tag', 'reduced_cov');
val = get(selCovH, 'value');

allNames = get(covH, 'string');
if (~isempty(val))
    [dd, ia, ib] = intersect(allNames, hd(val).name);
    allNames(ia) = [];
    
    selVals = [];
    try
        selVals = hd(val).val;
    catch
    end
    
    set(reducedH, 'value', selVals);
    set(reducedH, 'string', allNames);
    set(handles, 'userdata', hd);
end

function reducedCov(hObject, event_data, handles)
%% Reduced Model Callback

hd = get(handles, 'userdata');
covH = findobj(handles, 'tag', 'cov');
selCovH = findobj(handles, 'tag', 'sel_cov');
reducedH = findobj(handles, 'tag', 'reduced_cov');
val = get(selCovH, 'value');

if (~isempty(val))
    selStr = cellstr(get(reducedH, 'string'));
    selVal = get(reducedH, 'value');
    if (~isempty(selVal))
        hd(val).val = selVal;
        hd(val).str = selStr(selVal);
    end
    
    set(handles, 'userdata', hd);
end

function okCovCallback(hObject, event_data, handles)
%% OK button callback
%

hd = get(handles, 'userdata');


univH = findobj(handles, 'tag', 'chk_multivariate');
multStr = get(univH, 'string');
multVal = get(univH, 'value');

answer = lower(multStr{multVal});

univInfo = [];
if (strcmpi(answer, 'yes'))
    if (isempty(hd) || isempty(hd(1).name))
        error('No covariates selected');
    end
    univInfo = hd;
end

setappdata(0, 'UnivariateAppdata', univInfo);

delete(handles);


function helpCallback(hObject, event_data, handles)
%% Help Callback

msgStr = ['When you select "Yes" for "Do you want to skip multivariate tests?", you have the option to skip multivariate tests and proceed directly to univariate tests. Select covariates using the + button and this will appear in the selected covariates list. For each selected covariate, ', ...
    'you could remove variance associated with the other regressors by ctrl key and clicking on the regressors in the reduced model (nuisance covariates). If you accidentally selected covariate and wish to remove from the selected covariates listbox, use - button.'];
H = icatb_dialogBox('title', 'Covariates selection', 'textbody', msgStr, 'textType', 'large');
waitfor(H);
