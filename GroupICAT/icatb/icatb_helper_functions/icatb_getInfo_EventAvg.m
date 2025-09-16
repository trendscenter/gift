function appData = icatb_getInfo_EventAvg(listString, refAppData, numComp, handleVisibility)
% list dilaog box to get the subjects, sessions, components, regressor
% onsets
%
% Inputs:
% 1. numOfSub - Number of subjects
% 2. numOfSess - Number of sessions
% 3. numComp - Number of components
% 4. refAppData - reference data containing onsets, time course, names
% 5. spmMatFlag - spm mat flag is used to distinguish the regressors
% between subjects and sessions

% procedure
% first draw four listboxes three on top and one at the bottom.
% Regressor string will depend on how the subjects and sessions are
% selected

icatb_defaults;

if ~exist('handleVisibility', 'var')
    handleVisibility = 'on';
end

% Screen Color Defaults
global BG_COLOR;
global BG2_COLOR;
global BUTTON_COLOR;
global BUTTON_FONT_COLOR;
global FONT_COLOR;
global AXES_COLOR;
global FONT_COLOR;
global HELP_FONT_COLOR;

% font defaults
global UI_FONTNAME;
global UI_FONTUNITS;
global UI_FS;

for ii = 1:numComp
    compString(ii).name = num2str(ii);
end

compString = str2mat(compString.name);

titleFig = 'Select regressor and component for event average';

% Figure
[InputHandle] = icatb_getGraphics(titleFig, 'normal', titleFig, handleVisibility);

% user data for the figure
handles_data = struct('refAppData', refAppData, 'numComp', numComp);

set(InputHandle, 'menubar', 'none', 'userdata', handles_data);
% offsets
xOffset = 0.04; yOffset = 0.04;

% text width and height
textWidth = 0.35; textHeight = 0.05;

% button height and width
buttonHeight = 0.05; buttonWidth = 0.2;

% default list width and height
listHeight = (1 - 6*yOffset - buttonHeight - 2*textHeight) / 2;

%%%%%%%%%%%%%%%%%%%%%%%%%% Plot UI Controls %%%%%%%%%%%%%%%%%%%%%%%%

subjectTextPos = [xOffset (1 - yOffset - 0.5*textHeight) textWidth textHeight];
% Plot subject text
subjectTextH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', subjectTextPos, ...
    'BackgroundColor', BG2_COLOR, 'ForegroundColor', FONT_COLOR, 'fontunits', UI_FONTUNITS, 'fontname', ...
    UI_FONTNAME, 'FontSize', UI_FS, 'horizontalalignment', 'center', 'string', 'Datasets', 'userdata', ...
    [], 'tag', 'subject_text');

subjectPos = [subjectTextPos(1) subjectTextPos(2) - listHeight - yOffset textWidth listHeight];
% plot subject listbox
subjectH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'listbox', 'position', subjectPos, ...
    'horizontalalignment', 'left', 'min', 0,  'max', 1, 'string', [], 'tag', 'subjectListbox', 'value', 1, 'string', listString, ...
    'callback', {@subjectCallback, InputHandle});


spaceBetListBox = (1 - 2*xOffset - 2*textWidth) / 2;

% plot component listbox and text
compTextPos = [subjectTextPos(1) + subjectTextPos(3) + spaceBetListBox subjectTextPos(2) textWidth textHeight];
% Plot comp text
compTextH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', compTextPos, ...
    'horizontalalignment', 'center', 'string', 'Component', 'userdata', [], 'tag', 'component_text');

compPos = [compTextPos(1) compTextPos(2) - listHeight - yOffset textWidth listHeight];

% plot component listbox
componentH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'listbox', 'position', compPos, ...
    'horizontalalignment', 'left', 'string', [], 'tag', 'componentListbox', 'min', 0, 'max', 2, 'value', 1, 'string', compString);

textWidth = 0.4;

% plot text and listbox for regressors
regressTextPos = [xOffset, subjectPos(2) - yOffset - textHeight, textWidth textHeight];

% Plot regress text
regressTextH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', regressTextPos, ...
    'horizontalalignment', 'center', 'string', 'Regressor', 'userdata', ...
    [], 'tag', 'regress_text');

regressPos = [xOffset regressTextPos(2) - listHeight - yOffset textWidth listHeight];

% plot regress listbox
regressH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'listbox', 'position', regressPos, ...
    'horizontalalignment', 'left', 'min', 0,  'max', 2, 'string', ...
    [], 'tag', 'regressorListbox', 'value', 1, 'callback', {@regressorCallback, InputHandle});


buttonHeight = 0.05; buttonWidth = 0.12;

% plot button and help
okPos = [regressPos(1) + regressPos(3) + 2*xOffset regressPos(2) + 0.5*buttonHeight buttonWidth buttonHeight];

helpPos = [okPos(1) okPos(2) + okPos(4) + yOffset buttonWidth buttonHeight];

ResetPos = [helpPos(1) helpPos(2) + helpPos(4) + yOffset buttonWidth buttonHeight];

okH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', ...
    okPos, 'string', 'OK', 'tooltipstring', 'ok', 'tag', 'ok', 'callback', {@okCallback, InputHandle});

set(subjectH, 'value', 1);
set(componentH, 'value', 1);

% set ok data
regressorInfo = struct('selectedOnset', [], 'selectedRegressor', '');
okData = repmat(struct('regressorInfo', regressorInfo, 'TR', [], 'value', 1), 1, length(refAppData));
set(okH, 'userdata', okData); % set the userdata

for ii = size(listString, 1):-1:1
    set(subjectH, 'value', ii);
    % execute the callbacks
    subjectCallback(subjectH, [], InputHandle);
    if strcmpi(handleVisibility, 'on')
        set(regressH, 'value', okData(ii).value);
    else
        tempStr = get(regressH, 'string');
        set(regressH, 'value', 1:size(tempStr, 1));
    end
    % regressor callback
    regressorCallback(regressH, [], InputHandle);
end

% get the user data
if strcmpi(handleVisibility, 'on')
    
    try
        waitfor(InputHandle);
    catch
        delete(InputHandle);
    end
    
else
    % execute the callback
    okCallback(okH, [], InputHandle);
end
% end for getting the data

% get the application data
if isappdata(0, 'eventAppData')
    % get the application data
    appData = getappdata(0, 'eventAppData');
    rmappdata(0, 'eventAppData');
else
    error('figure window was quit');
end

% set callback only for the regressor listbox

% function callbacks
function subjectCallback(hObject, event_data, handles)

% get the handles data
handles_data = get(handles, 'userdata');
refAppData = handles_data.refAppData;

okH = findobj(handles, 'tag', 'ok');
okData = get(okH, 'userdata');

% value
subjectVal = get(hObject, 'value');

% regressor listbox
regressListH = findobj(handles, 'tag', 'regressorListbox');

% selected reference function names
refNames = refAppData(subjectVal).data.sel_refnames;

set(regressListH, 'string', str2mat(refNames), 'value', okData(subjectVal).value);


function regressorCallback(handleObj, event_data, handles)
% callback for regressor listbox
%

% get user data
handles_data = get(handles, 'userdata');

% get reference application data
refAppData = handles_data.refAppData;

clear handles_data;

okH = findobj(handles, 'tag', 'ok');
okData = get(okH, 'userdata');

% get the onset number and the selected reference function name

getString = get(handleObj, 'string');

getVal = get(handleObj, 'value');

% selected regressor
selectedRefName = deblank(getString(getVal, :));

% subject listbox
subjectH = findobj(handles, 'tag', 'subjectListbox');

% subject value
subjectVal = get(subjectH, 'value');

% spm variable
spmVar = refAppData(subjectVal).data.SPM;

% TR
TimeResponse = spmVar.xY.RT;

% UNITS
UNITS = spmVar.xBF.UNITS;

regressorInfo = repmat(struct('selectedOnset', [], 'selectedRegressor', ''), size(selectedRefName, 1), 1);
% Loop over selected regressors
for nR = 1:size(selectedRefName, 1)
    
    % get the onset number
    [data_sessionNumber, onset_number] = icatb_get_onsets(deblank(selectedRefName(nR, :)), spmVar);
    
    % get the session number
    if data_sessionNumber > length(spmVar.Sess)
        data_sessionNumber = length(spmVar.Sess);
    end
    
    %%%%% Modifying selected onset
    if length(spmVar.Sess(data_sessionNumber).U) ~= 0
        if strcmp(lower(UNITS), 'scans')
            selectedOnset = round(spmVar.Sess(data_sessionNumber).U(onset_number).ons);
        else
            selectedOnset = round((spmVar.Sess(data_sessionNumber).U(onset_number).ons) / TimeResponse);
        end
    else
        error('Cannot calculate event average as onset information is missing');
    end
    
    regressorInfo(nR).selectedOnset = selectedOnset;
    regressorInfo(nR).selectedRegressor = deblank(selectedRefName(nR, :));
    
end
% End loop over selected regressors

% set ok data
%okData(subjectVal).selectedOnset = selectedOnset;
%okData(subjectVal).selectedRegressor = selectedRefName;

okData(subjectVal).regressorInfo = regressorInfo;
okData(subjectVal).TR = TimeResponse;
% store the value in ok data
okData(subjectVal).value = getVal;

set(okH, 'userdata', okData);

function okCallback(hObj, event_data, handles)
% ok callback

% get the user data of the ok pushbutton
okData = get(hObj, 'userdata');

% component list box
compListH = findobj(handles, 'tag', 'componentListbox');
component_number = get(compListH, 'value');

% subject and session listboxes
subjectH = findobj(handles, 'tag', 'subjectListbox');
subject_str = get(subjectH, 'string');
% select for a particular component
try
    
    % check regressor info
    for ii = 1:length(okData)
        %selectedRegressor = okData(ii).selectedRegressor;
        if isempty(okData(ii).regressorInfo(1).selectedRegressor)
            error(['Regressor is not selected for ', deblank(subject_str(ii, :))]);
        end
    end
    % end for checking regressor info
    
    appData.eventData = okData;
    % get the component number
    appData.component_number = component_number;
    
    % set the application data
    setappdata(0, 'eventAppData', appData);
    delete(handles);
    
catch
    disp(lasterr);
    icatb_errorDialog(lasterr, 'Event Average Error', 'modal');
end