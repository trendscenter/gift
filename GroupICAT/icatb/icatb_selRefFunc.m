function [refAppdata] = icatb_selRefFunc(varargin)
% Select the regressors based on the given information
% case 1: When the regressor Names are given and Indices are not given
% case 2: When the regressor Nmes and Indices are given
%
% Input:
%
% Output:
% needed fields:
% selected timecourse
% selectedRegressorNames
% selectedIndices (actual)
% OriginalIndices
% OriginalRegressorNames

% defaults:
typeSelection = 'single';
helpDetails.str = sprintf(['Click on the left listbox and the corresponding reference functions will be displayed at', ...
    ' the bottom. The selected reference functions from the bottom listbox will be displayed at the right listbox.', ...
    ' After the selection process is done press the Ok button.']);
titleFig = 'Reference Function Window';
helpDetails.title = titleFig;
storeModelTimecourse = 'no';

% default for figure handle
handleVisibility = 'on';
% regressors for subjects and sessions
regress_subjects_sessions = [];
subjectNumber = [];
data_sessionNumber = [];
% loop over the number of args
for ii = 1:2:nargin
    if strcmp(lower(varargin{ii}), 'spmnames')
        spmNames = varargin{ii + 1}; % spm names
    elseif strcmp(lower(varargin{ii}), 'original_indices')
        originalIndices = varargin{ii + 1}; % original indices
    elseif strcmp(lower(varargin{ii}), 'substring')
        sub_string = varargin{ii + 1}; % subject listbox string
    elseif strcmp(lower(varargin{ii}), 'numsubjects')
        numSubjects = varargin{ii + 1}; % number of subjects
    elseif strcmp(lower(varargin{ii}), 'numsessions')
        numSessions = varargin{ii + 1}; % number of sessions
    elseif strcmp(lower(varargin{ii}), 'help_details')
        helpDetails = varargin{ii + 1}; % help string for the help button
    elseif strcmp(lower(varargin{ii}), 'count_time_points')
        % variable useful to compare the ICA timepoints with the model time
        % points
        countTimePoints = varargin{ii + 1}; % count for time points
    elseif strcmp(lower(varargin{ii}), 'store_model_timecourse')
        storeModelTimecourse = lower(varargin{ii + 1}); % flag for storing time course in the structure
    elseif strcmp(lower(varargin{ii}), 'title_figure')
        titleFig = varargin{ii + 1}; % Figure tag and figure title
    elseif strcmp(lower(varargin{ii}), 'typeselection')
        typeSelection = varargin{ii + 1};
    elseif strcmp(lower(varargin{ii}), 'data_sessionnumber')
        data_sessionNumber = varargin{ii + 1};
    elseif strcmp(lower(varargin{ii}), 'all_time_points')
        allTimePoints = varargin{ii + 1};
    elseif strcmp(lower(varargin{ii}), 'numofsub')
        numOfSub = varargin{ii + 1};
    elseif strcmp(lower(varargin{ii}), 'numofsess')
        numOfSess = varargin{ii + 1};

        %%%%%%%% Below vars are required for loading reference functions from a
        %%%%%%%% text file
    elseif strcmpi(varargin{ii}, 'handle_visibility')
        % check the  handle visibility
        handleVisibility = lower(varargin{ii + 1});
    elseif strcmpi(varargin{ii}, 'txtfile')
        % sorting text file
        sortingTextFile = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'subjectnumber')
        % subject number
        subjectNumber = varargin{ii + 1};
        %     elseif strcmpi(varargin{ii}, 'keywd')
        %         % keyword for data sets
        %         keywd = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'spmmatflag')
        spmMatFlag = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'flaguseregressors')
        flagUseRegressors = varargin{ii + 1};
    end
end

if isempty(data_sessionNumber)
    clear data_sessionNumber;
end

% Check the existence of the number of subjects and sessions
if ~exist('numSubjects', 'var')
    numSubjects = 1;
end

if ~exist('numSessions', 'var')
    numSessions = 1;
end

% caculate data sets
numDataSets = numSubjects*numSessions;

% check the spm names
if size(spmNames, 1) > numDataSets
    error('Number of spmNames should be the same as number of data sets');
end

% check the subject string
if ~exist('sub_string', 'var')
    count = 0;
    for ii = 1:numSubjects
        for jj = 1:numSessions
            count = count + 1;
            sub_string(count).name = ['Subject ', num2str(ii), ' Session ', num2str(jj)]; % default subject string
        end
    end
end

if strcmp(typeSelection, 'single')
    maxEntries = 1;
else
    maxEntries = 2;
end

% create SPM application data
for ii = 1:numSubjects*numSessions
    spmAppData(ii).data = [];
end

% Procedure: draw left top and bottom listboxes and right top listbox and
% beneath right top listbox draw textbox and disable the textbox but update
% the textbox with the number of entries selected.
% Fields included in the application data are
% as follows:
% selected timecourse - set as appdata
% selectedRegressorNames - set as appdata
% selectedIndices (actual) - set as appdata
% OriginalIndices - set as figure data
% OriginalRegressorNames - set as figure data

% Load defaults
icatb_defaults;

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

% Number of reference functions
dataCount = 0;
for numSub = 1:numSubjects
    startTp = (numSub - 1)*numSessions + 1;
    endTp = numSub*numSessions;
    for numSess = 1:numSessions
        dataCount = dataCount + 1;
        if exist('data_sessionNumber', 'var')
            sessNumber = data_sessionNumber;
            tempTimePoints = zeros(1, numOfSess);
            for nn = 1:numOfSess
                tempTimePoints(nn) = allTimePoints.sub(numSub).sess(nn);
            end
        else
            sessNumber = numSess;
            tempTimePoints = countTimePoints(startTp:endTp);
        end

        % get the subject number
        if isempty(subjectNumber)
            subNumber = numSub;
        else
            subNumber = subjectNumber;
        end

        % load SPM design matrix
        [spmAppData(dataCount).data] = icatb_loadSPM_new('spmName', deblank(spmNames(dataCount, :)), ...
            'countTimePoints', tempTimePoints, 'data_sessionNumber', sessNumber, 'typeStr', typeSelection, ...
            'check_spm_design_matrix', 'yes', 'check_regressors', 'yes');

        % figure visibility
        if ~strcmpi(handleVisibility, 'on')
            % flag for using regressors
            if strcmpi(flagUseRegressors, 'all')
                regress_subjects_sessions(dataCount).names = cellstr(spmAppData(dataCount).data.sel_refnames) ;
                if strcmpi(typeSelection, 'single')
                    regress_subjects_sessions(dataCount).names = {deblank(regress_subjects_sessions(dataCount).names(1, :))};
                    % print regressor selected
                    disp(['Regressor selected for data-set ', num2str(dataCount), ' is ', ...
                        str2mat(regress_subjects_sessions(dataCount).names{1})]);
                end
            else
                % read regressors from file
                [regress_sub_sess] = icatb_readRegressors_file('txtfile', deblank(sortingTextFile), 'numOfSub', numOfSub, ...
                    'numOfSess', numOfSess, 'spmmatflag', spmMatFlag, 'keywd', 'single_subject_session', 'sessionNumber', ...
                    sessNumber, 'subjectNumber', subNumber, 'typeSelection', typeSelection);
                % add it to a strcuture
                [regress_subjects_sessions(dataCount).names] = regress_sub_sess.sub.sess.names;
                clear regress_sub_sess;
            end
            % end for checking
        end
        % end for figure visibility
        clear tempTimePoints;
        spmAppData(dataCount).data.designMatrixName = deblank(spmNames(dataCount, :));
        % all the regressors
        allNames = spmAppData(dataCount).data.sel_refnames; %SPM.xX.name;

        % set the reference names and the indices of the regressors selected
        if exist('originalIndices', 'var')
            refNames = str2mat(allNames(originalIndices(dataCount, :), :));
            refData(dataCount).originalIndices = (1:size(originalIndices, 2)); %originalIndices(dataCount, :);
        else
            refNames = allNames;
        end
        % store the reference function names in a structure
        refData(dataCount).refNames = refNames;
        clear refNames;
    end
end

% create three listboxes one for subjects, one for reference functions,
% one for selected reference functions
subjectString = str2mat(sub_string.name);

clear sub_string;

% Figure
[InputHandle] = icatb_getGraphics(titleFig, 'normal', titleFig, handleVisibility);

% Set no menu bar for the figure
set(InputHandle, 'Menubar', 'None');

% tolerance
ytol = 0.03; xtol = 0.05;

%%%%%%%%%%% ok positions %%%%%%%%%%%%%%%
% define push button positions here
okWidth = 0.2; okHeight = 0.05;

okPos = [0.75 - 0.5*okWidth okHeight + 0.5*ytol  okWidth okHeight];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%% Help Position %%%%%%%%%%%%%%%%%%%
helpPos = [okPos(1) okPos(2) + okPos(4) + ytol okWidth okHeight];

%%%% selected Text Position %%%%%%%%%%%%%%%%%%%%%
selectedTextPos = helpPos;
selectedTextPos(2) = helpPos(2) + helpPos(4) + ytol;
selectedTextPos(4) = 0.05;

listHeight = 0.4;

listWidth = 0.42;

%%%%%%%%%%%%% Text height width %%%%%%%%%%%%%%%%%
textWidth = listWidth;
textHeight = 0.05;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%% Subject Listbox Positions %%%%%%%%%%%%%%%%%%%%%%%%%%%

% subject text position
subjectTextPos = [xtol 0.95 - 0.5*textHeight textWidth textHeight];

listboxPos = selectedTextPos(4) + selectedTextPos(2) + ytol; %okPos(4) + okPos(2) + ytol;

subjectPos = [xtol  subjectTextPos(2) - listHeight - ytol listWidth listHeight];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% Reference Listbox Positions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reference Text Position
refTextPos = [subjectPos(1) subjectPos(2) - ytol - textHeight textWidth textHeight];

referencePos = [xtol  ytol listWidth (refTextPos(2) - 2*ytol)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% Selected Reference Listbox Positions
% Selected Reference Text Position
selRefTextPos = [xtol + subjectPos(1) + subjectPos(3) subjectTextPos(2) textWidth textHeight];

% selected Text Position
selectedTextPos(1) = selRefTextPos(1);
selectedTextPos(3) = listWidth;

listHeight = (subjectTextPos(2) - listboxPos - ytol);

selReferencePos = [xtol + subjectPos(1) + subjectPos(3)  listboxPos listWidth listHeight];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%% Plot UI Controls %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot subject Text and listbox
subjectTextH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', ...
    subjectTextPos, 'horizontalalignment', 'center', 'string', 'Subjects', 'userdata', ...
    storeModelTimecourse, 'tag', 'subject_text');

subjectH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'listbox', 'position', subjectPos, ...
    'horizontalalignment', 'left', 'min', 0,  'max', 1, 'string', subjectString, 'callback', ...
    {@subjectCallback, InputHandle}, 'tag', 'subjectListbox');

% get the extent of the text under static text
getExtent = get(subjectTextH, 'extent');

if subjectTextPos(3) < getExtent(3)
    subjectTextPos(3) = getExtent(3);
end

set(subjectTextH, 'position', subjectTextPos);

% Plot Reference text and listbox
refTextH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', refTextPos, ...
    'horizontalalignment', 'center', 'string', 'Ref. Functions');
refListH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'listbox', 'position', referencePos, ...
    'horizontalalignment', 'left', 'min', 0, 'max', maxEntries, 'callback', {@refListCallback, InputHandle}, ...
    'tag', 'referenceListbox', 'userdata', regress_subjects_sessions);

% get the extent of the text under static text
getExtent = get(refTextH, 'extent');

if refTextPos(3) < getExtent(3)
    refTextPos(3) = getExtent(3);
end

set(refTextH, 'position', refTextPos);

% Plot Selected text and listbox
selrefTextH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', selRefTextPos, ...
    'horizontalalignment', 'center', 'string', 'Selected Ref Function');
selrefListH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'listbox', 'position', selReferencePos, ...
    'horizontalalignment', 'left', 'tag', 'selectListbox');
% selected text box
selectedTextH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', selectedTextPos, ...
    'horizontalalignment', 'center', 'tag', 'selectedTextBox');

% get the extent of the text under static text
getExtent = get(selrefTextH, 'extent');

if selRefTextPos(3) < getExtent(3)
    selRefTextPos(3) = getExtent(3);
end

set(selrefTextH, 'position', selRefTextPos);

% Draw push button
okH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', okPos, 'string', 'OK', ...
    'tooltipstring', 'Select Entries...', 'callback', {@okCallback, InputHandle}, 'tag', 'ok');

% Help push button
helpH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', helpPos, 'string', '?', 'BackgroundColor', BUTTON_COLOR,...
    'ForegroundColor', HELP_FONT_COLOR, 'fontweight', 'bold', 'fontunits', UI_FONTUNITS, 'fontname', UI_FONTNAME, 'FontSize', UI_FS, ...
    'tooltipstring', 'help', 'tag', 'help', 'callback', {@helpCallback, InputHandle}, ...
    'userdata', helpDetails);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set the userdata to the figure
set(InputHandle, 'userdata', refData);

% set the application data here
for ii = 1:numSubjects*numSessions
    refAppData(ii).name = [];
end

set(okH, 'userdata', refAppData);
set(subjectH, 'userdata', spmAppData); % store the spm information in subject listbox

% % perform the callback action here
% subjectCallback(subjectH, [], InputHandle);
%
% refListCallback(refListH, [], InputHandle)


% modify the callback here:
% store the regressors and selected indices here
% and set the listbox to the corresponding indices in case for
% batch files.


%set(InputHandle, 'visible', handleVisibility);

% if the figure is visible
if strcmpi(handleVisibility, 'on')

    for ii = 1:numSubjects*numSessions
        set(subjectH, 'value', ii);
        % perform the callback action here
        subjectCallback(subjectH, [], InputHandle);
        refListCallback(refListH, [], InputHandle);
    end

    set(subjectH, 'value', 1);
    % call the reference listbox again
    subjectCallback(subjectH, [], InputHandle);
    refListCallback(refListH, [], InputHandle);

    clear refData;

    clear refAppData;

    % get the user data from the figure
    try
        uiwait(InputHandle);
    catch
        if ishandle(InputHandle)
            delete(InputHandle);
        end
    end
    % end for try statement

else

    % get the reference functions using text file and store the values and regressors
    % in refAppData

    for ii = 1:numSubjects*numSessions
        set(subjectH, 'value', ii);
        % perform the callback action here
        subjectCallback(subjectH, [], InputHandle);
        refListCallback(refListH, [], InputHandle);
    end

    % execute the ok callback
    okCallback(okH, [], InputHandle);
    %delete(InputHandle);

end

% application data
if isappdata(0, 'refAppData')
    refAppdata = getappdata(0, 'refAppData');
    rmappdata(0, 'refAppData');
else
    error('Figure window was quit');
end

%%%%%%%%%%%%% Callbacks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define callbacks here

% Subject Listbox
function subjectCallback(handleObj, evd, handles)

% Purpose: get the string under selection and update the reference function
% listbox and the selected reference listbox

% Reference List box
refListH = findobj(handles, 'tag', 'referenceListbox');

% selection Listbox
selectH = findobj(handles, 'tag', 'selectListbox');

set(refListH, 'value', 1);
% Initialise reference Listbox
set(refListH, 'string', []);

% Initialise select String Listbox
set(selectH, 'value', 1);
set(selectH, 'string', []);

% userdata
refData = get(handles, 'userdata');

% get the values
getValue = get(handleObj, 'value');

% set the string of reference listbox to the selected subject string
set(refListH, 'string', str2mat(refData(getValue).refNames));

clear refData;

refAppData = get(findobj(handles, 'tag', 'ok'), 'userdata');
selListH = findobj(handles, 'tag', 'selectListbox');
set(selListH, 'string', refAppData(getValue).name);
set(findobj(handles, 'tag', 'selectedTextBox'), 'string', ...
    ['Selected Ref. Func: ', num2str(size(refAppData(getValue).name, 1))]);

clear refAppData;

% Reference Listbox
function refListCallback(handleObj, event_data, handles)

% Purpose: get the string under selection and update the selected reference
% function

% figure data
figureData = get(handles, 'userdata');

% object data
regress_subjects_sessions = get(handleObj, 'userdata');

% get the handle visibility
handleVisibility = get(handles, 'visible');

subListbox = findobj(handles, 'tag', 'subjectListbox');

% all the strings from the subject listbox
allStrings_subListbox = get(subListbox, 'string');

% get the subject under consideration
getValue = get(subListbox, 'value');

% get spm application data
spmAppData = get(subListbox, 'userdata');

% get the string
getString = get(handleObj, 'string');

% get the values
getselValue = get(handleObj, 'value');


% check the handle visibility
if ~strcmpi(handleVisibility, 'on')
    % get the data to compare the regressors from the text file with the
    % the current regressors

    % compare the regressors with the original refNames
    regressors_from_file = regress_subjects_sessions(getValue).names;

    %%%%%%%%% modify the get sel value %%%%%%%%%%%
    % Initialise selected values
    getselValue = zeros(1, length(regressors_from_file));
    % loop over selected regressors
    for ii = 1:length(regressors_from_file)
        % initialise index
        tempIndex = [];
        % loop over all the reference function names
        for jj = 1:size(getString, 1)
            % remove the trailing parts
            currentRefName = deblank(getString(jj, :));
            %tempIndex = strmatch(lower(refFunNames), lower(allNames), 'exact');
            if strcmpi(currentRefName, regressors_from_file{ii})
                %error(['Specified reference function with name ', lower(refFunNames{ii}), ' doesn''t exist']);
                tempIndex = jj;
                break;
            end
        end
        % end for loop for allnames
        if isempty(tempIndex)
            error(['Specified reference function with name ', lower(regressors_from_file{ii}), ' doesn''t exist', ...
                ' for ', deblank(allStrings_subListbox(getValue, :))]);
        end

        % store the indices
        getselValue(ii) = tempIndex;
    end

end

selectStrings = getString(getselValue, :);

% get the userdata
refAppData = get(findobj(handles, 'tag', 'ok'), 'userdata');

if ~isempty(refAppData)
    refAppData(getValue).name = selectStrings;
    if isfield(figureData, 'originalIndices')
        % get the actual value of the indices
        refAppData(getValue).value = figureData(getValue).originalIndices(getselValue);
    else
        refAppData(getValue).value = getselValue;
    end

    % store the following to the spm data as well
    spmAppData(getValue).data.getIndex = getselValue;
    spmAppData(getValue).data.selectedRegressors = spmAppData(getValue).data.sel_refnames(getselValue, :); % store the selected regressors
    spmAppData(getValue).data.nsess = length(spmAppData(getValue).data.nscan); % get the number of sessions in design matrix
    % set the string
    set(findobj(handles, 'tag', 'selectListbox'), 'string', refAppData(getValue).name);
    % show the number of selected entries in the text box
    set(findobj(handles, 'tag', 'selectedTextBox'), 'string', ['Selected Ref. Func: ', ...
        num2str(length(getselValue))]);
    set(findobj(handles, 'tag', 'ok'), 'userdata', refAppData);
    set(subListbox, 'userdata', spmAppData); % set the application data
end

% Ok callback
function okCallback(handleObj, evd, handles)

% get the selected reference functions and delete the current figure

try
    % get the userdata
    refAppData = get(handleObj, 'userdata');

    subListH = findobj(handles, 'tag', 'subjectListbox'); % subject listbox

    spmAppData = get(subListH, 'userdata'); % listbox user data

    maxEntries = get(findobj(handles, 'tag', 'referenceListbox'), 'max');

    subListStrings = get(subListH, 'string'); % get the strings of the subject listbox

    subjectTextH = findobj(handles, 'tag', 'subject_text'); % get the subject text

    % check if any of the selected model time courses are not selected
    for ii = 1:size(refAppData, 2)
        if isempty(refAppData(ii).name)
            error(['Reference Functions are not selected for data set ', subListStrings(ii, :)]);
        end
        % check if it is multiple regression criteria
        if ii > 1 & maxEntries > 1
            previousLength = size(refAppData(ii - 1).name, 1);
            currentLength = size(refAppData(ii).name, 1);
            if previousLength ~= currentLength
                error(['Select same number of regressors for data set ', subListStrings(ii, :)]);
            end
        end
    end

    storeModelTimecourse = get(subjectTextH, 'userdata'); % get the flag regarding storing time courses
    if ~strcmp(lower(storeModelTimecourse), 'no')
        % loop over all the data sets
        for numSub = 1:size(spmAppData, 2)
            spmData = spmAppData(numSub).data;
            spmName = spmData.designMatrixName;
            % get the specified time courses
            [spmAppData(numSub).data] = icatb_loadSPM_new('get_timecourses', 'yes', 'check_spm_design_matrix', 'no', ...
                'check_regressors', 'no', 'flag_selecting_regressors', 'no', 'spmData', spmAppData(numSub).data);
            % store time course, get index and other fields
            refAppData(numSub).spmData = spmAppData(numSub).data;
            clear tempTimecourse;
            clear timecourseIndices;
        end
    end
    % set the application data
    setappdata(0, 'refAppData', refAppData);
    set(handleObj, 'userdata', []);
    set(subListH, 'userdata', []); % remove the user data
    delete(handles);
catch
    icatb_errorDialog(lasterr, 'Reference Function Err');
end

function helpCallback(hObject, evd, handles)

% help button callback
helpDetails = get(hObject, 'userdata');

str = helpDetails.str;
titleString = helpDetails.title;

figHandle = icatb_dialogBox('title', titleString, 'textBody', str, 'textType', 'large');

waitfor(figHandle);
