function [varargout] = icatb_select_data(varargin)
% procedure: The entries whose data files to be selected are placed in a
% listbox. When clicked on each string a file selection window is opened. Inputs must be in pairs
%
% Input:
% 1. 'title' - specify the title of the figure
% 2. 'num_subjects' - Number of subjects
% 3. 'num_sessions' - Number of sessions
% 4. 'subject_string' - specify the subject string
% 5. 'files_specification' - 'equal' or 'unequal'
% 6. 'type_selection' - 'single' or 'multiple'
% 7. 'spm_check' - pass the time points
%
% Output:
% selectedFiles - structure containing the files of each data set

% check the number of arguments
if mod(nargin, 2) ~= 0
    error('Arguments must be in pairs');
end

windowStyle = 'normal';

% Load defaults
icatb_defaults;

% Screen Color Defaults
global BG_COLOR;
global BG2_COLOR;
global FONT_COLOR;
global BUTTON_COLOR;
global BUTTON_FONT_COLOR;
global AXES_COLOR;
global FONT_COLOR;

% font defaults
global UI_FONTNAME;
global UI_FONTUNITS;
global UI_FS;

% file namings
global READ_NAMING_COMPLEX_IMAGES;

complexType = 'real&imaginary';
complex_file_naming = READ_NAMING_COMPLEX_IMAGES.real_imag;

% Initialise Variables
titleFig = 'Select data files ';
num_subjects = 1;
num_sessions = 1;
files_specification = 'unequal';
type_file_selection = 'multiple';
spm_check = 'no';
filter_string = '*.img';
diffTimePoints = [];
startPath = pwd;
fileInfo = [];
figure_menu = [];
dataType = 'real';
fileType = 'any';

% Get the variables
for ii = 1:2:nargin
    if strcmp(lower(varargin{ii}), 'title')
        titleFig = varargin{ii + 1};
    elseif strcmp(lower(varargin{ii}), 'num_subjects')
        num_subjects = varargin{ii + 1};
    elseif strcmp(lower(varargin{ii}), 'num_sessions')
        num_sessions = varargin{ii + 1};
    elseif strcmp(lower(varargin{ii}), 'files_specification')
        files_specification = lower(varargin{ii + 1});
    elseif strcmp(lower(varargin{ii}), 'type_file_selection')
        type_file_selection = lower(varargin{ii + 1});
    elseif strcmp(lower(varargin{ii}), 'spm_check')
        spm_check = lower(varargin{ii + 1});
    elseif strcmp(lower(varargin{ii}), 'subject_string')
        subject_string = varargin{ii + 1};
    elseif strcmp(lower(varargin{ii}), 'filter_string')
        filter_string = varargin{ii + 1};
    elseif strcmp(lower(varargin{ii}), 'difftimepoints')
        diffTimePoints = varargin{ii + 1};
    elseif strcmp(lower(varargin{ii}), 'startpath')
        startPath = varargin{ii + 1};
    elseif strcmp(lower(varargin{ii}), 'fileinfo')
        fileInfo = varargin{ii + 1};
    elseif strcmp(lower(varargin{ii}), 'figure_menu')
        figure_menu = varargin{ii + 1};
    elseif strcmp(lower(varargin{ii}), 'datatype')
        dataType = lower(varargin{ii + 1});
    elseif strcmp(lower(varargin{ii}), 'complex_file_naming')
        complex_file_naming = varargin{ii + 1};
    elseif strcmp(lower(varargin{ii}), 'read_complex_images')
        complexType = varargin{ii + 1};
    elseif strcmp(lower(varargin{ii}), 'filetype')
        fileType = varargin{ii + 1};
    elseif strcmp(lower(varargin{ii}), 'windowstyle')
        windowStyle = varargin{ii + 1};
    end
end


% define subject string
if ~exist('subject_string', 'var')
    count = 0;
    for ii = 1:num_subjects
        for jj = 1:num_sessions
            count = count + 1;
            subject_string(count).name = ['Subject ', num2str(ii), ' Session ', num2str(jj)];
        end
    end
end

filter_answer = filter_string;

[pathstr2, name2, fileExtension] = fileparts(deblank(filter_string));

% depends on the length of the subject string
dataFiles = repmat(struct('name', []), 1, length(subject_string));
for ii = 1:length(subject_string)
    dataFiles(ii).name = [];
end

dataDim = zeros(length(dataFiles), 3);

% Setup figure for GUI
[InputHandle] = icatb_getGraphics(titleFig, 'normal', titleFig);

set(InputHandle, 'menubar', 'none', 'windowStyle', windowStyle);

[modalityType, dataTitle] = icatb_get_modality;

% Figure Data
figureData = struct('oldDir', pwd, 'pwd', startPath, 'num_subjects', num_subjects, 'num_sessions', num_sessions, ...
    'files_specification', files_specification, 'type_file_selection', type_file_selection, ...
    'spm_check', spm_check, 'filter_string', filter_string, 'diffTimePoints', diffTimePoints, 'dataCount', [], ...
    'changeFiles', 'no', 'dataType', dataType, 'dataDim', dataDim, 'filter_answer', filter_answer, 'modality', modalityType, 'dataTitle', ...
    dataTitle, 'fileType', fileType, 'file_ext', fileExtension);


figureData.complex_file_naming = complex_file_naming;
figureData.complexType = complexType;

if ~isempty(diffTimePoints)
    % Different time points options menu
    menuUserData.subject_string = subject_string;
    menuUserData.diffTimePoints = diffTimePoints;
    menuUserData.numSubjects = num_subjects;
    menuUserData.numSessions = num_sessions;
    guiOptionsMenu = uimenu('parent', InputHandle, 'label', 'GUI Options');
    diffTimePointsMenu = uimenu(guiOptionsMenu, 'label', 'Number of images for each dataset', 'callback', ...
        {@viewTimePointsCallback, InputHandle}, 'userdata', menuUserData);
end

if ~isempty(figure_menu)
    menu1H = uimenu('parent', InputHandle, 'label', 'GIFT-Help');
    if strcmp(lower(figure_menu), 'spm')
        menu2H = uimenu('parent', menu1H, 'label', 'SPM File Window', 'callback', 'icatb_openHTMLHelpFile(''icatb_SPMFileSelection.htm'');');
    else
        menu2H = uimenu('parent', menu1H, 'label', [modalityType, ' Data Window'], 'callback', 'icatb_openHTMLHelpFile(''icatb_gui_data_selection.htm'');');
    end
end

set(InputHandle, 'userdata', figureData);

%%%%%%%%%%%%% Draw Title here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% title color
titleColor = [0 0.9 0.9];

% fonts
titleFont = 13;

axisH = axes('Parent', InputHandle, 'position', [0 0 1 1], 'visible', 'off');

xPos = 0.5; yPos = 0.97;

text(xPos, yPos, titleFig, 'color', titleColor, 'FontAngle', 'italic', 'fontweight', 'bold', 'fontsize', ...
    titleFont, 'HorizontalAlignment', 'center', 'FontName', UI_FONTNAME, 'parent', axisH);

%%%%%%%%%%%%% end for drawing title %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% draw listbox on the left handside
xOffset = 0.04; yOffset = 0.05;

% draw a frame here
frameOrigin = yOffset;
frameWidth = 0.65; frameHeight = 0.06;
framePos = [0.5 - 0.5*frameWidth frameOrigin frameWidth frameHeight];

% plot frame
frameH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'frame', 'position', framePos);

% selected Text Position here
selectedTextPos = framePos;
selectedTextPos(2) = framePos(2) + 0.005; selectedTextPos(4) = framePos(4) - 0.01;
selectedTextPos(1) = framePos(1) + 0.005; selectedTextPos(3) = framePos(3) - 0.01;

% plot selected text here
selectedTextH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'text', 'position', ...
    selectedTextPos, 'string', 'Files Selected: ', 'tag', 'files_selected_text_box');

% Listbox positions
listWidth = 0.65;
listPos = framePos(2) + framePos(4) + yOffset;
listHeight = (yPos - listPos - 2*yOffset);
listPos = [xOffset listPos listWidth listHeight];

% plot the listbox here
listH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'listbox', 'position', listPos, ...
    'string', str2mat(subject_string.name), 'tag', 'data_listbox', 'userdata', dataFiles, 'callback', ...
    {@listCallback, InputHandle});
clear subject_string;

% calculate the positioning of the pushbuttons
pushButtonWidth = 0.2; pushButtonHeight = 0.05;
pushButtonXPos = listPos(1) + listPos(3) + 1.5*xOffset;
pushButtonYPos = listPos(2) + listPos(4) - yOffset;
% push button position
pushPos1 = [pushButtonXPos pushButtonYPos pushButtonWidth pushButtonHeight];
pushPos2 = pushPos1; pushPos2(2) = pushPos2(2) - pushPos2(4) - yOffset;
pushPos3 = pushPos2; pushPos3(2) = pushPos3(2) - pushPos3(4) - yOffset;
pushPos4 = pushPos3; pushPos4(2) = pushPos4(2) - pushPos4(4) - yOffset;
pushPos5 = pushPos4; pushPos5(2) = pushPos5(2) - pushPos5(4) - yOffset;

viewH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'position', pushPos1, 'style', 'pushbutton', ...
    'string', 'View', 'tag', 'view', 'callback', {@viewFilesCallback, listH});

changeH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'position', pushPos2, 'style', 'pushbutton', ...
    'string', 'Change', 'tag', 'change', 'callback', {@changeCallback, listH});

% saveH = uicontrol(InputHandle, 'units', 'normalized', 'position', pushPos3, 'style', 'pushbutton', ...
%     'BackgroundColor', BUTTON_COLOR, 'ForegroundColor', BUTTON_FONT_COLOR, 'fontunits', UI_FONTUNITS, 'fontname', ...
%     UI_FONTNAME, 'FontSize', UI_FS, 'string', 'Save', 'tag', 'save', 'callback', {@saveCallback, listH}, 'userdata', fileInfo);

helpString = sprintf(['Click on the left listbox and select the files for each dataset.', ...
    ' After the selection is done, press Ok button.', ' The selected files for that ' ...
    'data set can be viewed by clicking View button. The selected files can be changed by clicking on change button.']);
helpData.string = helpString; helpData.title =  'Data Selection';

helpH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'position', pushPos3, 'style', 'pushbutton', ...
    'BackgroundColor', BUTTON_COLOR, 'ForegroundColor', [1 1 0], 'fontunits', UI_FONTUNITS, 'fontname', ...
    UI_FONTNAME, 'FontSize', UI_FS, 'string', '?', 'tag', 'help', 'fontweight', 'bold', 'userdata', helpData, ...
    'callback', {@helpCallback, InputHandle});

okH = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'position', pushPos4, 'style', 'pushbutton', ...
    'string', 'OK', 'tag', 'OK', 'callback', {@okCallback, InputHandle});

% set figure visibility on
try
    set(InputHandle, 'visible','on');
    waitfor(InputHandle);
catch
    if ishandle(InputHandle)
        delete(InputHandle);
    end
end

% get the application data
if isappdata(0, 'okAppData')
    appData = getappdata(0, 'okAppData');
    varargout{1} = appData.dataFiles;
    varargout{2} = appData.dataCount;
    % remove application data
    rmappdata(0, 'okAppData');
    cd(figureData.oldDir);
else
    cd(figureData.oldDir);
    varargout = {};
    error('Figure Window was quit');
end

%%%%%%%%%%%%%%%% Function callbacks %%%%%%%%%%%%%%%%
function listCallback(hObject, evd, handles)
% 1. listbox callback
% purpose: open File selection window when clicked on the string
% use the defaults from the figure data
try
    set(handles, 'pointer', 'watch');
    % get the userdata
    figureData = get(handles, 'userdata');
    numOfSub = figureData.num_subjects; % number of subjects
    numOfSess = figureData.num_sessions; % number of sessions
    files_specification = figureData.files_specification; % files specification (equal or unequal)
    type_file_selection = figureData.type_file_selection; % type selection (single or multiple)
    spm_check = figureData.spm_check; % check spm design matrix
    filter_string = figureData.filter_string; % filter string
    dataCount = figureData.dataCount;
    % check the data type
    dataType = figureData.dataType;
    % start path for the select entry
    currentDir = figureData.pwd;
    diffTimePoints = figureData.diffTimePoints;
    % get string
    getValue = get(hObject, 'value');
    getString = get(hObject, 'string');
    selectedString = getString(getValue, :);
    % get listbox userdata
    dataFiles = get(hObject, 'userdata'); % data files information
    complexType = figureData.complexType;
    dataDim = figureData.dataDim; % data dimensions
    filter_answer = figureData.filter_answer;
    % data title
    dataTitle = lower(figureData.dataTitle);
    fileType = figureData.fileType;
    file_ext = figureData.file_ext;
    complex_file_naming = figureData.complex_file_naming;
    
    % real data
    if strcmp(dataType, 'real')
        
        if strcmp(lower(type_file_selection), 'multiple')
            promptString = ['Select files for ' deblank(getString(getValue, :))];
        else
            promptString = ['Select a file for ' deblank(getString(getValue, :))];
        end
        
    else
        % imaginary data
        if strcmp(lower(type_file_selection), 'multiple')
            if strcmpi(complexType, 'real&imaginary')
                promptString = ['Select imaginary data files for ' deblank(getString(getValue, :))];
            else
                promptString = ['Select phase data files for ' deblank(getString(getValue, :))];
            end
        else
            promptString = ['Select a file for ' deblank(getString(getValue, :))];
        end
    end
    
    if isempty(dataFiles(getValue).name) || ~strcmp(lower(figureData.changeFiles), 'no')
        try
            tempDir = fullfile(currentDir, '..');
            cd(tempDir);
            tempDir = pwd;
        catch
            tempDir = currentDir;
        end
        % New file selector window
        [dataFiles(getValue).name, filter_answer] = icatb_selectEntry('typeEntity', 'file', 'typeselection', type_file_selection, ...
            'num_datasets', 2, 'count_datasets', 1, 'filter', filter_answer, 'title', promptString, 'startPath', ...
            tempDir, 'fileType', fileType);
    end
    
    if isempty(dataFiles(getValue).name)
        error(['Data for dataset ', num2str(getValue), ' is not selected']);
    end
    
    % Convert files to cell array
    checkFiles = cellstr(dataFiles(getValue).name);
    
    numberFiles = length(checkFiles);
    
    % Check the file extension
    if strcmpi(fileType, 'image')
        % Check if the files are in IMG or NII format
        good_inds = findGoodCells(checkFiles, '\.nii|\.img');
        
        if length(good_inds) ~= length(checkFiles)
            error('Please check the selected image files as they should be in .img or .nii format.');
        end
        % End for checking the files in IMG or NII format
        
        % Get only the first file volume
        [nn, mm, dims] = icatb_get_countTimePoints(icatb_parseExtn(checkFiles{1}));
        
        % Check if the images start with I_ or end with _I
        if strcmp(dataType, 'complex')
            % get file names not full paths
            [startInd] = regexp(checkFiles, ['\', filesep]);
            fileNames = cell(length(checkFiles), 1);
            for nF = 1:length(fileNames)
                chk = checkFiles{nF}(startInd{nF}(end) + 1 :end);
                checkExtn = icatb_findstr(chk, '.');
                if ~isempty(checkExtn)
                    fileNames{nF} = chk(1:checkExtn(end)-1);
                else
                    fileNames{nF} = chk;
                end
            end
            
            tempVar = complex_file_naming{2};
            
            % Complex images should start either with I_ or end with _I
            good_inds = findGoodCells(fileNames, ['^', tempVar, '|', '.*', tempVar(end:-1:1), '$']);
            if length(good_inds) ~= length(fileNames)
                error('Please check the READ_NAMING_COMPLEX_IMAGES variable in icatb_defaults.m.');
            end
        end
        % End for checking the images
        
    else
        % Check other formats for files
        good_inds = findGoodCells(checkFiles, ['\', file_ext, '$']);
        if length(good_inds) ~= length(checkFiles)
            error(['Please check the file extension of the selected files. They should be in ', file_ext, ' format']);
        end
        dims = [];
        
    end
    % End for checking the data
    
    % find the number of files selected textbox
    selectedTextBoxH = findobj(handles, 'tag', 'files_selected_text_box');
    set(selectedTextBoxH, 'string', ['Selected Files: ', num2str(numberFiles)]);
    dataCount(getValue) = numberFiles;
    
    if ~isempty(dims)
        dataDim(getValue, :) = dims;
    end
    % end for selecting design matrix or spm design matrix
    
    if ~strcmp(spm_check, 'no')
        try
            startTp = (getValue - 1)*numOfSess + 1;
            endTp = getValue*numOfSess;
            for ii = 1:numOfSess
                % design information
                [spmData] = icatb_loadSPM_new('spmName', dataFiles(getValue).name, 'countTimePoints', ...
                    diffTimePoints(startTp:endTp), 'flag_selecting_regressors', 'no', 'data_sessionNumber', ii); % model Index doesn't contain actual indices
            end
        catch
            error(['SPM design matrix for ', selectedString, ' doesn''t match with the ', dataTitle, ' data.']);
        end
    end
    
    % store the current directory in figureData
    figureData.filter_answer = filter_answer;
    figureData.pwd = pwd;
    figureData.dataCount = dataCount;
    figureData.changeFiles = 'no';
    figureData.dataDim = dataDim;
    
    % figure data
    set(handles, 'userdata', figureData);
    
    % userdata for listbox
    set(hObject, 'userdata', dataFiles);
    
    set(handles, 'pointer', 'arrow');
    
catch
    set(handles, 'pointer', 'arrow');
    icatb_errorDialog(lasterr, 'File selection error');
    disp(lasterr);
end

% 2. Ok button callback
function okCallback(hObject, evd, handles)

try
    % get figureData
    figureData = get(handles, 'userdata');
    numOfSub = figureData.num_subjects; % number of subjects
    numOfSess = figureData.num_sessions; % number of sessions
    files_specification = figureData.files_specification; % files specification (equal or unequal)
    type_file_selection = figureData.type_file_selection; % type selection (single or multiple)
    spm_check = figureData.spm_check; % check spm design matrix
    filter_string = figureData.filter_string; % filter string
    diffTimePoints = figureData.diffTimePoints;
    % check the data type
    dataType = figureData.dataType;
    % file naming
    complex_file_naming = figureData.complex_file_naming;
    % get the listbox data
    listH = findobj(handles, 'tag', 'data_listbox');
    listString = get(listH, 'string');
    dataFiles = get(listH, 'userdata');
    
    if isempty(dataFiles)
        error('Data is not selected');
    end
    % Loop over the number of data sets
    %ii = 0;
    % Initialise number of files
    numFiles = zeros(1, length(dataFiles));
    % loop over subjects
    for ii = 1:length(dataFiles)
        % check data
        P = dataFiles(ii).name;
        % check if the data is empty
        if isempty(P)
            error(['Data for ', listString(ii, :), ' is not selected']);
        end
        
        % use the new function to get the time points
        %[numFiles(ii), extns] = icatb_get_countTimePoints(P); %size(P, 1);
        numFiles(ii) = figureData.dataCount(ii);
        
        %%%%%%%%%%%%%%%%%%%% check data dimensions %%%%%%%%%%%
        if ii == 1
            oldDim = figureData.dataDim(ii, :);
        else
            newDim = figureData.dataDim(ii, :); % new dimension
            checkDim = find((newDim == oldDim) > 0);
            if length(checkDim) ~= length(oldDim)
                error(['Data dimensions should be the same for all subjects. Please check the dimension for ', ...
                    deblank(listString(ii, :))]);
            end
        end
        
        %%%%%%%%%%%%%%%%%% end for checking the dimensions %%%%%%%%%
        
        
        %% check the count of the files
        if strcmp(lower(files_specification), 'equal')
            if ii > 1
                previousLength = numFiles(ii - 1);
                if previousLength ~= numFiles(ii)
                    error(['Same number of files should be selected for ', listString(ii, :)]);
                end
            end
        end
    end
    appData.dataFiles = dataFiles;
    clear dataFiles;
    appData.dataCount = figureData.dataCount;
    % set application data
    setappdata(0, 'okAppData', appData);
    delete(handles);
    
catch
    errorH = errordlg(lasterr);
    set(errorH, 'windowstyle', 'modal');
    disp(lasterr);
end

% 3. View Files callback
function viewFilesCallback(hObject, event_data, handles)

% View files

dataFiles = get(handles, 'userdata');
listValue = get(handles, 'value');
listString = get(handles, 'string');
selectedString = listString(listValue, :);
textBody = dataFiles(listValue).name;
if isempty(textBody)
    textBody = '';
end
titleFig = ['Files for ', selectedString];
% open dialog box
figHandle = icatb_dialogBox('title', titleFig, 'textBody', textBody, 'textType', 'large');
waitfor(figHandle);

% 4. Help callback
function helpCallback(hObject, evd, handles)

getUserdata = get(hObject, 'userdata');

% open dialog box
figHandle = icatb_dialogBox('title', getUserdata.title, 'textBody', getUserdata.string, 'textType', 'large');

waitfor(figHandle);

function changeCallback(hObject, event_data, handles)
% change the selected files
figH = get(hObject, 'parent'); %get the figure handle
figureData = get(figH, 'userdata'); % get the figure data
figureData.changeFiles = 'yes'; % pass the flag to change the files
set(figH, 'userdata', figureData); % set the figure data

listCallback(handles, [], figH); % use listbox callback to select the files

function saveCallback(hObject, event_data, handles)

% save the selected files information
figH = get(hObject, 'parent'); % get the figure
figData = get(figH, 'userdata'); % get the figure data
numOfSub = figData.num_subjects;
numOfSess = figData.num_sessions;
fileName = fullfile(figData.oldDir, 'Subject.mat');
fileFormat = '';
userData = get(hObject, 'userdata'); % get user data

if isfield(userData, 'format')
    % File informations
    fileFormat = userData.format;
end

if isfield(userData, 'fileName')
    % file name
    fileName = userData.fileName;
else
    prompt = {'Enter Valid File Name:'};
    dlg_title = 'Save Data as';
    num_lines = 1;
    def = {''};
    % save the file with the file name specified
    fileName = icatb_inputdlg2(prompt, dlg_title, num_lines, def);
    % make a full file
    fileName = fullfile(figData.oldDir, [fileName{1}, 'Subject', '.mat']);
end

% save the information in a file
if strcmp(lower(fileFormat), 'append')
    SPMFiles = get(handles, 'userdata');
    icatb_save(fileName, 'SPMFiles', '-append');
else
    files = get(handles, 'userdata');
    SPMFiles.name = [];
    icatb_save(fileName, 'files', 'SPMFiles', 'numOfSub', 'numOfSess');
end

% display the information
dispString = ['Data information is stored in ', fileName];
disp(dispString);
% open dialog box
figHandle = icatb_dialogBox('title', 'File Saved', 'textBody', dispString, 'textType', 'large');
waitfor(figHandle);

function viewTimePointsCallback(hObject, event_data, handles)
% view different time points
% set the figure style normal

menuUserData = get(hObject, 'userdata');
subject_string = menuUserData.subject_string;
diffTimePoints = menuUserData.diffTimePoints;
num_subjects = menuUserData.numSubjects;
num_sessions = menuUserData.numSessions;

% loop over subject string
for ii = 1:length(subject_string)
    startTp = (ii - 1)*num_sessions + 1; endTp = ii*num_sessions;
    subject_string(ii).name = [subject_string(ii).name, ': ', '[', num2str(diffTimePoints(startTp:endTp)), ']'];
end

dialogHandle = icatb_dialogBox('title', 'Information on number of images for each dataset', 'textBody', ...
    str2mat(subject_string.name), 'textType', 'large', 'windowstyle', 'normal');


function good_inds = findGoodCells(checkFiles, checkRegExp)
% Match regular expression

checkExtn = regexpi(checkFiles, checkRegExp);
good_inds = icatb_good_cells(checkExtn);
good_inds = find(good_inds ~= 0);