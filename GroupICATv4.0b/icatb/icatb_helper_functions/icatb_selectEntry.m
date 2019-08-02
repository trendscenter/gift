function [outputText, filterAnswer] = icatb_selectEntry(varargin)
%% Purpose: Select file/files or directory
% Inputs: Arguments must be in pairs
% 1. keyword: 'title'
% title - title of the figure
% 2. keyword:  'typentity'
% typeEntity - 'file' or 'directory'
% 3. keyword:  'typeSelection'
% typeSelection - 'single' or 'multiple'
% 4. keyword: 'filter'
% FilterText - Specify filter (*, *.img, etc).
% 5. keyword: 'num_datasets'
% numDataSets - Number of data sets(number of subjects * number of
% sessions) = 1 by default.
% 6. keyword: 'count_datasets'
% count_datasets - current number of data set: Count number for data set: useful when more than one
% % data file is to be selected
% 7. keyword: 'startpath'
% outputDir - specify starting path

% Output: Full path for the files or the directory selected

%% Initialise vars
typeEntity = 'file'; % type entity
typeSelection = 'single'; % selection type
figName = 'Directory Selection Window';  % Name for the figure
numDataSets = 1; % default data sets
count_datasets = 1; % default count data sets
outputDir = pwd; % default output directory
FilterText = '*'; % default filter text
selectMode = 1; % Selection mode is single

filterAnswer = [];

fileType = 'any';
fileNumbers = '';

%% Loading the defaults from the defaults file
icatb_defaults;

% Screen Color Defaults
global BG2_COLOR;
global FONT_COLOR;
global BUTTON_COLOR;
global BUTTON_FONT_COLOR;

% font size and name
global UI_FONTNAME;
global UI_FONTUNITS;
global UI_FS;

%% loop over the number of input arguments passed
for ii = 1:2:nargin
    if strcmpi(varargin{ii}, 'title')
        titleFig = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'typeentity')
        typeEntity = lower(varargin{ii + 1});
    elseif strcmpi(varargin{ii}, 'typeselection')
        typeSelection = lower(varargin{ii + 1});
    elseif strcmpi(varargin{ii}, 'filter')
        FilterText = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'num_datasets')
        numDataSets = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'count_datasets')
        count_datasets = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'startpath')
        outputDir = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'filetype')
        fileType = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'filenumbers')
        fileNumbers = varargin{ii + 1};
    end
end
%% End for loop over the number of input arguments passed

% Check the existence of typeEntity
if ~strcmpi(typeEntity, 'file') && ~strcmpi(typeEntity, 'directory')
    typeEntity = 'file';
end

% figure window name
if strcmpi(typeEntity, 'file')
    figName = 'File Selection Window';
end

% default title for figure
if ~exist('titleFig', 'var')
    switch typeEntity
        case 'file'
            titleFig = 'Select a file';
        case 'directory'
            titleFig = 'Select a directory';
    end
end

if isnumeric(fileNumbers)
    fileNumbers = num2str(fileNumbers);
end

% convert the title string to cell array
if ~iscell(titleFig)
    titleFig = {titleFig};
end

% If selection type is single selection mode will be 1
if strcmpi(typeSelection, 'multiple')
    selectMode = 2;
end

[modalityType] = icatb_get_modality;

if strcmpi(modalityType, 'fmri')
    helpLabel = 'GIFT-Help';
elseif strcmpi(modalityType, 'smri')
    helpLabel = 'SBM-Help';
else
    helpLabel = 'EEGIFT-Help';
end

%% Draw figure
handleFig = icatb_getGraphics(figName, 'normal', 'interactiveFileWindow', 'off');
set(handleFig, 'menubar', 'none');
set(handleFig, 'name', titleFig{1});
% add help menu to the file selection window
menu1H = uimenu('parent', handleFig, 'label', helpLabel);
% add
menu2H = uimenu(menu1H, 'label', 'Interactive Figure Window', 'callback', 'icatb_openHTMLHelpFile(''icatb_file_selection_window.htm'');');
% directory information
dirInfo.directories = []; dirInfo.prevFiles.file = [];

% Initialise the data with the files
filePrefixData = struct('numfiles', 0, 'files', [], 'dirInfo', dirInfo);

% set the figure data
figureData = struct('filePrefix', filePrefixData, 'currentDirectory', '', 'fileType', fileType);

set(handleFig, 'userdata', figureData);

%%%%%%%%%%%%% Draw Title here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% title color
titleColor = [0 0.9 0.9];

% fonts
titleFont = 13;

axisH = axes('Parent', handleFig, 'position', [0 0 1 1], 'visible', 'off');

xPos = 0.5; yPos = 0.97;

text(xPos, yPos, titleFig, 'color', titleColor, 'FontAngle', 'italic', 'fontweight', 'bold', 'fontsize', titleFont, 'HorizontalAlignment', ...
    'center', 'FontName', UI_FONTNAME, 'parent', axisH);

%%%%%%%%%%%%% end for drawing title %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% Draw User Interface Controls %%%%%%%%%%%%%%%%

%% Define Positions for the user interface controls
controlYOffset = 0.02; controlXOffset = 0.025;

% Browse Text box (shows the current directory)
browseTextPos = [controlXOffset yPos - 0.1 1 - 2*controlXOffset 0.055];

% Previous directories position
prevDirHeight = 0.05;
prevDirPos = [browseTextPos(1) browseTextPos(2) - browseTextPos(4) - controlYOffset browseTextPos(3) prevDirHeight];

% Define filter positions only if typeEntity is file
% Back button Position
backWidth = 0.12; backHeight = 0.05;
backPos = [browseTextPos(1) prevDirPos(2) - prevDirPos(4) - controlYOffset backWidth backHeight];

% Filter
filterPos = [backPos(1) + backPos(3) + controlXOffset backPos(2) 0.15 0.05];

% Type of Filter
editPos = [filterPos(1) + filterPos(3) + controlXOffset backPos(2) 0.25 0.05];


% File no. position
fileNoPos(1) = editPos(1) + editPos(3) + controlXOffset;
remSpace = 0.45*(1 - controlXOffset - fileNoPos(1));
fileNoPos = [fileNoPos(1) backPos(2) remSpace 0.05];

% Drives
% drivePos(1) = editPos(1) + editPos(3) + controlXOffset;
% drivePos = [drivePos(1) backPos(2) fileNoPos 0.05];

drivePos(1) = fileNoPos(1) + fileNoPos(3) + controlXOffset;
drivePos = [drivePos(1) backPos(2) remSpace 0.05];

% % Browse in a particular directory
browseTextH = icatb_uicontrol('parent', handleFig, 'units', 'normalized', 'style', 'edit', 'position', ...
    browseTextPos, 'string', outputDir, 'tag', 'TextForDirectory', 'horizontalalignment', 'left', ...
    'tooltipstring', 'Select the directory  you want to browse to...', 'callback', {@editboxCallback, handleFig});

% Draw History for directories
% show history of previous directories
[previousDir] = icatb_showcwdHistory;

prevDirH = icatb_uicontrol('parent', handleFig, 'units', 'normalized', 'style', 'popup', 'position', ...
    prevDirPos, 'string', str2mat('Directory History...', str2mat(previousDir)), 'tag', ...
    'directory history', 'tooltipstring', 'Directory History', 'userdata', str2mat(previousDir), ...
    'callback', {@doUpdatePopup, browseTextH, handleFig});

% Back button
backH = icatb_uicontrol('parent', handleFig, 'units', 'normalized', 'style', 'pushbutton', 'position', backPos, ...
    'string', 'Back', 'tooltipstring', 'Return to previous folder in the same drive...', 'callback', ...
    {@backbutton_Callback, browseTextH, handleFig});

if strcmpi(typeEntity, 'file')
    filterH = icatb_uicontrol('parent', handleFig, 'units', 'normalized', 'style', 'pushbutton', 'position', ...
        filterPos, 'string', 'Filter', 'tooltipstring', 'Leave blank or reset filter to *.* and press enter. Otherwise click on the filter pushbutton...', ...
        'callback', {@filterCallback, browseTextH, handleFig});
    editH = icatb_uicontrol('parent', handleFig, 'units', 'normalized', 'style', 'edit', 'position', editPos, ...
        'string', FilterText, 'horizontalalignment', 'left', 'tag', 'filter', 'tooltipstring', ...
        'Use filter to list files of a certain type...', 'callback', {@filterTextCallback, browseTextH, handleFig});
    fileNoH = icatb_uicontrol('parent', handleFig, 'units', 'normalized', 'style', 'edit', 'position', ...
        fileNoPos, 'string', fileNumbers, 'tooltipstring', 'Enter file numbers to include for nifti files.', ...
        'tag', 'file_numbers', 'callback', {@fileNoCallback, handleFig, browseTextH});
    
end

% Draw pop up bar for drawing the drives
driveStr = icatb_getdrives;

promptString = 'Drives';

driveH = icatb_uicontrol('parent', handleFig, 'units', 'normalized', 'style', 'popup', 'position', drivePos, ...
    'string', str2mat(promptString, driveStr{1:length(driveStr)}), 'tooltipstring', 'Drives in the system...', ...
    'tag', 'drives', 'userdata', str2mat(driveStr{1:length(driveStr)}), ...
    'callback', {@doUpdatePopup, browseTextH, handleFig});

% define push button positions here
okWidth = 0.2; okHeight = 0.05;

okPos = [0.75 - 0.5*okWidth 0.75*okHeight  okWidth okHeight];

refreshPos = okPos;

refreshPos(1) = 0.25 - 0.5*refreshPos(3);

% Origin of the listbox
listOrigin = okHeight + okPos(2) + 4*controlYOffset;


%%%%%%%%%%%%%%%%%%%%%%%%%%% Draw Listboxes Here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the position of the sub-folders listbox
subfolderPos = [browseTextPos(1) listOrigin 0.4 (backPos(2) - listOrigin - 2*controlYOffset)];

% Draw subfolders listbox here
subfolderH = icatb_uicontrol('parent', handleFig, 'units', 'normalized', 'style', 'listbox', 'position', ...
    subfolderPos, 'tag', 'folders', 'horizontalalignment', 'left', ...
    'userdata', outputDir, 'max', 1, 'min', 0, ...
    'callback', {@doUpdateTextbox_dir, browseTextH, handleFig});

% Entry Listbox Origin
entryOrigin = subfolderPos(1) + subfolderPos(3) + 2*controlXOffset;

% Define position of the expanded view of the files
entryPos = [entryOrigin listOrigin (1 - (entryOrigin + controlXOffset)) ...
    (backPos(2) - listOrigin - 2*controlYOffset)];

toolStr = 'Sub-folders in the current folder...';

% Tags are chaged depending upon whether the file or folder is to be selected
tag = 'folders_folder';

if strcmpi(typeEntity, 'file')
    toolStr = 'Files in the current folder...';
    tag = 'files_folder';
end

entryH = icatb_uicontrol('parent', handleFig, 'units', 'normalized', 'style', 'listbox', 'position', entryPos, 'tag', tag, 'BackgroundColor', BG2_COLOR,...
    'horizontalalignment', 'left', 'userdata', pwd, 'callback', {@filePrefixListbox, handleFig});

% Create a push button here
okH = icatb_uicontrol('parent', handleFig, 'units', 'normalized', 'style', 'pushbutton', 'position', okPos, 'string', 'OK', 'BackgroundColor', BUTTON_COLOR,...
    'tooltipstring', 'Select Entries...', 'callback', {@doOK, handleFig}, 'tag', 'doOk');

% Create a push button here
refreshH = icatb_uicontrol('parent', handleFig, 'units', 'normalized', 'style', 'pushbutton', 'position', ...
    refreshPos, 'string', 'Reset', 'tooltipstring', 'Refresh the drives and the folders or files...', ...
    'callback', {@refreshCallback, browseTextH, handleFig});

frameWidth = 0.4 + 2*controlXOffset;
fileTextHeight = okHeight + okPos(2) + 1.5*controlYOffset;
frameHeight = (listOrigin - fileTextHeight - 0.25*controlYOffset);
frameXOrigin = 0.5 - 0.5*frameWidth;
frameYOrigin = fileTextHeight - 0.5*controlYOffset;
framePos = [frameXOrigin frameYOrigin frameWidth frameHeight];

fileInfoTextPos(2) = framePos(2) + 0.005; fileInfoTextPos(4) = framePos(4) - 0.01;
fileInfoTextPos(1) = framePos(1) + 0.005; fileInfoTextPos(3) = framePos(3) - 0.01;

if strcmpi(typeEntity, 'file')
    buttonPos = fileInfoTextPos;
    buttonWidth = 0.12; buttonHeight = 0.05;
    buttonPos(1) = 0.9 - 0.5*buttonWidth ; buttonPos(3) = buttonWidth; buttonPos(4) = buttonHeight;
    buttonPos(2) = framePos(2) + 0.5*framePos(4) - 0.5*buttonPos(4);
    toolTip = 'Expand mode...';
    stringText = 'Expand';
    userData{1} = 'Prefix'; userData{2} = 'Expand';
    % Create a push button here
    buttonH = icatb_uicontrol('parent', handleFig, 'units', 'normalized', 'style', 'pushbutton', 'position', buttonPos, ...
        'string', stringText, 'tooltipstring', toolTip, 'callback', {@buttonHCallback, browseTextH, handleFig}, ...
        'tag', 'button_prefix_expand', 'userdata', userData);
    extentPos = get(buttonH, 'extent');
    buttonPos(3) = extentPos(3) + 0.001; buttonPos(4) = extentPos(4) + 0.0001;
    set(buttonH, 'position', buttonPos);
end

% set expand mode for single selection
if strcmpi(typeEntity, 'file') & strcmpi(typeSelection, 'single')
    set(buttonH, 'String', 'Expand');
    % execute the button callback
    buttonHCallback(buttonH, [], browseTextH, handleFig);
end


stringStaticText = 'Directories Selected:';
staticTag = 'directory-selected';

if strcmpi(typeEntity, 'file')
    stringStaticText = 'Files Selected:';
    staticTag = 'files-selected';
end

frametH = icatb_uicontrol('parent', handleFig, 'units', 'normalized', 'style', 'frame', 'position', framePos);
% Show the user number of files selected
staticTextH = icatb_uicontrol('parent', handleFig, 'units', 'normalized', 'style', 'text', 'position', fileInfoTextPos, 'BackgroundColor', BG2_COLOR,...
    'ForegroundColor', FONT_COLOR, 'fontunits', UI_FONTUNITS, 'fontname', UI_FONTNAME, 'FontSize', UI_FS - 1, 'string', stringStaticText, ...
    'horizontalalignment', 'left', 'tag', staticTag);

% Add select all button to select all files
if strcmpi(typeSelection, 'multiple') && strcmpi(typeEntity, 'file')
    
    select_all_pos = okPos;
    
    select_all_pos(1) = 0.5 - 0.5*select_all_pos(3);
    
    % Create a push button here
    select_all_H = icatb_uicontrol('parent', handleFig, 'units', 'normalized', 'style', 'pushbutton', 'position', select_all_pos, 'string', 'Select-all', ...
        'BackgroundColor', BUTTON_COLOR, 'ForegroundColor', BUTTON_FONT_COLOR, 'fontunits', UI_FONTUNITS, 'fontname', UI_FONTNAME, 'FontSize', UI_FS, ...
        'tooltipstring', 'Select_ALL', 'tag', 'select-all-button', 'callback', ...
        {@selectAllCallback, entryH, handleFig});
    
    
    % Get sub-folder listbox position
    sub_folderPos = get(subfolderH, 'position');
    
    % Plot Edit button to edit files
    editTextButtonWidth = 0.12; editTextButtonHeight = 0.04;
    editTextPos = [sub_folderPos(1), sub_folderPos(2) - controlYOffset - editTextButtonHeight, ...
        editTextButtonWidth, editTextButtonHeight];
    editTextH = icatb_uicontrol('parent', handleFig, 'units', 'normalized', 'style', 'pushbutton', 'position', ...
        editTextPos, 'string', 'Ed', 'tag', 'edit_files', 'tooltipstring', 'Edit files ...', ...
        'callback', {@editFilesCallback, handleFig});
end



%% Callback for the editbox
% Purpose: List entries in the listboxes
editboxCallback(browseTextH, [], handleFig);

% Selection mode: single or multiple selection
if strcmpi(typeEntity, 'file')
    if strcmpi(typeSelection, 'multiple')
        set(entryH, 'max', selectMode);
    end
end

%% set figure visibility on
try
    set(handleFig, 'visible','on');
    waitfor(handleFig);
catch
    if ishandle(handleFig)
        delete(handleFig);
    end
end

%% get the user data
if isappdata(0, 'listboxData')
    % get the application data
    ad = getappdata(0, 'listboxData');
    
    if strcmpi(typeEntity, 'directory')
        outputText = ad.directory;
    else
        numFiles = ad.numFiles;
        if numFiles == 0
            outputText = [];
        else
            outputText = ad.files;
        end
        filterAnswer = ad.filterAnswer;
    end
    % remove the application data
    rmappdata(0, 'listboxData');
    
    if count_datasets == numDataSets
        % change to old directory
        cd(outputDir);
    end
    
else
    
    % figure was deleted
    outputText = [];
    % change to old directory
    cd(outputDir);
    promptString = 'Figure window was quit';
    %error(promptString);
    disp(promptString);
end

drawnow;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % define function callbacks
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Edit box callback
function editboxCallback(handleObj, event_data, handles)
% Edit box callback
% purpose -
% 1. updates the strings in listbox
% 2. get the current directory information from the edit box
% 3. get the filter information from the filter edit box
% 4. List files with the file prefixes
% 5. If the click on the file prefix is identified then enable off and make
% file prefix window invisible, enable expanded view listbox on

% set the figure pointer watch
set(handles, 'pointer', 'watch');

try
    fd = get(handles, 'userdata');
    fileType = fd.fileType;
    clear fd;
    
    % get the current directory information from the edit box
    getEditString = get(handleObj, 'string');
    
    % Change directory
    cd(getEditString);
    
    % Initialise folder listbox
    % set the folder contents
    set(findobj(handles, 'tag', 'folders'), 'value', 1);
    
    % check if the object exists
    check_folder = isempty(findobj(handles, 'tag', 'folders_folder'));
    
    check_files = isempty(findobj(handles, 'tag', 'files_folder'));
    
    fileNoH = findobj(handles, 'tag', 'file_numbers');
    
    try
        fileNumbers = str2num(get(fileNoH, 'string'));
    catch
        fileNumbers = [];
    end
    
    % Initialising the folders and file listboxes to 1
    if check_folder == 0
        set(findobj(handles, 'tag', 'folders_folder'), 'value', 1);
    end
    
    if check_files == 0
        set(findobj(handles, 'tag', 'files_folder'), 'value', 1);
    end
    
    drawnow;
    
    % List directories
    tempDir = icatb_listDir(pwd);
    
    folderContents = []; % initialise vars
    
    % % All directories
    if ~isempty(tempDir)
        strIndex = strmatch('.', tempDir, 'exact');
        if isempty(strIndex)
            folderContents = str2mat('.', tempDir);
        else
            folderContents = tempDir;
        end
        clear strIndex;
    end
    
    clear tempDir;
    
    % set the folder contents
    set(findobj(handles, 'tag', 'folders'), 'string', folderContents);
    
    % get the filter information from the filter edit box
    filterTextString = get(findobj(handles, 'tag', 'filter'), 'string');
    
    stringVar = '*';
    
    if ~isempty(findobj(handles, 'tag', 'filter'))
        stringVar = filterTextString; % get the filter string
    end
    
    summaryContents = []; % Initialise contents of the right listbox
    
    if check_files == 0
        set(handles, 'pointer', 'watch');
        drawnow;
        % File prefix visibility and enable
        dirContents = icatb_listFiles_inDir(pwd, stringVar);
        if strcmpi(fileType, 'image')
            dirContents = icatb_rename_4d_file(dirContents, fileNumbers);
        end
        drawnow;
        listboxUserData.files = [];
        listboxUserData.indices_prefix = [];
        % Summary entries in a listbox: Classify files with the prefixes
        if ~isempty(dirContents)
            getUserData = get(findobj(handles, 'tag', 'button_prefix_expand'), 'userdata');
            % prefix mode
            if strcmpi(getUserData{1}, 'prefix')
                [summaryContents, indices_prefix] = icatb_classifyFiles(dirContents);
                listboxUserData.files = dirContents; % store the original files
                listboxUserData.indices_prefix = indices_prefix; % store the indices
            else
                % expand mode
                [summaryContents] = dirContents;
                listboxUserData.files = dirContents; % store the original files
            end
        end
        
        % Fix the selection highlight for R2008a
        myListH = findobj(handles, 'tag', 'files_folder');
        
        % Store the expanded view present directory
        set(myListH, 'userdata', listboxUserData);
        
        set(myListH, 'string', '');
        
        %% Set the listbox summary to the classified files
        set(myListH, 'string', summaryContents);
        
        drawnow;
        
        clear summaryContents;
        
    end
    
    %% If folders are to be displayed in the right side
    if check_folder == 0
        set(findobj(handles, 'tag', 'folders_folder'), 'string', folderContents);
    end
    
    %% Reset the figure pointer
    set(handles, 'pointer', 'arrow');
    
catch
    
    icatb_errorDialog(lasterr);
    set(handles, 'pointer', 'arrow');
    disp(lasterr);
end



%% popup callback
function doUpdatePopup(handleObj, event_data, handles, figHandle)
% Update the directory from the pop up
% set the text box with string from pop up

selectedValDrive = 0;

selectedValDrive = get(handleObj, 'value');

selectedValDrive = selectedValDrive - 1;

if selectedValDrive > 0
    driveString = get(handleObj, 'userdata');
    set(handles, 'string', deblank(driveString(selectedValDrive, :)));
    % Call editbox callback
    editboxCallback(handles, [], figHandle);
end

%% Sub folders callback
function doUpdateTextbox_dir(hObject, event_data, handles, figHandle)
% Update the : Get the string under selection and then update the
% filter string, and the contents in the left listbox for the sub folders

% Initialise the files listbox to 1
checkFileHandle = isempty(findobj(figHandle, 'tag', 'files_folder'));

% check Folders handles
checkFolderHandle = isempty(findobj(figHandle, 'tag', 'folders_folder'));

if checkFileHandle == 0
    set(findobj(figHandle, 'tag', 'files_folder'), 'value', 1);
end

if checkFolderHandle == 0
    set(findobj(figHandle, 'tag', 'folders_folder'), 'value', 1);
end

currentDir = deblank(get(handles, 'string')); % get the directory information from the browse text box

% Get the summary string
getString = get(hObject, 'string');

% get the selection under consideration
selection = get(hObject, 'value');

% get the string under consideration
%currentChar = fullfile(currentDir, deblank(getString(selection(length(selection)), :)));
currentChar = icatb_fullFile('directory', currentDir, 'files', deblank(getString(selection(length(selection)), :)));

cd(currentChar);

currentChar = deblank(pwd);

% set the filter the string selected
set(handles, 'string', currentChar);

editboxCallback(handles, [], figHandle);

%% Back button callback
function backbutton_Callback(handleObj, event_data, handles, figHandle)
% Back button used for browsing the directories in the same drive

% get the directory string
dirStr = deblank(get(handles, 'string'));

% Update the string
newStr = fullfile(dirStr, '..');

cd(newStr);

clear dirStr;

dirStr = deblank(pwd);

set(handles, 'string', dirStr);

editboxCallback(handles, [], figHandle);

% Filter button callback
function filterCallback(handleObj, event_data, handles, figHandle)
% reset filter to '*'
% update the list boxes

set(findobj(figHandle, 'tag', 'filter'), 'string', '*');

editboxCallback(handles, [], figHandle);

%% Filter Edit box Callback
function filterTextCallback(hObject, evd, handles, figHandle)

editboxCallback(handles, [], figHandle);

%% File Prefix Listbox Callback
function filePrefixListbox(handleObj, event_data, handles)
% gets the information from the selected strings in the right listbox and
% shows the number of files selected to the user

% if double click is encountered then don't update the files
if ~strcmpi(get(handles, 'selectionType'), 'open')
    % if folders listbox is detected
    check_folder = isempty(findobj(handles, 'tag', 'folders_folder'));
    % update the string in the listbox with the tag 'files_folder'
    check_files = isempty(findobj(handles, 'tag', 'files_folder'));
    
    % check select all
    checkSelectAll = isempty(findobj(handles, 'tag', 'select-all-button'));
    
    if checkSelectAll == 0
        if ~strcmpi(get(findobj(handles, 'tag', 'select-all-button'), 'enable'), 'on')
            set(findobj(handles, 'tag', 'select-all-button'), 'enable', 'on') ;
        end
    end
    
    if check_files == 0
        % Get the strings
        getString = get(handleObj, 'string');
        
        % Get the value under consideration
        getValue = get(handleObj, 'value');
        
        if checkSelectAll == 0
            % If the number of files selected is equal to the all files
            if size(getString, 1) == length(getValue)
                set(findobj(handles, 'tag', 'select-all-button'), 'enable', 'off') ;
            end
        end
        
        fileContents = [];
        % Store the vectors with different Indices
        % Get the figure user data
        figureData = get(handles, 'userdata');
        
        filePrefixData = figureData.filePrefix;
        
        numFiles = filePrefixData.numfiles;  % Number of files
        
        dirInfo = filePrefixData.dirInfo; % directory information
        
        temp = []; temp.files = [];
        
        getUserData = get(findobj(handles, 'tag', 'button_prefix_expand'), 'userdata');
        
        listboxUserData = get(handleObj, 'userdata'); % get the user data
        if ~isempty(listboxUserData.files)
            allFileNames = listboxUserData.files; % get the files
        end
        if ~isempty(listboxUserData.indices_prefix)
            indices_prefix = listboxUserData.indices_prefix; % get the prefix indices
        end
        
        % case for selecting all files
        if length(getValue) == size(getString, 1)
            fileContents = allFileNames;
            clear allFileNames;
            % selection mode is single
            if get(handleObj, 'max') == 1
                set(handleObj, 'enable', 'inactive');
                tempVariable = deblank(fileContents(1, :));
                clear fileContents;
                fileContents = tempVariable;
                clear tempVariable;
            end
            fullFileContents = icatb_fullFile('directory', pwd, 'files', fileContents); % full file contents
            if ~isempty(fullFileContents)
                % store the files
                if ~isempty(temp.files)
                    temp.files = str2mat(temp.files, fullFileContents);
                else
                    temp.files = fullFileContents;
                end
            end
            clear fileContents; clear fullFileContents;
        else
            % If prefix mode then list files
            if strcmpi(getUserData{1}, 'prefix')
                set(handles, 'pointer', 'watch')
                % Loop over each entry to get the full file path
                for ii = 1:length(getValue)
                    % List files with the prefixes
                    currentIndex = indices_prefix(getValue(ii), :); % current prefix indices
                    nonZeroIndices = find(currentIndex ~= 0); % non zero indices
                    fileContents = allFileNames(currentIndex(nonZeroIndices), :); % file contents
                    % check the select mode
                    if get(handleObj, 'max') == 1
                        set(handleObj, 'enable', 'inactive');
                        tempVariable = deblank(fileContents(1, :));
                        clear fileContents;
                        fileContents = tempVariable;
                        clear tempVariable;
                    end
                    
                    fullFileContents = icatb_fullFile('directory', pwd, 'files', fileContents);
                    if ~isempty(fullFileContents)
                        % store the files
                        if ~isempty(temp.files)
                            temp.files = str2mat(temp.files, fullFileContents);
                        else
                            temp.files = fullFileContents;
                        end
                    end
                    clear fileContents; clear fullFileContents;
                end
                set(handles, 'pointer', 'arrow')
            else
                % List files
                fileContents = icatb_fullFile('directory', pwd, 'files', allFileNames(getValue, :));
                % check the select mode
                if get(handleObj, 'max') == 1
                    set(handleObj, 'enable', 'inactive');
                    tempVariable = deblank(fileContents(1, :));
                    clear fileContents;
                    fileContents = tempVariable;
                    clear tempVariable;
                end
                
                if ~isempty(fileContents)
                    % store the files
                    if ~isempty(temp.files)
                        temp.files = str2mat(temp.files, fileContents);
                    else
                        temp.files = fileContents;
                    end
                end
                clear fileContents;
                
            end
            % End if condition for checking whether it is prefix mode or expand
            % mode
            
        end
        tempFiles = temp.files;
        clear temp;
        
        % directory Info contains the directories and the previousFiles
        if ~isempty(dirInfo.directories)
            numDirectories = size(dirInfo.directories, 1) + 1;
            dirInfo.directories = str2mat(dirInfo.directories, pwd);
            dirInfo.prevFiles(size(dirInfo.directories, 1)).file = tempFiles;
        else
            dirInfo.directories = pwd;
            dirInfo.prevFiles(1).file = tempFiles;
        end
        
        
        [newDirs, uniqueInd] = icatb_uniqueStr(cellstr(str2mat(dirInfo.directories)), 'dir');
        
        filePrefixData.files = str2mat(dirInfo.prevFiles(uniqueInd).file);
        
        filePrefixData.numfiles = size(filePrefixData.files, 1);
        
        % Report number of files
        checkStaticH = isempty(findobj(handles, 'tag', 'files-selected'));
        
        if checkStaticH == 0
            set(findobj(handles, 'tag', 'files-selected'), 'string', ['Files Selected: ', ...
                num2str(filePrefixData.numfiles)]);
        end
        
        filePrefixData.dirInfo = dirInfo;
        
        % Set the file prefix data to figure data
        figureData.filePrefix = filePrefixData;
        
        set(handles, 'userdata', figureData); % set the application data
        
    end
    
    if check_folder == 0
        
        figureData = get(handles, 'userdata');
        
        % Get the strings
        getString = get(handleObj, 'string');
        
        % Get the value under consideration
        getValue = get(handleObj, 'value');
        
        set(handleObj, 'enable', 'inactive');
        
        % Report number of files
        checkStaticH = isempty(findobj(handles, 'tag', 'directory-selected'));
        
        if checkStaticH == 0
            set(findobj(handles, 'tag', 'directory-selected'), 'string', ...
                ['Directories Selected: ', num2str(length(getValue))]);
        end
        
        % selected string
        selectedStr = deblank(getString(getValue, :));
        
        getTextboxString = deblank(get(findobj(handles, 'tag', 'TextForDirectory'), 'string'));
        
        getStr = fullfile(getTextboxString, selectedStr);
        
        cd(getStr);
        
        % Store the current directory
        figureData.currentDirectory = pwd;
        
        % set the list box userdata
        set(handles, 'userdata', figureData);
        
    end
    
end

%% Refresh Callback
function refreshCallback(handleObj, event_data, handles, figHandle)
% Refresh entries in text box and drives

% Initialise all the controls
set(figHandle, 'pointer', 'watch');

fd = get(figHandle, 'userdata');
fileType = fd.fileType;
clear fd;

% Initialise entries in listbox
filePrefixData.numfiles = 0; filePrefixData.files = [];

dirInfo.directories = []; dirInfo.prevFiles.file = [];

filePrefixData.dirInfo = dirInfo;

figureData.filePrefix = filePrefixData;

figureData.currentDirectory = '';

figureData.fileType = fileType;

set(figHandle, 'userdata', figureData);

% Draw pop up bar for drawing the drives
driveStr = icatb_getdrives;

promptString = 'Drives';

set(findobj(figHandle, 'tag', 'drives'), 'string', str2mat(promptString, driveStr{1:length(driveStr)}), 'userdata', ...
    str2mat(driveStr{1:length(driveStr)}));

check_folder = isempty(findobj(figHandle, 'tag', 'folders_folder'));

% update the string in the listbox with the tag 'files_folder'
check_files = isempty(findobj(figHandle, 'tag', 'files_folder'));

% check select all
check_select_all = isempty(findobj(figHandle, 'tag', 'select-all-button'));

if check_folder == 0
    set(findobj(figHandle, 'tag', 'folders_folder'), 'enable', 'on');
end

if check_files == 0
    set(findobj(figHandle, 'tag', 'files_folder'), 'enable', 'on');
end

if check_select_all == 0
    set(findobj(figHandle, 'tag', 'select-all-button'), 'enable', 'on');
end

% Call the edit box callback
editboxCallback(handles, [], figHandle);

% check if the static handle is present(report the number of files
% selected)
checkStaticH = isempty(findobj(figHandle, 'tag', 'files-selected'));

if checkStaticH == 0
    set(findobj(figHandle, 'tag', 'files-selected'), 'string', ['Files Selected: ', num2str(filePrefixData.numfiles)]);
end

% check if the static handle is present(report the number of files
% selected)
checkStaticH = isempty(findobj(figHandle, 'tag', 'directory-selected'));

if checkStaticH == 0
    set(findobj(figHandle, 'tag', 'directory-selected'), 'string', ['Directories Selected: 0']);
end

set(figHandle, 'pointer', 'arrow');

%% Select All callback
function selectAllCallback(handleObj, event_data, handles, figHandle)

% Find the selection mode
selectMode = get(handles,  'max');

% If selection mode is greater than 1
if selectMode > 1
    listString = get(handles, 'string');
    set(handles, 'value', 1:size(listString, 1));
    filePrefixListbox(handles,  [], figHandle);
    % Select all the entries
    set(handleObj, 'enable', 'off');
end

%% OK button callback
function doOK(ok_btn, event_data, handles)
% Callback for the ok push button
% Store the results in a structure

% Check the existence of the folders or files
check_folder = isempty(findobj(handles, 'tag', 'folders_folder'));

% update the string in the listbox with the tag 'files_folder'
check_files = isempty(findobj(handles, 'tag', 'files_folder'));

% store this in application data
if check_files == 0
    
    % Filter control
    filterControl = findobj(handles, 'tag', 'filter');
    
    % Get the filter answer
    filterAnswer = get(filterControl, 'string');
    
    % get the figure data
    figureData = get(handles, 'userdata');
    
    % file prefix files
    filePrefixData = figureData.filePrefix;
    
    filePrefixFiles = [];
    
    if ~isempty(filePrefixData.files)
        filePrefixFiles = filePrefixData.files;
    end
    
    numFiles = size(filePrefixFiles, 1);
    
    okData.files = filePrefixFiles;
    
    okData.numFiles = numFiles;
    
    okData.filterAnswer = filterAnswer;
    
    setappdata(0, 'listboxData', okData);
    
    set(handles, 'userdata', []);
    
    delete(handles);
    
end

% store this in application data
if check_folder == 0
    figureData = get(handles, 'userdata');
    okData.directory = deblank(figureData.currentDirectory);
    
    % force to select the directory
    if ~isempty(okData.directory)
        
        setappdata(0, 'listboxData', okData);
        delete(handles);
    else
        % show dialog box for selecting the directory
        icatb_errorDialog('Directory is not selected. Please select the directory in the right listbox.', ...
            'Select Directory', 'modal');
    end
    
end


function buttonHCallback(handleObj, event_data, handles, figHandle)
% Purpose: List files in expanded mode or prefix mode

getUserData = get(handleObj, 'userdata');
getString = get(handleObj, 'string');

temp = getUserData;
temp{1} = getString;
kk = 1;
% Comparing the userdata with the string of the button
% Swapping the userdata
for ii = 1:length(getUserData)
    if ~strcmpi(getUserData{ii}, getString)
        kk = kk + 1;
        temp{kk} = getUserData{ii};
    end
end

clear getUserData;
getUserData = temp;
clear temp;

if strcmpi(getString, 'expand')
    set(handleObj, 'userdata', getUserData);
    editboxCallback(handles, [], figHandle);
    set(handleObj, 'string', 'Prefix', 'tooltipstring', 'Prefix mode...');
else
    set(handleObj, 'userdata', getUserData);
    editboxCallback(handles, [], figHandle);
    set(handleObj, 'string', 'Expand', 'tooltipstring', 'Expand mode...');
end



%% EDIT Files Callback
function editFilesCallback(hObject, event_data, handles)
% Files can also be edited by typing the text

try
    
    figureData = get(handles, 'userdata');
    
    filePrefixData = figureData.filePrefix;
    
    filePrefixFiles = [];
    
    if ~isempty(filePrefixData.files)
        filePrefixFiles = filePrefixData.files;
    end
    
    set(handles, 'visible', 'off');
    
    if ~isempty(filePrefixFiles)
        editString = filePrefixFiles;
    else
        editString = '';
    end
    
    
    %%%% Plot Figure showing the selected text %%%%%%%%%%%%%
    newFigureHandle = icatb_getGraphics('Edit Files ...', 'normal', 'edit_files_figure');
    
    set(newFigureHandle, 'menubar', 'none');
    
    xOrigin = 0.05;
    yOrigin = 0.2;
    editWidth = 1 - 2*xOrigin;
    editHeight = 0.7;
    editPos = [xOrigin, yOrigin, editWidth, editHeight];
    
    editTextH = icatb_uicontrol('parent', newFigureHandle, 'units', 'normalized', 'style', 'edit', 'position', ...
        editPos, 'string', editString, 'tag', 'editbox_text', 'max', 2, 'min', 0, 'horizontalalignment', 'left');
    
    % Plot Ok and cancel buttons
    buttonWidth = 0.12; buttonHeight = 0.05;
    okPos = [0.75 - 0.5*buttonWidth, 0.1, buttonWidth, buttonHeight];
    cancelPos = okPos;
    cancelPos(1) = 0.25 - 0.5*buttonWidth;
    
    cancelH = icatb_uicontrol('parent', newFigureHandle, 'units', 'normalized', 'style', 'pushbutton', ...
        'position', cancelPos, 'string', 'Cancel', 'callback', {@setFilesCancelCallback, newFigureHandle, ...
        handles}, 'tooltipstring', 'Return to the file selection window ...');
    
    okH = icatb_uicontrol('parent', newFigureHandle, 'units', 'normalized', 'style', 'pushbutton', ...
        'position', okPos, 'string', 'OK', 'callback', {@setFilesOkCallback, newFigureHandle, ...
        handles, editTextH}, 'tooltipstring', 'Select the items in the text box and return the result ...');
    
catch
    set(handles, 'visible', 'on');
    icatb_displayErrorMsg;
end


%% SET Files Button Callback
function setFilesOkCallback(hObject, event_data, handles, mainFigHandle, editH)
% Edit files and return the result

try
    editString = get(editH, 'string');
    
    figureData = get(mainFigHandle, 'userdata');
    
    filePrefixData = figureData.filePrefix;
    
    files_in_folder = findobj(mainFigHandle, 'tag', 'files_folder');
    
    files = cellstr(editString);
    
    files = str2mat(deblank(cellstr(strvcat(files))));
    
    if isempty(files)
        files = [];
    end
    
    filePrefixData.files = files;
    
    filePrefixData.numFiles = size(files, 1);
    
    figureData.filePrefix  = filePrefixData;
    
    set(mainFigHandle, 'userdata', figureData);
    
    msgString = ['Number of selected files is ', num2str(filePrefixData.numFiles), ...
        '. Do you want to continue?'];
    
    [answerQuestion] = icatb_questionDialog('title', 'Number of files', 'textbody', msgString);
    
    okHandle = findobj(mainFigHandle, 'tag', 'doOk');
    
    if answerQuestion
        delete(handles);
        % Execute ok button callback
        doOK(okHandle, [], mainFigHandle);
    end
    
catch
    set(mainFigHandle, 'visible', 'on');
    icatb_displayErrorMsg;
end


%% SET Files Cancel Callback
function setFilesCancelCallback(hObject, event_data, handles, mainFigHandle)
% Cancel the edit box figure and return to the main figure

delete(handles);
set(mainFigHandle, 'visible', 'on');


%% File No. Callback
function fileNoCallback(hObject, event_data, handles, editH)
% File number callback

editboxCallback(editH, [], handles);