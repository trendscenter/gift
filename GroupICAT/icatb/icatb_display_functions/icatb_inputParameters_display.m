function parameters = icatb_inputParameters_display(displayType)
% figure window used to pass the necessary parameters for displaying
% components with the corresponding display type

if ~exist('displayType', 'var')
    displayType = 'component explorer';
end

% run defaults file to get global variables
icatb_defaults;
global PARAMETER_INFO_MAT_FILE; %Holds information for group session parameters
global FUNCTIONAL_DATA_FILTER;

% Screen Color Defaults
global BG_COLOR;
global BG2_COLOR;
global BUTTON_COLOR;
global BUTTON_FONT_COLOR;
global FONT_COLOR;
global AXES_COLOR;
global FONT_COLOR;

% Fonts
global UI_FONTNAME;
global UI_FONTUNITS;
global UI_FS;

% display_defaults
global SORT_COMPONENTS;
global IMAGE_VALUES;
global THRESHOLD_VALUE;
global CONVERT_Z;
global IMAGES_PER_FIGURE;
global ANATOMICAL_PLANE;

displayType = lower(displayType);

%%%%%%%%% Load structural Image %%%%%%%%%%
% select the structural image: components will be overlaid on the structural
P = icatb_selectEntry('typeEntity', 'file', 'title', 'Select Structural Image', 'filter', '*.img;*.nii', ...
    'typeselection', 'single', 'fileType', 'image', 'filenumbers', 1);

P = icatb_parseExtn(P);

[startPath, f_name, extn] = fileparts(P);

% check the image extension
if ~strcmpi(extn, '.nii') & ~strcmpi(extn, '.img')
    error('Structural image should be in NIFTI or Analyze format');
end

if isempty(startPath)
    startPath = pwd;
end

% structural image file
structuralFile = deblank(P);

% get the data
[figureData] = getCompFiles(displayType, startPath);

% get the necessary vars
compFiles = figureData.compFiles;
filesOutputDir = figureData.outputDir;
numComp = figureData.numComp;
compNumbers = figureData.compNumbers;

%--Image Values
options(1).str = 'Positive and Negative';
options(2).str = 'Positive';
options(3).str = 'Absolute Value';
options(4).str = 'Negative'; % Added negative image values

% Apply defaults from icatb_defaults
imageOptions = {'Positive and Negative', 'Positive', 'Absolute Value', 'Negative'};
options(1).str =  IMAGE_VALUES;

% Apply defaults from icatb_defaults
imageOptions = checkDefaults_gui(options(1).str, imageOptions, 'exact');

for jj = 1:length(imageOptions)
    options(jj).str = imageOptions{jj};
end

% needed for all the visualization methods
numParameters = 1;
inputText(numParameters).promptString = 'Image Values';
inputText(numParameters).answerString = str2mat(options.str);
inputText(numParameters).uiType = 'popup';
inputText(numParameters).answerType = 'string';
inputText(numParameters).tag = 'image_values';
inputText(numParameters).enable = 'on';
inputText(numParameters).tooltipString = 'Only use voxels.';
inputText(numParameters).userData = '';
inputText(numParameters).uiPos = [0.4 0.05];
clear options;


%--Convert To Z Scores
% Apply defaults from icatb_defaults
zOptions = {'No', 'Yes'};
options(1).str =  CONVERT_Z;
% Apply defaults from icatb_defaults
zOptions = checkDefaults_gui(options(1).str, zOptions, 'exact');
for jj = 1:length(zOptions)
    options(jj).str = zOptions{jj};
end

numParameters = numParameters + 1;
% used for all the visualization methods
inputText(numParameters).promptString = 'Convert To Z Scores';
inputText(numParameters).answerString = str2mat(options.str);
inputText(numParameters).uiType = 'popup';
inputText(numParameters).answerType = 'string';
inputText(numParameters).tag = 'Convert_To_ZScores';
inputText(numParameters).enable = 'on';
inputText(numParameters).tooltipString = 'Convert Images to Z score.';
inputText(numParameters).userData = '';
inputText(numParameters).uiPos = [0.25 0.05];
clear options;

numParameters = numParameters + 1;
inputText(numParameters).promptString = 'Threshold Value';
inputText(numParameters).answerString = THRESHOLD_VALUE;
inputText(numParameters).uiType = 'edit';
inputText(numParameters).answerType = 'numeric';
inputText(numParameters).tag = 'Threshold_Value';
inputText(numParameters).enable = 'on';
inputText(numParameters).tooltipString = 'Value to use as threshold';
inputText(numParameters).userData = '';
inputText(numParameters).uiPos = [0.25 0.05];

% options for number of images per figure
otherOptions = {'1', '4', '9', '16', '25'};
options(1).str =  IMAGES_PER_FIGURE;
% Apply defaults from icatb_defaults
otherOptions = checkDefaults_gui(options(1).str, otherOptions, 'exact');

for jj = 1:length(otherOptions)
    options(jj).str = otherOptions{jj};
end

numParameters = numParameters + 1;
inputText(numParameters).promptString = 'Images Per Figure';
inputText(numParameters).answerString = str2mat(options.str);
inputText(numParameters).uiType = 'popup';
inputText(numParameters).answerType = 'numeric';
inputText(numParameters).tag = 'images_per_figure';
inputText(numParameters).enable = 'on';
inputText(numParameters).tooltipString = 'The number of spatial maps display on each figure';
inputText(numParameters).userData = '';
inputText(numParameters).uiPos = [0.25 0.05];
clear options;

slicePlane = lower(ANATOMICAL_PLANE);
imagVol = icatb_returnHInfo(structuralFile);
% get the slices in mm for the corresponding plane
[sliceParameters] = icatb_get_slice_def(imagVol, slicePlane);
% get the slices in mm
slices_in_mm = sliceParameters.slices;
clear sliceParameters;
% construct string
slices_in_mm = icatb_constructString(slices_in_mm);
options(1).str = slices_in_mm;

numParameters = numParameters + 1;
inputText(numParameters).promptString = 'Slice Range';
inputText(numParameters).answerString = options.str;
inputText(numParameters).uiType = 'edit';
inputText(numParameters).answerType = 'numeric';
inputText(numParameters).tag = 'slice_range';
inputText(numParameters).enable = 'on';
inputText(numParameters).tooltipString = ['Enter the slices you wish to look at.' ...
    'You can enter a vector of slices, ex [1 5 10] or [1:skip:34].'];
inputText(numParameters).userData = '';
inputText(numParameters).uiPos = [0.4 0.05];
clear options;

planeOptions = {'Axial', 'Sagittal', 'Coronal'};
planeOptions = checkDefaults_gui(lower(ANATOMICAL_PLANE), lower(planeOptions), 'optional');

for jj = 1:length(planeOptions)
    options(jj).str = planeOptions{jj};
end

numParameters = numParameters + 1;
inputText(numParameters).promptString = 'Anatomical Plane';
inputText(numParameters).answerString = str2mat(options.str);
inputText(numParameters).uiType = 'popup';
inputText(numParameters).answerType = 'string';
inputText(numParameters).tag = 'anatomical_plane';
inputText(numParameters).enable = 'on';

%%%%%%%%% user data %%%%%%%%%%%%%%%
%anatomicalUserData.structHInfo = structHInfo;
anatomicalUserData.targetTag = 'slice_range';
inputText(numParameters).userData = anatomicalUserData;
inputText(numParameters).tooltipString = ['Enter the slices you wish to look at.', ...
    ' You can enter a vector of slices, ex [1 5 10] or [1:skip:34].'];
inputText(numParameters).uiPos = [0.4 0.05];
%%%%%%%% done with defining the parameter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% select the corresponding parameters only
switch lower(displayType)
    case 'component explorer'
        tagsNeeded = {'image_values', 'Convert_To_ZScores', 'Threshold_Value', 'images_per_figure', ...
            'anatomical_plane', 'slice_range'};
        %paramVec = (1:1:6);
    case 'orthogonal viewer'
        tagsNeeded = {'image_values', 'Convert_To_ZScores', 'Threshold_Value'};
        %paramVec = [1 2 3];
    case 'composite viewer'
        tagsNeeded = {'image_values', 'Convert_To_ZScores', 'Threshold_Value', ...
            'anatomical_plane', 'slice_range'};

        %paramVec = [1 2 3 6 5];
end

% Initialise vars
paramVec = zeros(1, length(tagsNeeded));
for ii = 1:length(tagsNeeded)
    paramVec(ii) = strmatch(tagsNeeded{ii}, str2mat(inputText.tag), 'exact');
end

numUIcontrols = length(paramVec);

figureData.pwd = pwd;
figureData.inputText = inputText(paramVec);
figureData.displayType = displayType;
figureData.structFile = structuralFile;

% Setup figure for GUI
[InputHandle] = icatb_getGraphics('Display Options', 'normal', 'figure');

% Set no menu bar for the figure
set(InputHandle, 'Menubar', 'None', 'tag', 'disp_parameters_gui', 'userdata', figureData);

menu1H = uimenu('parent', InputHandle, 'label', 'GIFT-Help'); % create a help menu
callbackStr = 'icatb_openHTMLHelpFile(''icatb_helpManual.htm'');';
menuName = 'GIFT';
switch lower(displayType)
    case 'component explorer'
        callbackStr = 'icatb_openHTMLHelpFile(''icatb_component_explorer.htm'');';
        menuName = 'Component Explorer';
    case 'composite viewer'
        callbackStr = 'icatb_openHTMLHelpFile(''icatb_composite_viewer.htm'');';
        menuName = 'Composite Viewer';
    case 'orthogonal viewer'
        callbackStr = 'icatb_openHTMLHelpFile(''icatb_ortho_viewer.htm'');';
        menuName = 'Orthogonal Viewer';
end

menu2H = uimenu(menu1H, 'label', menuName, 'callback', callbackStr);

[InputHandle] = icatb_plot_controls_fig(inputText, get(InputHandle, 'tag'), 'on', 'Done', 'Cancel', ...
    cellstr(str2mat(inputText(paramVec).tag)), InputHandle);

cancelHandle = findobj(InputHandle, 'tag', 'Cancel');
okHandle = findobj(InputHandle, 'tag', 'Done');

%%%%% set function callbacks %%%%%

% Cancel callback
set(cancelHandle, 'callback', {@closeCallback, InputHandle});

% Done callback
set(okHandle, 'callback', {@applyCallback, InputHandle});

% anatomical plane Callback
set(findobj(InputHandle, 'tag', 'anatomical_plane'), 'callback', {@anatomicalCallback, InputHandle});

try
    set(InputHandle, 'visible', 'on');
    waitfor(InputHandle);
catch
    if ishandle(InputHandle)
        delete(InputHandle);
    end
end

% get the parameters from the figure
if isappdata(0, 'inputParaDisplay')
    parameters = getappdata(0, 'inputParaDisplay');
    % remove the application data
    rmappdata(0, 'inputParaDisplay');
    % load the image files here
    %numComp = size(compFiles, 1);
    parameters.numComp = numComp;
    % structural file
    parameters.structFile = structuralFile;
    %parameters.structHInfo = structHInfo;
    parameters.compFiles = compFiles;
    parameters.filesOutputDir = filesOutputDir;
    parameters.compNumbers = compNumbers;
end

function optionsString = checkDefaults_gui(choiceString, optionsString, stringChar)
% put the defaults at the Top

temp = optionsString; % Assign options string to a temporary variable

if strcmp(stringChar, 'exact')
    matchIndex = strmatch(lower(choiceString), lower(temp), stringChar); % find the index of the choice string
else
    matchIndex = strmatch(lower(choiceString), lower(temp)); % find the index of the choice string
end

jj = 1; optionsString{1} = temp{matchIndex}; % choice string

for numOptions = 1:length(optionsString)
    if numOptions ~= matchIndex
        jj = jj + 1;
        optionsString{jj} = temp{numOptions}; % collect other options below the choice string
    end
end

function anatomicalCallback(hObject, evd, handles)
% anatomical plane callback

% purpose: updates the editbox slice plane
getString = lower(get(hObject, 'string'));
getValue = get(hObject, 'value');
slicePlane = deblank(getString(getValue, :));

figureData = get(handles, 'userdata');
% get the name of the structural file
structFile = figureData.structFile;
imagVol = icatb_returnHInfo(structFile);
% get the slices in mm for the corresponding plane
[sliceParameters] = icatb_get_slice_def(imagVol, slicePlane);
% get the slices in mm
slices_in_mm = sliceParameters.slices;
clear sliceParameters;
% construct string
slices_in_mm = icatb_constructString(slices_in_mm);

% set the editbox to the display the default size of slice plane
set(findobj(handles, 'tag', 'slice_range'), 'string', slices_in_mm);

function applyCallback(handleObj, event_data, handles)
% get the required parameters

figureData = get(handles, 'userdata');
% get the input parameters
inputText = figureData.inputText;

parameters = figureData; %struct;
% loop over number of input args
for ii = 1:length(inputText)
    answerType = inputText(ii).answerType;
    answerTag = [inputText(ii).tag];
    objHandle = findobj(handles, 'tag', answerTag);
    if strcmp(lower(inputText(ii).uiType), 'edit')
        % check the answer type
        if strcmp(lower(answerType), 'numeric')
            answerVal = str2num(get(objHandle, 'string'));
        else
            answerVal = get(objHandle, 'string');
        end
    elseif strcmp(lower(inputText(ii).uiType), 'popup')
        % check the answer type
        getString = get(objHandle, 'string');
        getVal = get(objHandle, 'value');
        answerVal = deblank(getString(getVal, :));
        % convert to numeric
        if strcmp(lower(answerType), 'numeric')
            answerVal = str2num(answerVal);
        end
    end
    % set field to the parameters
    parameters = setfield(parameters, strrep(deblank(lower(inputText(ii).tag)), '_', ''), answerVal);
end

% get the return Value
if strcmp(lower(parameters.imagevalues), 'positive')
    returnValue = 2;
elseif strcmp(lower(parameters.imagevalues), 'absolute value')
    returnValue = 3;
elseif strcmp(lower(parameters.imagevalues), 'negative')
    returnValue = 4;
else
    returnValue = 1;
end

if strcmp(lower(parameters.converttozscores), 'yes')
    convertToZ = 1;
else
    convertToZ = 0;
end

% pass the return values and convert to z options
parameters.returnValue = returnValue;
parameters.convertToZ = convertToZ;

% set application data
setappdata(0, 'inputParaDisplay', parameters);

delete(handles);

% close callback
function closeCallback(hObject, event_data, handles)

delete(handles);

function [figureData] = getCompFiles(displayType, startPath)
% get the component files and the component numbers
%  load the ICA Time course and get the number of components in case of
%  nifti and zip files


try

    icatb_defaults;
    % component naming
    global COMPONENT_NAMING;
    files_in_zip = {};
    zipFileName = {};

    msgString = 'Are the component images compressed in a zip file?';

    % Check if the files are zipped or not
    [answerZipFile] = icatb_questionDialog('title', msgString, 'textbody', msgString);

    % unzip component files to temporary directory
    if answerZipFile == 1
        % select zip file
        compZipFile = icatb_selectEntry('typeEntity', 'file', 'title', 'Select a component zip file', ...
            'filter', '*.zip', 'typeselection', 'single', 'startPath', startPath);
        drawnow;
        % zip file path and name
        [zipFilePath, zipFileName, extn] = fileparts(compZipFile);
        zipFileName = compZipFile;

        cd(zipFilePath);

        tempDirName = ['temp_comp_zip'];

        helpString = ['Unzipping ', zipFileName, ' ...'];
        % unzip the files to current location
        disp(helpString);
        helpHandle = helpdlg('Unzipping component file. Please wait ...', 'Unzipping files ...');

        % make a file pattern for deleting files inside temporary directory
        if strcmp(zipFilePath(end), filesep)
            deleteFilePattern = [zipFilePath, tempDirName, filesep, '*'];
            checkDirName = [zipFilePath, tempDirName];
        else
            deleteFilePattern = [zipFilePath, filesep, tempDirName, filesep, '*'];
            checkDirName = [zipFilePath, filesep, tempDirName];
        end

        % unzip zip file to a temporary directory get the file contents
        if ~exist(checkDirName, 'dir')
            mkdir(zipFilePath, tempDirName);
        end
        icatb_unzip(zipFileName, checkDirName);
        files_in_zip = icatb_listFiles_inDir(checkDirName, '*');
        %files_in_zip = icatb_fullFile('directory', zipFilePath, 'files', files_in_zip);
        cd(zipFilePath);

        % check if files are present in the zip file
        if isempty(files_in_zip)
            rmdir(tempDirName);
            error(['No files present in zip file :', zipFileName]);
        end

        % make a cell array
        files_in_zip = cellstr(files_in_zip);


        % delete all the files inside that temporary directory
        delete(deleteFilePattern);

        % remove the temporary directory
        rmdir(tempDirName);

        icatb_unzip(zipFileName, zipFilePath);

        if ishandle(helpHandle)
            delete(helpHandle);
        end
        % change start path
        startPath = zipFilePath;
    end

    cname = regexprep(COMPONENT_NAMING, '^(_)', '');
    % component file pattern
    filePattern = ['*', cname, '*.img', ';', '*', cname, '*.nii'];

    %%%%%% Load component files %%%%%%%%%%%
    % use the options for the particular display method
    switch displayType
        case 'component explorer'
            selectionMode = 'multiple'; % selection mode
            titleStr = 'Select component map files'; % title
            % select the component files
            compFiles = icatb_selectEntry('typeEntity', 'file', 'title', titleStr, ...
                'filter', filePattern, 'typeselection', selectionMode, 'startPath', startPath);
        case 'composite viewer'
            % select a maximum of five
            selectionMode = 'multiple'; % selection mode
            titleStr = 'Select component map files'; % file
            compFiles = icatb_selectEntry('typeEntity', 'file', 'title', titleStr, ...
                'filter', filePattern, 'typeselection', selectionMode, 'startPath', startPath);
        case 'orthogonal viewer'
            % fmri data
            fmriFiles = icatb_selectEntry('typeEntity', 'file', 'title', 'Select original functional files', ...
                'filter', '*.img;*.nii', 'typeselection', 'multiple', 'startPath', startPath, 'fileType', 'image');
            if isempty(fmriFiles)
                error('Functional data is not selected');
            end
            figureData.fmriFiles = fmriFiles;

            selectionMode = 'single'; % selection mode
            titleStr = 'Select a component'; % title
            % select component
            compFiles = icatb_selectEntry('typeEntity', 'file', 'title', titleStr, ...
                'filter', filePattern, 'typeselection', selectionMode, 'startPath', startPath);
        otherwise
            error('Please check the display type');
    end

    % check if component files are selected or not
    if isempty(compFiles)
        error('Component files are not selected');
    end

    % Loop over components
    for ii = 1:size(compFiles, 1)
        % temporary file
        temp = deblank(compFiles(ii, :));
        temp = icatb_parseExtn(temp);
        [pathstr, fName, extn] = fileparts(temp);
        imData(ii).extns = extn; %store extensions
        imData(ii).file_name = fName; % file name
        imData(ii).pathstr = pathstr; % directory

        if ~strcmpi(imData(ii).extns, '.nii') &  ~strcmpi(imData(ii).extns, '.img')
            error('Images should be in Nifti or analyze format');
        end

        % if more than one file is selected
        if ii > 1
            % files should have same extensions
            if ~strcmpi(imData(ii).extns, imData(ii-1).extns)
                error('Component files should be of same extensions');
            end

            % files should be in the same directory
            if ~strcmpi(imData(ii).pathstr, imData(ii-1).pathstr)
                error(['Component files should be in the same directory']);
            end

            % can't have multiple files for NII and zip
            if strcmpi(imData(ii).extns, '.nii')
                % select the first file for zip and Nii
                compFiles = deblank(compFiles(1, :));
                disp(['Can''t select more than one file for files of extension ', imData(ii).extns, ...
                    '. Selecting component file: ', compFiles, ' to display.']);
                break;
            end
        end
        % end for checking
    end
    % end for loop over componet files

    % handle nifti data here
    if strcmpi(imData(1).extns, '.nii')
        % load time course here to get the number of components
        [A] = icatb_loadICATimeCourse(compFiles);
        % number of components
        numComp = size(A, 2);
        % Initialise variables
        temp = repmat(struct('name', []), 1, numComp);
        for ii = 1:numComp
            temp(ii).name = ['Component ', num2str(ii)];
        end
        answerString = str2mat(temp.name); % form component strings
        clear temp;
        [answerQuestion3, name_button] = icatb_listdlg('PromptString', titleStr, 'SelectionMode', selectionMode,...
            'ListString', str2mat(answerString), 'movegui', 'center', 'windowStyle', 'modal');
        if isempty(answerQuestion3)
            error('Component/components is/are not selected');
        end
        % component numbers
        compNumbers = answerQuestion3;

    else
        % component numbers
        for ii = 1:size(compFiles, 1)
            temp = imData(ii).file_name;
            underScorePos = icatb_findstr(temp, '_');

            if ii == 1
                oldCompName = temp(1:underScorePos(end));
            end
            currentCompName = temp(1:underScorePos(end));

            % check viewing set
            if ~strcmpi(oldCompName, currentCompName)
                error(['Components of the same viewing set must be selected where first component naming is ', oldCompName, ...
                    ' and ', 'file number ', num2str(ii), ' naming is ', currentCompName]);
            end
            % end for checking

            % get the component number
            compNumbers(ii) = str2num(temp(underScorePos(end)+1:end));
        end
        clear compFiles;
        compFiles = icatb_listFiles_inDir(pathstr, ['*', temp(1:underScorePos(end)), '*.img']);
        compFiles = icatb_fullFile('files', compFiles, 'directory', pathstr);
    end
    % end for checking the cases


    % set the fields to the data
    if strcmpi(displayType, 'composite viewer')
        if length(compNumbers) > 5
            disp('Atmost five components with different colorbars are allowed ...');
            compNumbers = compNumbers(1:5);
        end
    end
    % check composite viewer

    % set the required fields to the figure data
    figureData.compNumbers = compNumbers;
    figureData.numComp = length(compNumbers);
    figureData.compFiles = compFiles;
    figureData.outputDir = pathstr;
    figureData.zipFileName = zipFileName;
    figureData.files_in_zip = files_in_zip;

catch
    if exist('helpHandle', 'var')
        if ishandle(helpHandle)
            delete(helpHandle);
        end
    end
    % delete the unzipped files
    if ~isempty(files_in_zip)
        icatb_delete_file_pattern(files_in_zip, zipFilePath);
    end
    icatb_displayErrorMsg;
end
