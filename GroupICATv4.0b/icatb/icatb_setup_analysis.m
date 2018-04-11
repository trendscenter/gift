function varargout = icatb_setup_analysis(inputFile)
%% Setup ICA analysis
%
% Inputs are:
% inputFile - Input file for batch analysis.
%
% Outputs:
% ICA Parameter file.

if exist('inputFile', 'var')
    [pathstr, fileN, extn] = fileparts(inputFile);
    if ~strcmpi(extn, '.mat')
        varargout{1} = icatb_read_batch_file(inputFile);
        return;
    else
        sub_file = inputFile;
    end
end

% Define input parameters
inputText = icatb_define_parameters;

[modalityType] = icatb_get_modality;

matchedIndex = strmatch('prefix', str2mat(inputText.tag), 'exact');
output_prefix = inputText(matchedIndex).answerString;

% Tags to be plotted in main figure window
if strcmpi(modalityType, 'fmri')
    tagsTobePlotted = {'prefix', 'files', 'estimate_components', 'numComp', 'autofill', 'algorithm', 'which_analysis', 'parallel_info'};
elseif strcmpi(modalityType, 'eeg')
    tagsTobePlotted = {'prefix', 'files', 'numComp', 'autofill', 'algorithm', 'which_analysis', 'parallel_info'};
else
    tagsTobePlotted = {'prefix', 'files', 'estimate_components', 'numComp', 'algorithm', 'which_analysis'};
end

if strcmpi(modalityType, 'fmri')
    helpLabel = 'GIFT-Help';
    figTitle = 'GIFT Setup ICA GUI';
    htmlFile = 'icatb_setup_ica.htm';
elseif strcmpi(modalityType, 'smri')
    helpLabel = 'SBM-Help';
    figTitle = 'SBM Setup ICA GUI';
    htmlFile = 'icatb_setup_ica.htm';
else
    helpLabel = 'EEGIFT-Help';
    figTitle = 'EEGIFT Setup ICA GUI';
    htmlFile = 'icatb_eeg_setup_ica.htm';
end


if ~exist('sub_file', 'var')
    
    % Select analysis output directory
    [oldDir] = icatb_selectEntry('typeEntity', 'directory', 'title', 'Select Analysis Output Directory');
    
else
    % set the prefix automatically
    formStr = 'Subject.mat';
    [oldDir, fileNN, file_extn] = fileparts(sub_file);
    if (isempty(oldDir))
        oldDir = pwd;
    end
    inputString = [fileNN, file_extn];
    matchIndex = icatb_findstr(inputString, formStr);
    output_prefix = inputString(1:matchIndex(1)-1);
    clear fileNN; clear file_extn;
    clear matchIndex; clear inputString;
end

handle_visibility = 'on';

drawnow;

% put a check for selecting the output directory
if isempty(oldDir)
    error('Output directory is not selected');
end

% Change working directory
cd(oldDir);

%%%%%%%%% Draw figure and set up the controls %%%%%%%%

figureTag = 'setup_ICAGUI';

figHandle = findobj('tag', figureTag);

if ~isempty(figHandle)
    delete(figHandle);
end

% Setup figure for GUI
[InputHandle] = icatb_getGraphics(figTitle, 'normal', figureTag, handle_visibility);

% store the directory information to sesInfo
sesInfo.userInput.pwd = oldDir;

handles_data.sesInfo = sesInfo;

% set figure data
set(InputHandle, 'userdata', handles_data, 'Menubar', 'none');

% set up menus (on click display an input dialog box)
setupDefaults = uimenu('parent', InputHandle, 'label', 'SetupICA-Defaults', 'tag', 'SetupICA-Defaults');

% help on setup ICA
giftHelpTitle = uimenu('parent', InputHandle, 'label', helpLabel);
giftHelpMenu = uimenu(giftHelpTitle, 'label', 'Setup-ICA', 'callback', ...
    ['icatb_openHTMLHelpFile(''', htmlFile, ''');']);

algoIndex = strmatch('algoType', char(inputText.tag), 'exact');

inputText(algoIndex).tag = 'algorithm';


countTags_menu = 0;

% total tags in menu
tagsVec = zeros(1, length(inputText) - length(tagsTobePlotted));

% store the tags that need to be plotted in a menu
for ii = 1:length(inputText)
    strIndex = strmatch(inputText(ii).tag, tagsTobePlotted, 'exact');
    if isempty(strIndex)
        countTags_menu = countTags_menu + 1;
        tagsInMenu{countTags_menu} = inputText(ii).tag;
        tagsVec(countTags_menu) = ii;
    end
end

% plot the input parameters to the parameter figure
icatb_plotInputPara('input_prefix', inputText, 'controls_to_plot', tagsTobePlotted, 'handles', ...
    InputHandle);

%%%%%%%%% End for drawing the figure and setting up the controls %%%%


%%%%%%%%%% Set up function callbacks %%%%%%%%%%%%%%%%%

% setup defaults data
setupDefaultsData.tagsTobePlotted = tagsInMenu;  setupDefaultsData.inputText = inputText;
setupDefaultsData.tagsVec = tagsVec;

%% Open figure window when SetupICA-Defaults menu is clicked
set(setupDefaults, 'callback', {@setupDefaultsCallback, InputHandle}, 'userdata', setupDefaultsData);

%% Answer function callbacks
set(findobj(InputHandle, 'tag', 'prefix'), 'callback', {@editboxCallback, InputHandle});

%% Load the functional files (data callback)
set(findobj(InputHandle, 'tag', 'files'), 'callback', ...
    {@dataCallback, InputHandle, setupDefaults});

%% Set callback for cancel button
set(findobj(InputHandle, 'tag', 'Cancel'), 'callback', {@closeCallback, InputHandle});

estimateCompHandle = findobj(InputHandle, 'tag', 'estimate_components');

if ~isempty(estimateCompHandle)
    % set callback for estimate components callback
    set(estimateCompHandle, 'callback', {@estimateCompCallback, InputHandle});
end

componentHandle = findobj(InputHandle, 'tag', 'numComp');
autofillHandle = findobj(InputHandle, 'tag', 'autofill');

%% Component callback
set(componentHandle, 'callback', {@setCompCallback, InputHandle, autofillHandle});

%% Autofill handle
set(autofillHandle, 'callback', {@autoFillCallback, InputHandle, componentHandle});

%% Done callback
okHandle = findobj(InputHandle, 'tag', 'Done');
set(okHandle, 'callback', {@applyCallback, InputHandle, setupDefaults});

%% Which Analysis Callback
WhichAnalysisH = findobj(InputHandle, 'tag', 'which_analysis');
set(WhichAnalysisH, 'callback', {@WhichAnalysisCallback, InputHandle});

%% ICA Algorithm callback
set(findobj(InputHandle, 'tag', 'algorithm'), 'callback', {@icaTypeCallback, InputHandle});

% %% Group ica type callback
% groupICATypeTag = 'group_ica_type';
% groupICAAnalysisHandle = findobj(InputHandle, 'tag', groupICATypeTag);
% if (~isempty(groupICAAnalysisHandle))
%     set(groupICAAnalysisHandle, 'callback', {@group_ica_typeCallback, InputHandle});
% end

%% Parallel analysis callback
parallelH = findobj(InputHandle, 'tag', 'parallel_info');
if (~isempty(parallelH))
    set(parallelH, 'callback', {@parallelCallback, InputHandle});
end

%%%%%%%%%% End for setting up the function callbacks %%%%%%%%%%%%%%%%%

if exist('output_prefix', 'var')
    prefixH = findobj(InputHandle, 'tag', 'prefix');
    set(prefixH, 'string', output_prefix);
end



% Execute this function for GUI
editboxCallback(findobj(InputHandle, 'tag', 'prefix'), [], InputHandle);


%%%%%%%%% Sub Functions %%%%%%%%%

function  varargout = icatb_plotInputPara(varargin)
% sub function to plot the input parameters

tagsTobePlotted = {};
figureTag = 'Input Parameters';
windowStyle = 'normal';
varargout = {};

% loop over arguments
for ii = 1:2:nargin
    if strcmpi(varargin{ii}, 'input_prefix')
        inputText = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'controls_to_plot')
        tagsTobePlotted = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'handles')
        InputHandle = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'tag_figure')
        figureTag = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'windowstyle')
        windowStyle = varargin{ii + 1};
    end
    % end for checking the input args
end
% end for loop

if ~exist('inputText', 'var')
    error('inputText variable does not exist');
end

% tags to plot
if isempty(tagsTobePlotted)
    disp('No controls to plot');
    return;
end

if ~exist('InputHandle', 'var')
    % Setup figure for GUI
    [InputHandle] = icatb_getGraphics(figureTag, 'normal', figureTag);
    
    % Set no menu bar for the figure
    set(InputHandle, 'Menubar', 'none');
end
% end for checking the input handle

if ~ispc
    windowStyle = 'normal';
end

set(InputHandle, 'windowStyle', windowStyle);

[InputHandle] = icatb_plot_controls_fig(inputText, figureTag, 'on', 'Done', 'Cancel', tagsTobePlotted, ...
    InputHandle);


function set_data_para(handles, menuH)
% transfer figure data to menu

% sesInfo structure
handles_data = get(handles, 'userdata');
sesInfo = handles_data.sesInfo;

% menu data
menuData = get(menuH, 'userdata');

% input text structure
inputText = menuData.inputText;

numOfSub = 0; numOfSess = 0;

if isfield(sesInfo.userInput, 'numOfSub')
    % Number of subjects and sessions
    numOfSub = sesInfo.userInput.numOfSub;
end

if isfield(sesInfo.userInput, 'numOfSess')
    % Number of subjects and sessions
    numOfSess = sesInfo.userInput.numOfSess;
end

algoH = findobj(handles, 'tag', 'algorithm');
icaStr = get(algoH, 'string');
icaVal = get(algoH, 'value');

%% group ica type
groupICATypeTag = 'group_ica_type';
group_icaIndex = strmatch(groupICATypeTag, cellstr(char(inputText.tag)), 'exact');
useTemporalICA = 0;
if (~isempty(group_icaIndex))
    groupICATypeOpts = lower(cellstr(inputText(group_icaIndex).answerString));
    groupICATypeVal = inputText(group_icaIndex).value;
    useTemporalICA = strcmpi('temporal', groupICATypeOpts{groupICATypeVal});
end


backReconTag = 'backReconType';
backReconIndex = strmatch(backReconTag, cellstr(char(inputText.tag)), 'exact');
if (~isempty(backReconIndex))
    inputText(backReconIndex).enable = 'on';
    if (strcmpi(deblank(icaStr(icaVal, :)), 'iva-gl') || strcmpi(deblank(icaStr(icaVal, :)), 'iva-l') ...
            || strcmpi(deblank(icaStr(icaVal, :)), 'moo-icar') || strcmpi(deblank(icaStr(icaVal, :)), 'gig-ica') || strcmpi(deblank(icaStr(icaVal, :)), 'constrained ica (spatial)'))
        inputText(backReconIndex).enable = 'inactive';
    end
end

if numOfSub*numOfSess > 0
    % get the number of components from the textbox in main figure window
    numCompH = findobj(handles, 'tag', 'numComp');
    compIndex = strmatch('numComp', char(inputText.tag), 'exact');
    % set component number to sesInfo
    sesInfo.userInput.numComp = str2num(deblank(get(numCompH, 'string')));
    % set the component number to the input text
    inputText(compIndex).answerString = num2str(sesInfo.userInput.numComp);
    
    % get the data reduction tag
    dataRedIndex = strmatch('numReductionSteps', char(inputText.tag), 'exact');
    
    inputText(dataRedIndex).enable = 'off';
    
    numReductionSteps = 1;
    if (~strcmpi(deblank(icaStr(icaVal, :)), 'iva-gl') && ~strcmpi(deblank(icaStr(icaVal, :)), 'iva-l') && ~strcmpi(deblank(icaStr(icaVal, :)), 'moo-icar') ...
            && ~strcmpi(deblank(icaStr(icaVal, :)), 'gig-ica') && ~strcmpi(deblank(icaStr(icaVal, :)), 'constrained ica (spatial)'))
        % Number of data reduction steps
        if(numOfSub == 1 && numOfSess == 1)
            numReductionSteps = 1;
        else
            if (~useTemporalICA)
                inputText(dataRedIndex).enable = 'on';
            end
            % get the value from the data reduction control
            getStr = inputText(dataRedIndex).answerString; getVal = inputText(dataRedIndex).value;
            if isfield(sesInfo.userInput, 'numReductionSteps')
                numReductionSteps = sesInfo.userInput.numReductionSteps;
                if (numReductionSteps == 3)
                    numReductionSteps = 2;
                end
                getVal = strmatch(num2str(numReductionSteps), getStr, 'exact');
                % set the value to the data reduction tag
                inputText(dataRedIndex).value = getVal;
            else
                numReductionSteps = str2num(deblank(getStr(getVal, :)));
            end
        end
        % end for checking data reduction steps
    end
    
    % set reduction steps to sesInfo
    sesInfo.userInput.numReductionSteps = numReductionSteps;
    
    prefix1 = 'numOfPC1'; prefix2 = 'numOfPC2';
    
    % get the index of the PC1, PC2
    PC1Index = strmatch(prefix1, char(inputText.tag), 'exact');
    PC2Index = strmatch(prefix2, char(inputText.tag), 'exact');
    
    PCBefore = sesInfo.userInput.numComp;
    
    [minTp, minTpInd] = min(sesInfo.userInput.diffTimePoints);
    
    if (minTp == 1)
        error('Error:TimePoints', 'Please re-select the data as the no of selected files is found to be 1 (%s)\n', deblank(sesInfo.userInput.files(minTpInd).name(1, :)));
    end
    
    PCBefore = round(min([minTp, 1.5*PCBefore]));
    
    if numReductionSteps == 1
        %% PC1
        inputText(PC1Index).promptString = 'Number Of PC/IC (Step 1)';
        inputText(PC1Index).enable = 'inactive';
        sesInfo.userInput.numOfPC1 = sesInfo.userInput.numComp;
        inputText(PC1Index).answerString = num2str(sesInfo.userInput.numOfPC1);
        
        %% PC2
        inputText(PC2Index).promptString = 'Number Of PC (Step 2)';
        inputText(PC2Index).enable = 'inactive';
        sesInfo.userInput.numOfPC2 = 0;
        inputText(PC2Index).answerString = num2str(sesInfo.userInput.numOfPC2);
        
    elseif numReductionSteps == 2
        % display the object with PC1 and set the string to PC/IC
        %% PC1
        inputText(PC1Index).promptString = 'Number Of PC (Step 1)';
        inputText(PC1Index).enable = 'on';
        if ~isfield(sesInfo.userInput, 'numOfPC1')
            
            sesInfo.userInput.numOfPC1 = PCBefore; %sesInfo.userInput.numComp;
            
        else
            if isempty(sesInfo.userInput.numOfPC1)
                sesInfo.userInput.numOfPC1 = PCBefore; %sesInfo.userInput.numComp;
            end
        end
        inputText(PC1Index).answerString = num2str(sesInfo.userInput.numOfPC1);
        
        %% PC2
        inputText(PC2Index).promptString = 'Number Of PC/IC (Step 2)';
        inputText(PC2Index).enable = 'inactive';
        
        sesInfo.userInput.numOfPC2 = sesInfo.userInput.numComp;
        inputText(PC2Index).answerString = num2str(sesInfo.userInput.numOfPC2);
        
    end
    
    menuData.inputText = inputText;
    
    % set the updated data
    set(menuH, 'userdata', menuData);
    
    handles_data.sesInfo = sesInfo;
    set(handles, 'userdata', handles_data);
    
end
% end for checking the data

function [varOut] = check_var_integer(currentVar)
% check whether the variable in the string is integer or not

[varOut] = icatb_check_variable_integer(currentVar);


function [varOut, status, message] = check_var_character(currentVar)
% check whether the variable is a valid character or not

[varOut] = icatb_check_char(currentVar);
[status, message] = icatb_errorCheck(currentVar, 'output_prefix', 'prefix');


%%%%%%%%%% Function Callbacks %%%%%%%%

function editboxCallback(handleObj, event_data, handles)
% callback for output prefix text box
% Execute this function when prefix is entered

icatb_defaults;

global PARAMETER_INFO_MAT_FILE; %Holds information for group session parameters

% get the string for the text box
getString = get(handleObj, 'string'); %(output file prefixes)

if ~isempty(getString)
    
    % check if it is a valid prefix
    [getString, status, message] = check_var_character(getString);
    
    if status == 0
        error(message);
    end
    
end

% Get the figure data
handles_data = get(handles, 'userdata');

% sesInfo structure
sesInfo = handles_data.sesInfo;

% directory where the information will be stored
oldDir = sesInfo.userInput.pwd;

% get the prefix of the output string
sesInfo.userInput.prefix = getString;

% All the entered parameters are going to be stored in this file
sesInfo.userInput.param_file = [sesInfo.userInput.prefix, PARAMETER_INFO_MAT_FILE, '.mat'];
sesInfo.userInput.param_file = fullfile(oldDir, sesInfo.userInput.param_file);

% estimateComponents callback
set(findobj(handles, 'tag', 'estimate_components'), 'enable', 'off');

% number of components callback
set(findobj(handles, 'tag', 'numComp'), 'enable', 'inactive');

% subject file
subjectFile = [getString, 'Subject.mat']; % check if the subjects file exist or not

% full file path for the subject file
subjectFile = fullfile(oldDir, subjectFile);

% find the setup ICA defaults menu
menuH = findobj(handles, 'tag', 'SetupICA-Defaults');

handles_data.sesInfo = sesInfo;

% set the figure data
set(handles, 'userdata', handles_data);

% check if the subjects file exists
if exist(subjectFile)
    set(findobj(handles, 'tag', 'files'), 'style', 'popup', 'string', ...
        char('Yes', 'No'));
    % data callback
    dataCallback(findobj(handles, 'tag', 'files'), [], handles, menuH);
else
    set(findobj(handles, 'tag', 'files'), 'style', 'pushbutton', 'string', ...
        'Select');
end


% Data callback for selecting the functional files
function dataCallback(handleObj, event_data, handles, menuH)

try
    
    % handles data
    handles_data = get(handles, 'userdata');
    
    sesInfo = handles_data.sesInfo;
    
    oldDir = sesInfo.userInput.pwd;
    
    icatb_defaults;
    
    global PARAMETER_INFO_MAT_FILE; %Holds information for group session parameters
    global FUNCTIONAL_DATA_FILTER;
    
    % Screen Color Defaults
    global BG_COLOR;
    global FONT_COLOR;
    global AXES_COLOR;
    
    % disable the estimate components and the data reduction parameters
    set(findobj(handles, 'tag', 'estimate_components'), 'enable', 'off');
    
    set(findobj(handles, 'tag', 'numReductionSteps'), 'enable', 'off');
    
    % number of components callback
    set(findobj(handles, 'tag', 'numComp'), 'enable', 'inactive');
    
    % get the prefix of the output string
    sesInfo.userInput.prefix = get(findobj(handles, 'tag', 'prefix'), 'string');
    
    % get the subject matrix file
    subjectFile = [sesInfo.userInput.prefix, 'Subject.mat'];
    subjectFile = fullfile(oldDir, subjectFile);
    
    % All the entered parameters are going to be stored in this file
    sesInfo.userInput.param_file = [sesInfo.userInput.prefix, PARAMETER_INFO_MAT_FILE, '.mat'];
    sesInfo.userInput.param_file = fullfile(oldDir, sesInfo.userInput.param_file);
    
    newParamFile = sesInfo.userInput.param_file; % get the new parameter file
    
    % get the style of the UIcontrol
    getStyle = get(handleObj, 'style');
    
    if strcmp(lower(getStyle), 'pushbutton')
        getSubjects = 'no';
    else
        getPopValue = get(handleObj, 'value');
        getPopString = get(handleObj, 'string');
        getSubjects = lower(deblank(getPopString(getPopValue, :)));
    end
    
    % First determine whether the subjects file exists or not
    % If it doesn''t exist then show the user the error dialog
    % else load the subject file
    if strcmp(getSubjects, 'yes') %getSubjects == 2
        if ~exist(subjectFile, 'file')
            icatb_errorDialog('Subject File doesn''t exist. Please use the other option', 'File not found');
        else
            load(subjectFile);
        end
    end
    %%%%%% Check if the files doesn't exist %%%%%%%%%%%%%%%%%
    
    
    if ~exist('modalityType', 'var')
        [modalityType] = icatb_get_modality;
    else
        [modalityType2] = icatb_get_modality;
        if ~strcmpi(modalityType, modalityType2)
            delete(handles);
            if strcmpi(modalityType, 'fmri')
                disp('fMRI subject file already exists. Use GIFT toolbox to set up analysis');
            elseif strcmpi(modalityType, 'smri')
                disp('sMRI subject file already exists. Use SBM toolbox to set up analysis');
            else
                disp('EEG subject file already exists. Use EEGIFT toolbox to set up analysis');
            end
            return;
        end
    end
    
    % store these fields regarding dataType, complex naming
    dataType = 'real'; read_complex_images = 'real&imaginary'; write_complex_images = 'real&imaginary';
    
    sesInfo.userInput.dataType = lower(dataType);
    sesInfo.userInput.read_complex_images = lower(read_complex_images);
    sesInfo.userInput.write_complex_images = lower(write_complex_images);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    if ~strcmp(getSubjects, 'yes')
        
        [sesInfo] = icatb_name_complex_images(sesInfo, 'read');
        
        set(handles, 'pointer', 'watch');
        
        [files, designMatrix, numOfSub, numOfSess, dataSelMethod, diffTimePoints, spmMatFlag] = icatb_dataSelection(...
            [], sesInfo.userInput.pwd, sesInfo.userInput.prefix, ...
            sesInfo.userInput.read_complex_file_naming, sesInfo.userInput.read_complex_images);
        setappdata(0, 'create_mask_gica',1);
        drawnow;
        sesInfo.userInput.files = files;
        SPMFiles = designMatrix;
        sesInfo.userInput.dataSelMethod = dataSelMethod;
        sesInfo.userInput.designMatrix = designMatrix;
        sesInfo.userInput.spmMatFlag = spmMatFlag;
        sesInfo.userInput.diffTimePoints = diffTimePoints;
        
        icatb_save(subjectFile, 'files', 'numOfSub', 'numOfSess', 'SPMFiles', 'modalityType');
        
        %%%%%%%%%% Adjust the autofill popup box %%%%%%%%%%%
        % set the autofill popup to yes
        autoFillH = findobj(handles, 'tag', 'autofill');
        % auto fill string
        autoFillString = get(autoFillH, 'string');
        
        % auto fill value
        matchIndex = strmatch('yes', lower(autoFillString), 'exact');
        
        
        % set the autofill popup to yes
        set(autoFillH, 'value', matchIndex);
        %%%%%%%%%% End for adjusting the autofill popup box %%%%%%%%%%%
        
        % Create mask
        %[sesInfo] = create_mask(sesInfo);
        
        set(handles, 'pointer', 'arrow');
        
    end
    % end for selecting the data
    
    
    % Initialise all the variables
    sesInfo.userInput.modality = modalityType;
    sesInfo.userInput.files = files;
    sesInfo.userInput.designMatrix = SPMFiles;
    sesInfo.userInput.numOfSub = numOfSub;
    sesInfo.userInput.numOfSess = numOfSess;
    numOfDataSets = numOfSub * numOfSess;
    sesInfo.userInput.numOfGroups1 = numOfDataSets;
    
    % If the subjects data exist
    % Update the files and the numbers for the data reduction steps
    if strcmp(getSubjects, 'yes') %getSubjects == 2
        % Update the parameter field
        if exist(newParamFile, 'file')
            % load the parameter file
            load(newParamFile);
            % make sure to update the new parameter file
            sesInfo.userInput.param_file = newParamFile;
        end
        
        % check time points
        if ~isfield(sesInfo.userInput, 'diffTimePoints')
            set(handles, 'pointer', 'watch');
            if (strcmpi(modalityType, 'fmri') || strcmpi(modalityType, 'smri'))
                % get the count for time points
                diffTimePoints = icatb_get_countTimePoints(files);
            else
                diffTimePoints = icatb_get_num_electrodes(files);
            end
            % end for getting time points
            sesInfo.userInput.diffTimePoints = diffTimePoints;
            set(handles, 'pointer', 'arrow');
        else
            diffTimePoints = sesInfo.userInput.diffTimePoints;
        end
        % end for checking the time points
    end
    
    % number of scans
    numberOfScans = diffTimePoints(1);
    
    if ~isfield(sesInfo.userInput, 'numComp')
        numComp = 20;
    else
        numComp = sesInfo.userInput.numComp;
    end
    
    % can't extract more components than time points
    if min(diffTimePoints) < numComp
        numComp = min(diffTimePoints);
    end
    % set the number of components
    sesInfo.userInput.numComp = numComp;
    
    
    algoHandle = findobj(handles, 'tag', 'algorithm');
    tmpAlgVal = get(algoHandle, 'value');
    tmpAlgStr = cellstr(get(algoHandle, 'string'));
    
    numCompH = findobj(handles, 'tag', 'numComp');
    set(numCompH, 'enable', 'on', 'string', num2str(sesInfo.userInput.numComp));
    
    if (strcmpi(tmpAlgStr{tmpAlgVal}, 'moo-icar') || strcmpi(tmpAlgStr{tmpAlgVal}, 'gig-ica') || strcmpi(tmpAlgStr{tmpAlgVal}, 'constrained ica (spatial)'))
        set(numCompH, 'enable', 'inactive');
    end
    
    %% make sure these two vars exist
    sesInfo.userInput.pwd = oldDir;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    handles_data.sesInfo = sesInfo;
    
    % set the user data
    set(handles, 'userdata', handles_data);
    
    clear inputText;
    
    menuData = get(menuH, 'userdata');
    inputText = menuData.inputText;
    
    if strcmp(getSubjects, 'yes')
        
        %% ICA algorithm
        algoTag = 'algorithm';
        if isfield(sesInfo.userInput, algoTag)
            algoHandle = findobj(handles, 'tag', algoTag);
            algoVal = getfield(sesInfo.userInput, algoTag);
            algoStr = cellstr(get(algoHandle, 'string'));
            set(algoHandle, 'value', algoVal);
            if (strcmpi(algoStr{algoVal}, 'iva-gl') || strcmpi(algoStr{algoVal}, 'iva-l') || strcmpi(algoStr{algoVal}, 'moo-icar')  || ....
                    strcmpi(algoStr{algoVal}, 'gig-ica') || strcmpi(algoStr{algoVal}, 'constrained ica (spatial)'))
                icaTypeCallback(algoHandle, [], handles);
            end
        end
        
        
        %% stability Type
        whichAnalysisTag = 'which_analysis';
        if isfield(sesInfo.userInput, whichAnalysisTag)
            whichAnalysisHandle = findobj(handles, 'tag', whichAnalysisTag);
            whichAnalysisVal = getfield(sesInfo.userInput, whichAnalysisTag);
            set(whichAnalysisHandle, 'value', whichAnalysisVal);
        end
        
        %% group ica type
        groupICATypeTag = 'group_ica_type';
        group_icaIndex = strmatch(groupICATypeTag, cellstr(char(inputText.tag)), 'exact');
        if ~isempty(group_icaIndex)
            if (isfield(sesInfo.userInput, groupICATypeTag))
                group_ica_type = getfield(sesInfo.userInput, groupICATypeTag);
            else
                group_ica_type = 'spatial';
            end
            group_ica_val = strmatch(group_ica_type, lower(cellstr(inputText(group_icaIndex).answerString)), 'exact');
            inputText(group_icaIndex).value = group_ica_val;
        end
        
        %         groupICAAnalysisHandle = findobj(handles, 'tag', groupICATypeTag);
        %
        %         if (~isempty(groupICAAnalysisHandle))
        %             if isfield(sesInfo.userInput, groupICATypeTag)
        %                 groupICATypeOpts = cellstr(get(groupICAAnalysisHandle, 'string'));
        %                 tmpgroupICAVal = 'spatial';
        %                 try
        %                     tmpgroupICAVal = lower(getfield(sesInfo.userInput, groupICATypeTag));
        %                 catch
        %                 end
        %
        %                 if (~isnumeric(tmpgroupICAVal))
        %                     tmpgroupICAVal = strmatch(tmpgroupICAVal, lower(groupICATypeOpts), 'exact');
        %                 end
        %                 set(groupICAAnalysisHandle, 'value', tmpgroupICAVal);
        %             end
        %         end
        
        %% Run Parallel?
        parallelTag = 'parallel_info';
        if isfield(sesInfo.userInput, parallelTag)
            parallelH = findobj(handles, 'tag', parallelTag);
            opts = {'serial', 'parallel'};
            parallelMode = 'serial';
            try
                parallelMode = sesInfo.userInput.parallel_info.mode;
            catch
            end
            parallelVal  = strmatch(parallelMode, opts, 'exact');
            %parallelVal = lower(getfield(sesInfo.userInput, parallelTag));
            set(parallelH, 'value', parallelVal);
        end
        
        
        
        %% Back-reconstruction
        backReconTag = 'backReconType';
        if (isfield(sesInfo.userInput, backReconTag))
            backReconType = getfield(sesInfo.userInput, backReconTag);
        else
            backReconType = 'regular';
        end
        
        if (strcmpi(backReconType, 'str'))
            backReconType = 'spatial-temporal regression';
        end
        
        backReconType = lower(backReconType);
        
        backReconIndex = strmatch(backReconTag, cellstr(char(inputText.tag)), 'exact');
        if ~isempty(backReconIndex)
            backReconVal = strmatch(backReconType, lower(inputText(backReconIndex).answerString), 'exact');
            if (~isempty(backReconVal))
                inputText(backReconIndex).value = backReconVal;
            else
                backReconStrings = cellstr(inputText(backReconIndex).answerString);
                sesInfo.userInput.backReconType = lower(deblank(backReconStrings{inputText(backReconIndex).value}));
                handles_data.sesInfo = sesInfo;
                % set the user data
                set(handles, 'userdata', handles_data);
            end
        end
        
        %% Data Pre-processing
        preprocTag = 'preproc_type';
        if (isfield(sesInfo.userInput, preprocTag))
            preprocType = getfield(sesInfo.userInput, preprocTag);
        else
            preprocType = 'remove mean per timepoint';
        end
        
        preprocType = lower(preprocType);
        
        preprocIndex = strmatch(preprocTag, cellstr(char(inputText.tag)), 'exact');
        if ~isempty(preprocIndex)
            
            preprocVal = strmatch(preprocType, lower(inputText(preprocIndex).answerString), 'exact');
            inputText(preprocIndex).value = preprocVal;
        end
        
        %% Group PCA Type
        groupPCATag = 'group_pca_type';
        if (isfield(sesInfo.userInput, groupPCATag))
            group_pca_type = getfield(sesInfo.userInput, groupPCATag);
        else
            group_pca_type = 'subject specific';
        end
        group_pca_type = lower(group_pca_type);
        group_pcaIndex = strmatch(groupPCATag, cellstr(char(inputText.tag)), 'exact');
        if ~isempty(group_pcaIndex)
            group_pca_val = strmatch(group_pca_type, lower(inputText(group_pcaIndex).answerString), 'exact');
            inputText(group_pcaIndex).value = group_pca_val;
        end
        
        %% PCA Type
        pcaTag = 'pcaType';
        if (isfield(sesInfo.userInput, pcaTag))
            pcaType = getfield(sesInfo.userInput, pcaTag);
        else
            pcaType = 'standard';
        end
        
        pcaType = lower(pcaType);
        
        pcaIndex = strmatch(pcaTag, cellstr(char(inputText.tag)), 'exact');
        if ~isempty(pcaIndex)
            
            pcaVal = strmatch(pcaType, lower(inputText(pcaIndex).answerString), 'exact');
            inputText(pcaIndex).value = pcaVal;
        end
        
        if isfield(sesInfo.userInput, 'convertToZ')
            sesInfo.userInput = rmfield(sesInfo.userInput, 'convertToZ');
        end
        
        % check the calibration step
        if isfield(sesInfo.userInput, 'scaleType') %& isfield(sesInfo.userInput, 'convertToZ')
            adjustVal = sesInfo.userInput.scaleType + 1;
            scale_opts = icatb_scaleICA;
            scaleStr = deblank(scale_opts(adjustVal, :));
            % get the adjust information
            adjustTag = 'scaleType';
            adjustIndex = strmatch(adjustTag, char(inputText.tag), 'exact');
            adjustVal = strmatch(lower(scaleStr), lower(inputText(adjustIndex).answerString), 'exact');
            if (isempty(adjustVal))
                adjustVal = 1;
            end
            inputText(adjustIndex).value = adjustVal;
        end
        
        % set mask information also
        maskTag = 'maskFile';
        if isfield(sesInfo.userInput, maskTag)
            maskVal = getfield(sesInfo.userInput, maskTag);
            if ~isempty(maskVal)
                maskVal = 2;
            else
                maskVal = 1;
            end
            maskIndex = strmatch(maskTag, char(inputText.tag), 'exact');
            if ~isempty(maskIndex)
                inputText(maskIndex).value = maskVal;
            end
        end
        
    end
    
    % set input text to menu
    menuData.inputText = inputText;
    % set menu data
    set(menuH, 'userdata', menuData);
    
    
    % set data reduction parameters
    set_data_para(handles, menuH);
    
    set(findobj(handles, 'tag', 'files'), 'style', 'popup', 'string', ...
        char('Yes', 'No'));
    
    set(findobj(handles, 'tag', 'files'), 'value', 1);
    
    % disable the estimate components and the data reduction parameters
    set(findobj(handles, 'tag', 'estimate_components'), 'enable', 'on');
    
catch
    set(handles, 'pointer', 'arrow');
    % Prints line number information along with the error message
    icatb_displayErrorMsg;
    % change to old directory
    %cd(oldDir);
end


function setupDefaultsCallback(hObject, event_data, handles)
% setup ICA defaults callback

% transfer data from sesInfo to menu
set_data_para(handles, hObject);


% get the user data for the menu
menuData = get(hObject, 'userdata');

% get the figure data
handles_data = get(handles, 'userdata');
sesInfo = handles_data.sesInfo;
inputText = menuData.inputText;


handles_visibility = 'on';


% tags to be plotted
tagsTobePlotted = menuData.tagsTobePlotted;
% tag vector
tagVec = menuData.tagsVec;

figureTag = 'SetupDefaults';
% Setup figure for GUI
InputH = icatb_getGraphics(figureTag, 'normal', figureTag, handles_visibility);

% Set no menu bar for the figure
set(InputH, 'Menubar', 'none');

modalityType = icatb_get_modality;

if strcmpi(modalityType, 'fmri')
    helpLabel = 'GIFT-Help';
    htmlFile = 'icatb_setup_ica_defaults.htm';
elseif strcmpi(modalityType, 'smri')
    helpLabel = 'SBM-Help';
    htmlFile = 'icatb_setup_ica_defaults.htm';
else
    helpLabel = 'EEGIFT-Help';
    htmlFile = 'icatb_eeg_setup_ica_defaults.htm';
end

% help on setup ICA
giftHelpTitle = uimenu('parent', InputH, 'label', helpLabel);
giftHelpMenu = uimenu(giftHelpTitle, 'label', 'Setup ICA Defaults', 'callback', ...
    ['icatb_openHTMLHelpFile(''', htmlFile, ''');']);

% plot the input parameters to the parameter figure
icatb_plotInputPara('input_prefix', inputText, 'controls_to_plot', tagsTobePlotted, ...
    'windowStyle', 'normal', 'handles', InputH);

dataRedControl = findobj(InputH, 'tag', 'numReductionSteps');

% Pre-processing data callback
preprocH = findobj(InputH, 'tag', 'preproc_type');
set(preprocH, 'callback', {@preprocCallback, handles});

% Back reconstruction callback
backReconH = findobj(InputH, 'tag', 'backReconType');
set(backReconH, 'callback', {@backReconCallback, handles});

% Group PCA Callback
groupPCAH = findobj(InputH, 'tag', 'pcaType');
set(groupPCAH, 'callback', {@pcaTypeCallback, handles});


% Group ICA type callback
groupICATypeTag = 'group_ica_type';
groupICAAnalysisHandle = findobj(InputH, 'tag', groupICATypeTag);
if (~isempty(groupICAAnalysisHandle))
    set(groupICAAnalysisHandle, 'callback', {@group_ica_typeCallback, handles});
end

% Covariance options callback
%set(findobj(InputH, 'tag', 'covariance_opts'), 'callback', {@covOptsCallback, handles, InputH});

% set the callback
set(dataRedControl, 'callback', {@dataRedCallback, InputH, handles, hObject});
set(findobj(InputH, 'tag', 'Done'), 'callback', {@doneCb, InputH, hObject, handles});
set(findobj(InputH, 'tag', 'Cancel'), 'callback', {@closeCallback, InputH});

% put a callback for selecting mask
maskTag = findobj(InputH, 'tag', 'maskFile');
if ~isempty(maskTag)
    set(maskTag, 'callback', {@selectMaskCallback, handles, hObject});
end

try
    set(InputH, 'visible', 'on');
    waitfor(InputH);
catch
    delete(InputH);
end


function doneCb(hObject, event_data, handles, menuH, mainHandle)
% get the menu data and transfer to main figure data (sesInfo)

% menu data
menuData = get(menuH, 'userdata');

% session information
handles_data = get(mainHandle, 'userdata');
sesInfo = handles_data.sesInfo;

tagVec = menuData.tagsVec;
inputText = menuData.inputText;

% set menu data to input text
for ii = 1:length(tagVec)
    currentTag = findobj(handles, 'tag', inputText(tagVec(ii)).tag);
    inputText(tagVec(ii)).value = get(currentTag, 'value');
    promptTag = findobj(handles, 'tag', ['prompt', inputText(tagVec(ii)).tag]);
    inputText(tagVec(ii)).promptString = get(promptTag, 'string');
    inputText(tagVec(ii)).enable = get(currentTag, 'enable');
    % get the current string
    if strcmpi(get(currentTag, 'style'), 'edit')
        inputText(tagVec(ii)).answerString = deblank(get(currentTag, 'string'));
    end
    % check if the value is empty or not
    if icatb_findstr(get(currentTag, 'style'), 'popup')
        if isempty(inputText(tagVec(ii)).value)
            inputText(tagVec(ii)).value = 1;
        end
    end
end

menuData.inputText = inputText;
set(menuH, 'userdata', menuData);

% get the information for the following tags to attach to sesInfo
tagsNeeded = {'numOfPC1', 'numOfPC2', 'numReductionSteps'};
% loop over required tags
for ii = 1:length(tagsNeeded)
    currentTag = tagsNeeded{ii};
    currentHandle = findobj(handles, 'tag', currentTag);
    enableVal = get(currentHandle, 'enable');
    % check if the control is edit or popup
    if strcmpi(get(currentHandle, 'style'), 'edit')
        currentVal = str2num(deblank(get(currentHandle, 'string')));
    else
        currentVal = str2num(deblank(get(currentHandle, 'string')));
        currentVal = currentVal(get(currentHandle, 'value'));
    end
    % end for checking
    % if data reduction steps control is enabled
    if strcmp(currentTag, 'numReductionSteps')
        if strcmpi(enableVal, 'on')
            sesInfo.userInput = setfield(sesInfo.userInput, currentTag, currentVal);
        end
    else
        sesInfo.userInput = setfield(sesInfo.userInput, currentTag, currentVal);
    end
end
% end for checking

% set sesInfo structure
handles_data.sesInfo = sesInfo;
set(mainHandle, 'userdata', handles_data);
delete(handles);


function closeCallback(handleObj, event_data, handles)
% closes the figure window

% Close the current figure
try
    delete(handles);
catch
    %rethrow(lasterror);
end

function dataRedCallback(hObject, event_data, inputH, handles, menuH)
% data reduction callback
% Enable PC edit boxes as needed

menuData = get(menuH, 'userdata');

% main figure data
handles_data = get(handles, 'userdata');
sesInfo = handles_data.sesInfo;

prefix1 = 'numOfPC1';
prefix2 = 'numOfPC2';

PC1 = findobj(inputH, 'tag', prefix1);
PC2 = findobj(inputH, 'tag', prefix2);

prompt1 = findobj(inputH, 'tag', ['prompt', prefix1]);
prompt2 = findobj(inputH, 'tag', ['prompt', prefix2]);

getStr = get(hObject, 'string');
getVal = get(hObject, 'value');
numRedSteps = str2num(deblank(getStr(getVal, :)));
if (numRedSteps == 2)
    
    % display the object with PC1 and set the string to PC/IC
    set(prompt1, 'string', 'Number Of PC1');
    set(PC1, 'enable', 'on');
    
    set(prompt2, 'string', 'Number Of PC/IC (Step 2)');
    set(PC2, 'enable', 'inactive', 'string', num2str(sesInfo.userInput.numComp));
    
else
    % display the object with PC1 and set the string to PC/IC
    set(prompt1, 'string', 'Number Of PC/IC (Step 1)');
    set(PC1, 'string', num2str(sesInfo.userInput.numComp), 'enable', 'inactive');
    
    set(prompt2, 'string', 'Number Of PC2');
    set(PC2, 'enable', 'inactive', 'string', '0');
end

function selectMaskCallback(hObject, event_data, handles, menuH)

icatb_defaults;

% select mask callback
handles_data = get(handles, 'userdata');
sesInfo = handles_data.sesInfo;
getString = get(hObject, 'string');
getMask = get(hObject, 'value');
maskTag = get(hObject, 'tag');

% Get modality type
[modalityType, dataTitle] = icatb_get_modality;

maskFilter = '*.img;*.nii';
maskFigTitle = 'Select a mask file in analyze or nifti format';
fileType = 'image';

try
    
    maskFile = [];
    
    % select mask in 3D Analyze data
    if getMask == 2
        [P] = icatb_selectEntry('filter', maskFilter, 'title', maskFigTitle, 'fileType', fileType, ...
            'fileNumbers', 1);
        if ~isempty(P)
            maskFile = P;
        else
            disp('No mask file selected. Setting to default mask ...');
            set(hObject, 'value', 1);
        end
    end
    
    
    getMask = get(hObject, 'value');
    
    if getMask == 2
        if ~isfield(sesInfo.userInput, 'files')
            error('Data must be selected first before selecting mask');
        end
        [dT, extns, maskDim] = icatb_get_countTimePoints(icatb_parseExtn(maskFile));
        [dT, extns, dims] = icatb_get_countTimePoints(icatb_parseExtn(deblank(sesInfo.userInput.files(1).name(1, :))));
        if length(find(maskDim == dims)) ~= length(maskDim)
            error('Error:MaskDim', 'Mask dimensions ([%s]) doesn''t match that of data dimensions ([%s])', ...
                num2str(maskDim), num2str(dims));
        end
    end
    
    % Detect change in mask
    if (isfield(sesInfo.userInput, 'maskFile'))
        oldVal = 1;
        if (~isempty(sesInfo.userInput.maskFile))
            oldVal = 2;
        end
        
        if (getMask ~= oldVal)
            setappdata(0, 'create_mask_gica', 1);
        else
            if (~isempty(sesInfo.userInput.maskFile))
                if (ispc && ~strcmpi(sesInfo.userInput.maskFile, maskFile))
                    setappdata(0, 'create_mask_gica', 1);
                elseif (~ispc && ~strcmp(sesInfo.userInput.maskFile, maskFile))
                    setappdata(0, 'create_mask_gica', 1);
                end
            end
        end
    end
    
    % create Mask
    sesInfo.userInput.maskFile = maskFile;
    handles_data.sesInfo = sesInfo;
    set(handles, 'userdata', handles_data);
    
    % get the menu user data
    menuData = get(menuH, 'userdata');
    % input text
    inputText = menuData.inputText;
    % mask index
    maskIndex = strmatch('maskFile', char(inputText.tag), 'exact');
    % input text structure
    inputText(maskIndex).value = getMask;
    % set the menu data
    menuData.inputText = inputText;
    set(menuH, 'userdata', menuData);
    
catch
    icatb_errorDialog(lasterr, 'Mask Error', 'modal');
end

function estimateCompCallback(handleObj, event_data, handles)
% Estimating components callback

% Do estimation
[ec, c, m, a, sesInfo] = icatb_estimateCompCallback(handleObj, handles);
handles_data = get(handles, 'userdata');
handles_data.sesInfo = sesInfo;
set(handles, 'userdata', handles_data);

function applyCallback(handleObj, event_data, handles, menuH)
% apply callback

icatb_defaults;
global WRITE_ANALYSIS_STEPS_IN_DIRS;
global CONSERVE_DISK_SPACE;

write_analysis_steps_in_dirs = WRITE_ANALYSIS_STEPS_IN_DIRS;
conserve_disk_space = CONSERVE_DISK_SPACE;

% get the userdata
%sesInfo = get(handles, 'userdata');
handles_data = get(handles, 'userdata');
set(handles, 'pointer', 'watch');

try
    
    sesInfo = handles_data.sesInfo;
    
    menuData = get(menuH, 'userdata');
    inputText = menuData.inputText;
    
    oldDir = sesInfo.userInput.pwd;
    
    % By default the data type is real
    if ~isfield(sesInfo.userInput, 'dataType')
        sesInfo.userInput.dataType = 'real';
    end
    
    % name the complex images
    [sesInfo, complexInfo] = icatb_name_complex_images(sesInfo, 'read');
    
    % See if the files exist or not
    if ~isfield(sesInfo.userInput, 'files')
        error('Functional data is not selected');
    end
    
    if ~isfield(sesInfo.userInput, 'designMatrix')
        sesInfo.userInput.designMatrix.name = [];
    end
    
    if ~isfield(sesInfo.userInput, 'maskFile')
        sesInfo.userInput.maskFile = [];
    end
    
    
    %% group ica type
    groupICATypeTag = 'group_ica_type';
    group_icaIndex = strmatch(groupICATypeTag, cellstr(char(inputText.tag)), 'exact');
    
    if (~isempty(group_icaIndex))
        groupICATypeOpts = lower(cellstr(inputText(group_icaIndex).answerString));
        groupICATypeVal = inputText(group_icaIndex).value;
        groupICAType = lower(groupICATypeOpts{groupICATypeVal});
        sesInfo.userInput = setfield(sesInfo.userInput, groupICATypeTag, groupICAType);
    end
    
    % ICA algorithm
    %[sesInfo.userInput.algorithm] = get(icaAlgoH, 'value');
    
    % Number of subjects and sessions
    numOfSub = sesInfo.userInput.numOfSub;
    numOfSess = sesInfo.userInput.numOfSess;
    
    icaAlgoH =  findobj(handles, 'tag', 'algorithm');
    icaVal = get(icaAlgoH, 'value');
    icaStr = get(icaAlgoH, 'string');
    
    if (strcmpi(deblank(icaStr(icaVal, :)), 'iva-gl') || strcmpi(deblank(icaStr(icaVal, :)), 'iva-l')  || strcmpi(deblank(icaStr(icaVal, :)), 'moo-icar') ...
            || strcmpi(deblank(icaStr(icaVal, :)), 'gig-ica') || strcmpi(deblank(icaStr(icaVal, :)), 'constrained ica (spatial)'))
        sesInfo.userInput.numReductionSteps = 1;
    else
        if (sesInfo.userInput.numReductionSteps > 2)
            sesInfo.userInput.numReductionSteps = 2;
        end
        %
        %         if (sesInfo.userInput.numReductionSteps == 1)
        %             sesInfo.userInput.numOfGroups1 = 1;
        %         end
        
    end
    
    
    % number of reduction steps
    numReductionSteps = sesInfo.userInput.numReductionSteps;
    
    
    % Initialise groups 2 and 3
    sesInfo.userInput.numOfGroups2 = 0;
    sesInfo.userInput.numOfGroups3 = 0;
    sesInfo.userInput.numOfPC3 = 0;
    
    % check if the components is an integer or not
    [currentVar] = check_var_integer(num2str(sesInfo.userInput.numComp));
    
    
    if isempty(currentVar)
        error(['Not a valid Number for IC.']);
    end
    
    % Pre-processing data
    if (~isfield(sesInfo.userInput, 'preproc_type'))
        try
            preProcTag = strmatch('preproc_type', char(inputText.tag), 'exact');
            options = cellstr(inputText(preProcTag).answerString);
            selectedVal = inputText(preProcVisuoTag).value;
            sesInfo.userInput.preproc_type = lower(options{selectedVal});
        catch
            sesInfo.userInput.preproc_type = 'remove mean per timepoint';
        end
    end
    
    
    gPCAVal = strmatch('group_pca_type', char(inputText.tag), 'exact');
    if (~isempty(gPCAVal))
        options = cellstr(inputText(gPCAVal).answerString);
        selectedVal = inputText(gPCAVal).value;
        sesInfo.userInput.group_pca_type = lower(options{selectedVal});
        
        if (strcmpi(sesInfo.userInput.group_pca_type, 'grand mean') && any(sesInfo.userInput.diffTimePoints ~= sesInfo.userInput.diffTimePoints(1)))
            if (strcmpi(icatb_get_modality, 'fmri'))
                error('Please select the same no. of timepoints between subjects if you want to select the grand mean group PCA');
            elseif(strcmpi(icatb_get_modality, 'eeg'))
                error('Please select the same no. of electrodes between subjects if you want to select the grand mean group PCA');
            end
        end
    else
        sesInfo.userInput.group_pca_type = 'subject specific';
    end
    
    % PCA Type
    if (~isfield(sesInfo.userInput, 'pcaType'))
        pcaTag = strmatch('pcaType', char(inputText.tag), 'exact');
        options = cellstr(inputText(pcaTag).answerString);
        selectedVal = inputText(pcaTag).value;
        sesInfo.userInput.pcaType = lower(options{selectedVal});
    end
    
    sesInfo.userInput = icatb_check_pca_opts(sesInfo.userInput);
    
    
    %% Back reconstruction type
    if (~isfield(sesInfo.userInput, 'backReconType'))
        try
            backReconTag = strmatch('backReconType', char(inputText.tag), 'exact');
            options = cellstr(inputText(backReconTag).answerString);
            selectedVal = inputText(backReconTag).value;
            sesInfo.userInput.backReconType = lower(options{selectedVal});
        catch
            sesInfo.userInput.backReconType = 'regular';
        end
    end
    
    % get the strings associated with the calibrate and Z-score
    adjustTag = strmatch('scaleType', char(inputText.tag), 'exact');
    options = cellstr(inputText(adjustTag).answerString);
    selStr = deblank(options{inputText(adjustTag).value});
    selectedVal = strmatch(lower(selStr), lower(icatb_scaleICA), 'exact');
    
    if (isempty(selectedVal))
        selectedVal = 1;
    end
    
    
    sesInfo.userInput.scaleType = selectedVal - 1; % scale type
    
    %--Acknowledge creaters of ICA algorithms
    %icatb_acknowledgeCreators;
    
    if sesInfo.userInput.numComp < 2
        error(['Select Number of IC to be more than or equal to 2.']);
    end
    
    if (sesInfo.userInput.numComp > min(sesInfo.userInput.diffTimePoints))
        error(['No. of components selected cannot be greater than ', num2str(min(sesInfo.userInput.diffTimePoints))]);
    end
    
    
    % check data reduction steps
    if numReductionSteps == 1
        % 1 reduction step
        sesInfo.userInput.numOfPC2 = 0;
        sesInfo.userInput.numOfPC1 = sesInfo.userInput.numComp;
    else
        % 2 reduction steps
        
        % Error check for PC1
        % check if the components is an integer or not
        [currentVar] = check_var_integer(num2str(sesInfo.userInput.numOfPC1));
        
        if isempty(currentVar)
            error(['Not a valid Number for PC (Step 1). Please click on SetupICA-Defaults ', ...
                'menu to see Number of PC (step 1).']);
        end
        
        if sesInfo.userInput.numOfPC1 < 2
            error(['Select Number for PC (Step 1) to be more than or equal to 2. Please click on SetupICA-Defaults menu', ...
                'to see Number of PC (step 1).']);
        end
        
        sesInfo.userInput.numOfPC2 = sesInfo.userInput.numComp;
        numOfGroups2 = 1;
        sesInfo.userInput.numOfGroups2 = numOfGroups2;
    end
    % end for checking data reduction steps
    
    drawnow;
    
    if (~isfield(sesInfo.userInput, 'default_mask_opts'))
        sesInfo.userInput.default_mask_opts = icatb_default_mask_opts;
    end
    
    icatb_chk_changes_in_mask(sesInfo);
    
    if (isappdata(0, 'create_mask_gica'))
        % Create mask
        sesInfo = icatb_update_mask(sesInfo);
        rmappdata(0, 'create_mask_gica');
    end
    
    drawnow;
    
    % Check the data reduction numbers
    parameter_err_chk(sesInfo);
    
    set(handles, 'pointer', 'arrow');
    
    %% Execute the ICA Algorithm callback
    sesInfo = icaAlgorithmCallback(sesInfo, handles);
    
    %write_analysis_steps_in_dirs = WRITE_ANALYSIS_STEPS_IN_DIRS;
    if (isempty(write_analysis_steps_in_dirs))
        write_analysis_steps_in_dirs = 0;
    end
        
    if (isempty(conserve_disk_space))
        conserve_disk_space = 0;
    end
    
    sesInfo.userInput.write_analysis_steps_in_dirs = write_analysis_steps_in_dirs;
    sesInfo.userInput.conserve_disk_space = conserve_disk_space;
    
    %--keeps track whether variables have been initialized
    sesInfo.isInitialized = 0;
    
    % check if subjects file exist or not
    subFile = fullfile(sesInfo.userInput.pwd, [sesInfo.userInput.prefix, 'Subject.mat']);
    
    % check if subject file exists or not
    if ~exist(subFile)
        clear files;  clear SPMFiles;
        files = sesInfo.userInput.files;
        SPMFiles.name = sesInfo.userInput.designMatrix.name;
        icatb_save(subFile, 'numOfSub', 'numOfSess', 'files', 'SPMFiles');
        clear files;
    end
    
    icatb_save(sesInfo.userInput.param_file, 'sesInfo');
    disp('');
    displayString = [' Parameters are saved in ', sesInfo.userInput.param_file];
    % display the parameters
    disp(displayString);
    fprintf('\n');
    delete(handles);
    
    disp('Please run the analysis using the same parameter file');
    fprintf('\n');
    
catch
    set(handles, 'pointer', 'arrow');
    icatb_errorDialog(lasterr, 'Setup ICA Error');
    % rethrow last error
    icatb_displayErrorMsg;
    
end

function setCompCallback(hObject, event_data, handles, autoFillH)
% set user data
% check if the parameter file exists or not:
handles_data = get(handles, 'userdata');
sesInfo = handles_data.sesInfo;

numComp = get(hObject, 'string');

outputDir = sesInfo.userInput.pwd;
if isfield(sesInfo.userInput, 'param_file')
    [p, param_file] = fileparts(sesInfo.userInput.param_file);
    param_file = fullfile(outputDir, [param_file, '.mat']);
    numRedSteps = sesInfo.userInput.numReductionSteps;
    
    sesInfo.userInput.numComp = str2num(numComp);
    
    answers = get(autoFillH, 'string');
    getVal = get(autoFillH, 'value');
    
    % selected str
    selectedStr = lower(deblank(answers(getVal, :)));
    
    PCBefore = sesInfo.userInput.numComp;
    
    [minTp, minTpInd] = min(sesInfo.userInput.diffTimePoints);
    PCBefore = round(min([minTp, 1.5*PCBefore]));
    
    %     if (min(sesInfo.userInput.diffTimePoints) > 30)
    %         PCBefore = 30;
    %     end
    
    % Autofill components
    if numRedSteps == 1
        sesInfo.userInput.numOfPC1 = sesInfo.userInput.numComp;
        sesInfo.userInput.numOfPC2 = 0;
    else
        if strcmpi(selectedStr, 'yes')
            sesInfo.userInput.numOfPC1 = PCBefore; %sesInfo.userInput.numComp;
        end
        sesInfo.userInput.numOfPC2 = sesInfo.userInput.numComp;
    end
    % end for doing autofill
    
    handles_data.sesInfo = sesInfo;
    set(handles, 'userdata', handles_data);
    
end


% object callback for the popup handle for autofilling components
function autoFillCallback(hObject, event_data, handles, compH)

% execute the component callback
setCompCallback(compH, [], handles, hObject);


function verSliderCallback(handleObj, event_data, handles)
% vertical slider callback

% execute the slider callback
icatb_verSliderCallback(handleObj, handles);


function sesInfo = icaAlgorithmCallback(sesInfo, handles)
% ICA algorithm Callback

icatb_defaults;
global ICAOPTIONS_WINDOW_DISPLAY;

menuH = findobj(handles, 'tag', 'SetupICA-Defaults');
set_data_para(handles, menuH);

hObject = findobj(handles, 'tag', 'algorithm');
sel_algorithm = get(hObject, 'value');
allAlgorithms = cellstr(get(hObject, 'string'));
mask_ind = sesInfo.userInput.mask_ind;
numComp = sesInfo.userInput.numComp;


%% group ica type
groupICATypeTag = 'group_ica_type';
groupICAAnalysisHandle = findobj(handles, 'tag', groupICATypeTag);
useTemporalICA = 0;
if (~isempty(groupICAAnalysisHandle))
    groupICATypeOpts = cellstr(get(groupICAAnalysisHandle, 'string'));
    groupICATypeVal = get(groupICAAnalysisHandle, 'value');
    useTemporalICA = strcmpi('temporal', lower(groupICATypeOpts{groupICATypeVal}));
end

% Size of the reduced data
if (~useTemporalICA)
    dataSize = [numComp, length(mask_ind)];
else
    dataSize = [numComp, sum(sesInfo.userInput.diffTimePoints)];
end


modalityType = icatb_get_modality;


ica_options_visibility = 'on';

% Update selected ICA algorithm
sesInfo.userInput.algorithm = sel_algorithm;


algorithmName = lower(allAlgorithms{sel_algorithm});

if strcmpi(algorithmName, 'moo-icar')
    algorithmName = 'gig-ica';
end

if (useTemporalICA)
    if (strcmpi(algorithmName, 'gig-ica') || strcmpi(algorithmName, 'constrained ica (spatial)')  || ...
            strcmpi(algorithmName, 'semi-blind infomax'))
        error(['Temporal ICA cannot be run using algorithm ', algorithmName]);
    end
end

if (strcmpi(algorithmName, 'semi-blind infomax'))
    % if the number of data-sets is greater than 1
    if sesInfo.userInput.numOfSub*sesInfo.userInput.numOfSess > 1
        error('Error:SBICA', ['Presently SBICA works with only one data set.', ...
            '\nPlease select another algorithm if you want to run multiple data-sets']);
    end
end


% Select the options of the ICA algorithm other than Erica ICA algorithm
matchedInd = strmatch(algorithmName, lower({'Infomax', 'Fast ICA', 'SDD ICA', 'Semi-blind Infomax', 'Constrained ICA (Spatial)', 'FBSS', 'ERBM', 'IVA-GL', 'IVA-L'}), 'exact');

sesInfo.userInput.ICA_Options = {};

if (~isempty(matchedInd))
    
    if strcmpi(ICAOPTIONS_WINDOW_DISPLAY, 'on')
        sesInfo.userInput.ICA_Options = icatb_icaOptions(dataSize, sesInfo.userInput.algorithm, ica_options_visibility);
    end
    % Check if cancel button for input dialog is selected
    if isempty(sesInfo.userInput.ICA_Options)
        disp('ICA Options are not Selected');
    end
    
    if (strcmpi(algorithmName, 'semi-blind infomax'))
        
        % design matrix name
        checkFile = sesInfo.userInput.designMatrix(1).name;
        
        % Check design matrix and regressor information
        % get the spm2 design matrix file
        if isempty(checkFile)
            checkFile = icatb_selectEntry('typeSelection', 'single', 'typeEntity', 'file', 'filter', ...
                '*.mat', 'title', 'Select SPM2/SPM5/SPM8 design matrix');
        end
        % edit box tag
        editboxTag = 'prefix';
        
        editPrefix = deblank(get(findobj(handles, 'tag', editboxTag), 'string'));
        
        % Pull out the time courses from the spm design matrix
        prompt = ['Select subject ', num2str(1), ' timecourse constraints'];
        % use the new function to load the spm design matrix
        [spmData] = icatb_loadSPM_new('spmName', checkFile, 'countTimePoints', ...
            sesInfo.userInput.diffTimePoints(1), 'typeStr', 'multiple', 'data_sessionNumber', 1, ...
            'check_spm_design_matrix', 'yes', 'check_regressors', 'yes', 'flag_selecting_regressors', ...
            'yes', 'get_timecourses', 'yes', 'maxSelection', 2, 'prompt', prompt);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % save the design matrix information in Subject file and
        % parameter file
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        SPMFiles.name = checkFile;
        icatb_save(fullfile(sesInfo.userInput.pwd, [editPrefix, 'Subject.mat']), 'SPMFiles', '-append');
        
        sesInfo.userInput.designMatrix.name = checkFile;
        timecourse = spmData.timecourse; % get the time course from the spm data
        clear spmData;
        % Number of ICA options
        numICAOptions = length(sesInfo.userInput.ICA_Options);
        % ICA Options
        ICAOptions = sesInfo.userInput.ICA_Options;
        % Append the ICA Options
        ICAOptions{numICAOptions + 1} = 'TC';
        ICAOptions{numICAOptions + 2} = timecourse;
        % Pass the updated ICA Options and append the structure
        % with all the time courses that are present
        sesInfo.userInput.ICA_Options = ICAOptions;
        clear modelT;
    end
    
end

if (strcmpi(algorithmName, 'moo-icar'))
    algorithmName = 'gig-ica';
end

%%%% Multi-fixed ICA algorithm for constraining spatial maps %%%
if (strcmpi(algorithmName, 'constrained ica (spatial)') || strcmpi(algorithmName, 'gig-ica'))
    
    if ~strcmpi(algorithmName, 'gig-ica')
        spatial_references = icatb_selectEntry('title', 'Select spatial reference files ...', ...
            'typeEntity', 'file', 'typeSelection', 'multiple', 'filter', '*.img;*.nii', 'fileType', 'image');
        if isempty(spatial_references)
            error('Spatial reference files need to be selected inorder to use constrained ICA');
        end
    else
        if (strcmpi(modalityType, 'fmri'))
            gigFileP = '*agg*comp*.img;*agg*comp*.nii';
            gigFileTitle = 'Select aggregate images from group ica ...';
        else
            gigFileTitle = 'Select reference images ...';
            gigFileP = '*.img;*.nii';
        end
        spatial_references = icatb_selectEntry('title', gigFileTitle, 'typeEntity', 'file', 'typeSelection', 'multiple', 'filter', gigFileP, 'fileType', 'image');
        if isempty(spatial_references)
            error('Aggregate images need to be selected in order to use GIG-ICA');
        end
    end
    
    
    % Get the count for spatial reference files
    numSpatialFiles = icatb_get_countTimePoints(spatial_references);
    
    fileNumbers = (1:numSpatialFiles);
    
    % get the spatial reference data
    [images, imHInfo] = icatb_loadData(spatial_references);
    
    imDims = imHInfo(1).DIM(1:3);
    funcDims = sesInfo.userInput.HInfo.DIM(1:3);
    
    if (length(find((imDims == funcDims) ~= 0)) ~= length(funcDims))
        error('Error:Dimensions', 'Spatial reference image dimensions ([%s]) are not equal to functional image dimensions ([%s])', ...
            num2str(imDims), num2str(funcDims));
    end
    
    images = reshape(images, prod(imDims), numSpatialFiles);
    images = (images(sesInfo.userInput.mask_ind, :))';
    
    ICAOptions = sesInfo.userInput.ICA_Options;
    
    % Update ICA Options
    sesInfo.userInput.ICA_Options = {'ref_data', images, ICAOptions{:}};
    sesInfo.userInput.numComp = length(fileNumbers);
    sesInfo.userInput.numOfPC1 = sesInfo.userInput.numComp;
    
end
% End for implementing the Multi-fixed ICA algorithm

function parameter_err_chk(sesInfo)
% check the parameters

tempSess.userInput = sesInfo.userInput;
tempSess.inputFiles = sesInfo.userInput.files;
tempSess.numReductionSteps = sesInfo.userInput.numReductionSteps;
tempSess.diffTimePoints = sesInfo.userInput.diffTimePoints;
if (~isfield(sesInfo.userInput, 'mask_ind'))
    sesInfo = icatb_update_mask(sesInfo);
end
mask_ind = sesInfo.userInput.mask_ind;
tempSess.mask_ind = mask_ind;

fprintf('\n');

[tempSess] = icatb_dataReductionSetup(tempSess, 'noprint');

% check to make sure valid parameters
icatb_parameterErrorCheck(tempSess, 'noprint');

clear tempSess;


function preprocCallback(hObject, event_data, handles)
%% Data pre-processing callback
%

handles_data = get(handles, 'userdata');

preProcOptions = cellstr(get(hObject, 'string'));

selectedOption = get(hObject, 'value');

%% Set figure data
handles_data.sesInfo.userInput.preproc_type = lower(preProcOptions{selectedOption});

if (strcmpi(handles_data.sesInfo.userInput.preproc_type, 'intensity normalization'))
    sH = findobj(gcbf, 'tag', 'scaleType');
    if (~isempty(sH))
        val = get(sH, 'value');
        if ((val == 2) || (val == 3))
            dialogH = icatb_dialogBox('title', 'Scaling Type', 'textbody', 'Don''t use component conversion to data units or z-scores when using intensity normalization.', ...
                'texttype', 'large');
            waitfor(dialogH);
        end
    end
end

set(handles, 'userdata', handles_data);

function pcaTypeCallback(hObject, event_data, handles)
%% Group PCA Callback

handles_data = get(handles, 'userdata');

pcaOptions = cellstr(get(hObject, 'string'));

selectedOption = get(hObject, 'value');

pcaType = lower(pcaOptions{selectedOption});

%% Set figure data
handles_data.sesInfo.userInput.pcaType = pcaType;

handles_data.sesInfo.userInput = icatb_check_pca_opts(handles_data.sesInfo.userInput);

handles_data.sesInfo.userInput.pca_opts = icatb_pca_options(pcaType, handles_data.sesInfo.userInput.pca_opts, 'on');

set(handles, 'userdata', handles_data);


function backReconCallback(hObject, event_data, handles)
%% Back reconstruction type callback
%

handles_data = get(handles, 'userdata');

backReconOptions = cellstr(get(hObject, 'string'));

selectedOption = get(hObject, 'value');

%% Set figure data
handles_data.sesInfo.userInput.backReconType = lower(backReconOptions{selectedOption});
set(handles, 'userdata', handles_data);


function WhichAnalysisCallback(hObject, event_data, handles)
%% Analysis Type callback
%

icatb_defaults;
global NUM_RUNS_GICA;


handles_data = get(handles, 'userdata');

selObject = get(hObject, 'value');

handles_data.sesInfo.userInput.which_analysis = selObject;

%% Open ICASSO options
if (selObject == 2)
    % ICASSO
    if (isfield(handles_data.sesInfo.userInput, 'icasso_opts'))
        icasso_opts = icatb_get_icasso_opts(handles_data.sesInfo.userInput.icasso_opts);
    else
        icasso_opts = icatb_get_icasso_opts;
    end
    
    handles_data.sesInfo.userInput.icasso_opts = icasso_opts;
    
elseif (selObject == 3)
    % MST
    
    % open input dialog box
    prompt = {'Enter no. of ICA Runs:'};
    dlg_title = 'ICA Runs';
    num_lines = 1;
    
    if (isfield(handles_data.sesInfo.userInput, 'mst_opts'))
        num_ica_runs = handles_data.sesInfo.userInput.mst_opts.num_ica_runs;
    else
        num_ica_runs = NUM_RUNS_GICA;
    end
    
    if (num_ica_runs <= 1)
        num_ica_runs = 2;
    end
    
    def = {num2str(num_ica_runs)};
    num_ica_runs = icatb_inputdlg2(prompt, dlg_title, num_lines, def);
    
    if (isempty(num_ica_runs))
        error('!!!No. of runs is not selected');
    end
    
    num_ica_runs = str2num(num_ica_runs{1});
    
    handles_data.sesInfo.userInput.mst_opts.num_ica_runs = num_ica_runs;
    
end

set(handles, 'userdata', handles_data);

function icaTypeCallback(hObject, event_data, handles)
%% ICA Type callback
%

menuH = findobj(handles, 'tag', 'SetupICA-Defaults');
numCompH = findobj(handles, 'tag', 'numComp');
set(numCompH, 'enable', 'on');

icaTypes = cellstr(get(hObject, 'string'));
icaVal = get(hObject, 'value');
whichAnalysisH = findobj(handles, 'tag', 'which_analysis');
set(whichAnalysisH, 'enable', 'on');
% if (strcmpi(icaTypes{icaVal}, 'iva-gl'))
%     set(whichAnalysisH, 'enable', 'inactive');
% end

if (strcmpi(icaTypes{icaVal}, 'moo-icar') || strcmpi(icaTypes{icaVal}, 'gig-ica') || strcmpi(icaTypes{icaVal}, 'constrained ica (spatial)'))
    set(numCompH, 'enable', 'inactive');
    set(whichAnalysisH, 'enable', 'inactive');
end

set_data_para(handles, menuH);


function group_ica_typeCallback(hObject, event_data, handles)
%% Group ica type callback
%

handles_data = get(handles, 'userdata');
menuH = findobj(handles, 'tag', 'SetupICA-Defaults');

opts = cellstr(get(hObject, 'string'));
val = get(hObject, 'value');
menuH = findobj(handles, 'tag', 'SetupICA-Defaults');
menuData = get(menuH, 'userdata');
inputText = menuData.inputText;
chk = strmatch('numReductionSteps', cellstr(char(menuData.inputText.tag)), 'exact');

sH = findobj(gcbf, 'tag', 'numReductionSteps');

if (strcmpi(opts{val}, 'temporal'))
    if (~isempty(chk))
        redVal = strmatch('1', cellstr(menuData.inputText(chk).answerString), 'exact');
        if (menuData.inputText(chk).value ~= redVal)
            dialogH = icatb_dialogBox('title', 'Temporal ICA', 'textbody', 'Number of reduction steps is set to 1 when using temporal ica.', ...
                'texttype', 'large');
            waitfor(dialogH);
        end
        menuData.inputText(chk).value = redVal;
        menuData.inputText(chk).enable = 'inactive';
        handles_data.sesInfo.userInput.numReductionSteps = 1;
        set(handles, 'userdata', handles_data);
        set(sH, 'value', redVal);
        set(sH, 'enable', 'off');
    end
else
    menuData.inputText(chk).enable = 'on';
    set(sH, 'enable', 'on');
end


set(menuH, 'userdata', menuData);
set_data_para(handles, menuH);


function parallelCallback(hObject, event_data, handles)
%% Parallel callback
%


handles_data = get(handles, 'userdata');
val = get(hObject, 'value');
opts = cellstr(get(hObject, 'string'));
selectedOption = lower(opts{val});

num_workers = 4;
try
    num_workers = handles_data.sesInfo.userInput.parallel_info.num_workers;
catch
end


if (strcmpi(selectedOption, 'parallel'))
    prompt = {'Enter number of sessions/workers:'};
    dlg_title = 'Number of sessions/workers';
    num_lines = 1;
    def = {num2str(num_workers)};
    % save the file with the file name specified
    answer = icatb_inputdlg2(prompt, dlg_title, num_lines, def);
    if (~isempty(answer) && ~isempty(answer{1}))
        num_workers = str2num(answer{1});
    else
        set(hObject, 'value', 1);
        selectedOption = lower(opts{1});
    end
end

drawnow;

parallel_info.mode = selectedOption;
parallel_info.num_workers = num_workers;
handles_data.sesInfo.userInput.parallel_info = parallel_info;

set(handles, 'userdata', handles_data);
