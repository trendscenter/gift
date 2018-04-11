function icatb_displaySesInfo

% displays analysis information
% get the input from the structure sesInfo
% plot figure
% three buttons
% define three callbacks

try
    % show directions
    %icatb_new_directions('analysis-info');
    icatb_defaults;
    
    % Screen Color Defaults
    global BG_COLOR;
    global BG2_COLOR;
    global BUTTON_COLOR;
    global FG_COLOR;
    global AXES_COLOR;
    global FONT_COLOR;
    global BUTTON_FONT_COLOR;
    
    % fonts
    global UI_FONTNAME;
    global UI_FONTUNITS;
    global UI_FS;
    global PARAMETER_INFO_MAT_FILE;
    
    % load parameters file
    filterP = ['*', PARAMETER_INFO_MAT_FILE, '*.mat'];
    [P] = icatb_selectEntry('typeEntity', 'file', 'title', 'Select Parameter File', 'filter', filterP);
    if isempty(P)
        error('Parameter file is not selected');
    end
    oldDir = fileparts(P);
    cd(oldDir);
    load(P);
    if ~exist('sesInfo', 'var')
        error('Please select the ICA parameter file');
    end
    
    if sesInfo.isInitialized == 0
        error('Please run the analysis to view the analysis information.');
    end
    
    [modalityType, dataTitle, compSetFields] = icatb_get_modality;
    
    if isfield(sesInfo, 'modality')
        if ~strcmpi(sesInfo.modality, modalityType)
            if strcmpi(sesInfo.modality, 'fmri')
                error('You have selected the fMRI parameter file. Use GIFT toolbox to analysis information.');
            elseif strcmpi(sesInfo.modality, 'smri')
                error('You have selected the sMRI parameter file. Use SBM toolbox to analysis information.');
            else
                error('You have selected the EEG parameter file. Use EEGIFT toolbox to analysis information.');
            end
        end
    else
        sesInfo.modality = 'fmri';
    end
    
    % Figure for displaying analysis information
    [InputHandle] = icatb_getGraphics('Analysis Information', 'normal', 'SesInfo Display');
    
    set(InputHandle, 'menubar', 'none');
    
    if strcmpi(modalityType, 'fmri')
        helpLabel = 'GIFT-Help';
    else
        helpLabel = 'SBM-Help';
    end
    
    menu1H = uimenu('parent', InputHandle, 'label', helpLabel);
    menu2H = uimenu(menu1H, 'label', 'Analysis Info', 'callback', 'icatb_openHTMLHelpFile(''icatb_analysis_info.htm'');');
    
    % figure data
    set(InputHandle, 'userdata', sesInfo);
    
    
    % button positions
    xOffset = 0.02;
    yOffset = 0.02;
    buttonOrigin = yOffset;
    buttonWidth = (1 - 4*xOffset) / 3;
    buttonHeight = 0.05;
    
    % Parameter Info
    parameterPos = [xOffset yOffset buttonWidth buttonHeight];
    paramButton = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', ...
        parameterPos, 'string', 'Parameter Info', 'tag', 'parameterButton', 'callback', ...
        {@paramInfoCallback, InputHandle});
    
    % Data Reduction Info
    dataReductionPos = [parameterPos(1) + parameterPos(3) + xOffset yOffset buttonWidth buttonHeight];
    dataReductionButton = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', ...
        'position', dataReductionPos, 'string', 'Data Reduction Info', 'tag', 'dataReductionButton', 'callback', ...
        {@dataReductionInfoCallback, InputHandle});
    
    % File Output Info
    fileOutputPos = [dataReductionPos(1) + dataReductionPos(3) + xOffset yOffset buttonWidth buttonHeight];
    fileOutputButton = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'pushbutton', 'position', ...
        fileOutputPos, 'string', 'File Output Info', 'tag', 'fileOutputButton', 'callback', ...
        {@fileOutputInfoCallback, InputHandle});
    
    % Listbox Origin
    listOrigin = buttonHeight + 2*yOffset;
    
    % listbox width
    listWidth = (1 - 2*xOffset);
    
    % listbox height
    listHeight = (1 - 2*yOffset - listOrigin);
    
    % draw listbox
    % Listbox automatically takes care of scroll bar
    handle_scroll = icatb_uicontrol('parent', InputHandle, 'units', 'normalized', 'style', 'listbox', 'position', ...
        [xOffset listOrigin listWidth listHeight], 'string', [], 'horizontalalignment', 'center', 'tag', ...
        'analysisInfoListbox');
    
    % display parameter information
    paramInfoCallback(paramButton, [], InputHandle);
    
catch
    if exist('InputHandle', 'var')
        if ishandle(InputHandle)
            delete(InputHandle);
        end
    end
    icatb_displayErrorMsg;
end

%%%%%%%%%%%%%%%%%%%%%%% Function callbacks %%%%%%%%%%%%%%%%%%

% parameter information callback

function paramInfoCallback(hObject, evd, handles)

set(handles, 'pointer', 'watch');

% purpose: displays parameter information in the listbox

sesInfo = get(handles, 'userdata');

% Initialise listbox
handle_scroll = findobj(handles, 'tag', 'analysisInfoListbox');

set(handle_scroll, 'string', []);

% available ICA algorithms
icaAlgo = icatb_icaAlgorithm;
% selected ICA algorithm
selected_ica_algorithm = deblank(icaAlgo(sesInfo.algorithm, :));

D(1).string = 'Group ICA Parameter Infomation';

D(size(D,2)+1).string = '';

% Status
D(size(D,2)+1).string = 'Got Parameters From User: YES';
D(size(D,2)+1).string = 'Parameters Initialized: YES';

% Parameter Info
modalityType = 'fmri';
try
    modalityType = sesInfo.modality;
catch
end

if (strcmpi(modalityType, 'fmri'))
    D(size(D,2)+1).string = ['Number of Subjects : ', num2str(sesInfo.numOfSub)];
    D(size(D,2)+1).string = ['Number of Sessions : ', num2str(sesInfo.numOfSess)];
end

D(size(D,2)+1).string = ['Number of Independent Components : ', num2str(sesInfo.numComp)];
D(size(D,2)+1).string = ['ICA Algorithm : ', selected_ica_algorithm];

if (~strcmpi(modalityType, 'smri'))
    D(size(D,2)+1).string = ['Number Of Scans/Timepoints : ', num2str(sesInfo.numOfScans)];
else
    D(size(D,2)+1).string = ['Number Of Subjects : ', num2str(sesInfo.numOfScans)];
end

[modality, dataTitle] = icatb_get_modality;

if(isempty(sesInfo.userInput.maskFile))
    D(size(D,2)+1).string = ['Mask File : Default Mask Created From ', dataTitle, ' Data'];
else
    [pathstr,name] = fileparts(sesInfo.userInput.maskFile);
    D(size(D,2)+1).string = ['Mask File : ',name];
end

% Data pre-processing Type
preproc_options = icatb_preproc_data;
preproc_options = cellstr(preproc_options);
preproc_type = 'remove mean per timepoint';
if (isfield(sesInfo, 'preproc_type'))
    preproc_type = sesInfo.preproc_type;
end

preprocInd = strmatch(lower(preproc_type), lower(preproc_options), 'exact');

D(size(D,2)+1).string = ['Data Pre-processing Type : ', preproc_options{preprocInd}];

% PCA Type
pcaTypes = icatb_pca_options;
pcaTypes = cellstr(pcaTypes);
pcaType = 'standard';
if (isfield(sesInfo, 'backReconType'))
    pcaType = sesInfo.pcaType;
end

if (isnumeric(pcaType))
    pcaType = pcaTypes{pcaType};
end

pcaInd = strmatch(lower(pcaType), lower(pcaTypes), 'exact');

D(size(D,2)+1).string = ['PCA Type : ', pcaTypes{pcaInd}];

% Group PCA Type
if (~strcmpi(modalityType, 'smri'))
    
    multiSubGroupPCA = ((sesInfo.numReductionSteps == 1) && (sesInfo.numOfSub*sesInfo.numOfSess > 1));
    
    if (~multiSubGroupPCA || strcmpi(selected_ica_algorithm, 'iva-gl'))
        groupPCAOpts = char('Subject Specific', 'Grand Mean');
        group_pca_type = 'Subject Specific';
        if (isfield(sesInfo, 'group_pca_type'))
            group_pca_type = sesInfo.group_pca_type;
        end
        groupPCAInd = strmatch(lower(group_pca_type), lower(groupPCAOpts), 'exact');
        group_pca_type = deblank(groupPCAOpts(groupPCAInd, :));
        
        D(size(D,2)+1).string = ['Group PCA Type : ', group_pca_type];
    end
end

useTemporalICA = 0;
if (strcmpi(modalityType, 'fmri'))
    group_ica_type = 'spatial';
    try
        group_ica_type = sesInfo.group_ica_type;
    catch
    end
    useTemporalICA = strcmpi(group_ica_type, 'temporal');
    D(size(D,2)+1).string = ['Group ICA Type : ', upper(group_ica_type(1)), group_ica_type(2:end)];
end

if strcmpi(selected_ica_algorithm, 'moo-icar')
    selected_ica_algorithm = 'gig-ica';
end


if (~(strcmpi(selected_ica_algorithm, 'iva-gl') || strcmpi(selected_ica_algorithm, 'gig-ica') || strcmpi(selected_ica_algorithm, 'constrained ica (spatial)') || useTemporalICA))
    % Back Reconstruction Type
    backReconType = 'Regular';
    if (isfield(sesInfo, 'backReconType'))
        backReconType = sesInfo.backReconType;
    end
    
    backReconOptions = cellstr(icatb_backReconOptions);
    backReconInd = strmatch(lower(backReconType), lower(backReconOptions), 'exact');
    
    backReconType = backReconOptions{backReconInd};
    if ((sesInfo.numReductionSteps == 1) && (sesInfo.numOfSub*sesInfo.numOfSess > 1))
        backReconType = 'Spatial-temporal Regression';
    end
    
    D(size(D, 2)+1).string = ['Back Reconstruction Type : ', backReconType];
end

calibrateOptions = icatb_scaleICA;
D(size(D,2)+1).string = ['Scaling Components : ', deblank(calibrateOptions(sesInfo.scaleType + 1, :))];

stability_analysis = 'none';
whichAnalysis = 1;
try
    whichAnalysis = sesInfo.which_analysis;
catch
end

if (whichAnalysis == 2)
    stability_analysis = 'ICASSO';
end

if (whichAnalysis == 3)
    stability_analysis = 'MST';
end

D(size(D,2)+1).string = ['Stability analysis type : ', stability_analysis];

checkAnalysisMode = 'Serial';
try
    checkAnalysisMode = sesInfo.parallel_info.mode;
catch
end

D(size(D,2)+1).string = ['Group analysis mode: ', upper(checkAnalysisMode(1)), checkAnalysisMode(2:end)];

% Adding a scroll bar to the parameterInfo
Initial_number = 1;

% Initial text in parameter info
newText = cellstr(str2mat(D.string));

newText = textwrap(handle_scroll, newText);
set(handle_scroll, 'string', newText);

% Apply conditions for dialog box differently for different platforms
if ispc
    set(handle_scroll, 'enable', 'inactive');
else
    set(handle_scroll, 'enable', 'on');
end


set(handle_scroll, 'min', 0, 'max', 2);

% make no selection
set(handle_scroll, 'value', []);

set(handles, 'pointer', 'arrow');

% data reduction information callback
function dataReductionInfoCallback(hObject, evd, handles)

% purpose: displays data reduction information in the listbox

set(handles, 'pointer', 'watch');

sesInfo = get(handles, 'userdata');

% Initialise listbox
handle_scroll = findobj(handles, 'tag', 'analysisInfoListbox');

set(handle_scroll, 'string', []);

%Data Reduction Info
for i = 1: sesInfo.numReductionSteps
    D(1).string= 'Reduction Step Info';
    D(size(D,2)+1).string = '';
    D(size(D,2)+1).string = ['Reduction Step ',num2str(i) ,' Info:'];
    D(size(D,2)+1).string = ['Number of Groups Before Concatenation : ',num2str(sesInfo.reduction(i).numOfGroupsBeforeCAT)];
    D(size(D,2)+1).string = ['Number of Groups After Concatenation : ',num2str(sesInfo.reduction(i).numOfGroupsAfterCAT)];
    D(size(D,2)+1).string = ['Number of Previous Groups in New Groups : ',num2str(sesInfo.reduction(i).numOfPrevGroupsInEachNewGroupAfterCAT)];
    D(size(D,2)+1).string = ['Number of PC Before Data Reduction : ',num2str(sesInfo.reduction(i).numOfPCInEachGroupAfterCAT)];
    D(size(D,2)+1).string = ['Number of PC After Data Reduction : ',num2str(sesInfo.reduction(i).numOfPCAfterReduction)];
end


% Initial text in parameter info
newText = cellstr(str2mat(D.string));

newText = textwrap(handle_scroll, newText);
set(handle_scroll, 'string', newText);

% Apply conditions for dialog box differently for different platforms
if ispc
    set(handle_scroll, 'enable', 'inactive');
else
    set(handle_scroll, 'enable', 'on');
end

set(handle_scroll, 'min', 0, 'max', 2);

% make no selection
set(handle_scroll, 'value', []);

set(handles, 'pointer', 'arrow');

% file output information callback
function fileOutputInfoCallback(hObject, evd, handles)

% purpose: displays file output information information in the listbox

set(handles, 'pointer', 'watch');

icatb_defaults;

% Output File Indices
global GROUP_ICA_INDEX;
global SUBJECT_ICA_INDEX;
global MEAN_ALL_INDEX;
global MEAN_INDEX;
global TMAP_INDEX;
global STD_INDEX;

sesInfo = get(handles, 'userdata');

% Initialise listbox
handle_scroll = findobj(handles, 'tag', 'analysisInfoListbox');

set(handle_scroll, 'string', []);

conserve_disk_space = 0;
if (isfield(sesInfo, 'conserve_disk_space'))
    conserve_disk_space = sesInfo.conserve_disk_space;
end

%File Ouput Info
D(1).string = 'Output File Info';
D(size(D,2)+1).string = '';
D(size(D,2)+1).string = 'Parameter Infomation Stored as Matlab File in:';
D(size(D,2)+1).string = ['     ', sesInfo.param_file];

if (conserve_disk_space ~= 2)
    D(size(D,2)+1).string = 'Data Reduction Results Stored as Matlab File in: ';
    D(size(D,2)+1).string = ['     ', sesInfo.data_reduction_mat_file];
end

D(size(D,2)+1).string = 'ICA Results Stored as Matlab File in:';
D(size(D,2)+1).string = ['     ', sesInfo.ica_mat_file];

if (~conserve_disk_space)
    D(size(D,2)+1).string = 'Back Reconstructed Results Stored as Matlab File in:';
    D(size(D,2)+1).string = ['     ', sesInfo.back_reconstruction_mat_file];
    D(size(D,2)+1).string = 'Calibrated Components Stored as Matlab File in:';
    D(size(D,2)+1).string=['     ', sesInfo.calibrate_components_mat_file];
end

% check the number of output files
% check the naming of the component files
% get the first component name.

% number of output files
numOutputFiles = length(sesInfo.icaOutputFiles);

% loop over number of output files
for ii = 1:numOutputFiles
    numFiles = length(sesInfo.icaOutputFiles(ii).ses);
    for jj = 1:numFiles
        componentName = sesInfo.icaOutputFiles(ii).ses(jj).name(1, :);
        
        %         for kk = 1:size(sesInfo.icaOutputFiles(ii).ses(jj).name, 1)
        %             % check existence of the file
        %             %                 if ~exist(deblank(sesInfo.icaOutputFiles(ii).ses(jj).name(kk, :)), 'file')
        %             %                     error(['File ', deblank(sesInfo.icaOutputFiles(ii).ses(jj).name(kk, :)), ' doesn''t exist.']);
        %             %                 end
        %         end
        if (ispc)
            componentName = lower(componentName);
        end
        
        if icatb_findstr(componentName, 'mean') & icatb_findstr(componentName, 'all')
            D(size(D,2)+1).string = 'Mean Component for all Subjects and Sessions Stored as:';
            D(size(D,2)+1).string = ['     ', componentName];
        elseif icatb_findstr(componentName, 'mean')
            D(size(D,2)+1).string = ['Mean Component Results for Session ', num2str(jj), ' Stored as:'];
            D(size(D,2)+1).string = ['     ', componentName];
        elseif icatb_findstr(componentName, 'tmap')
            D(size(D,2)+1).string = ['TMap Component Results for Session ', num2str(jj), ' Stored as:'];
            D(size(D,2)+1).string =['     ', componentName];
        elseif icatb_findstr(componentName, 'std')
            D(size(D,2)+1).string = ['Std Component Results for Session ', num2str(jj), ' Stored as:'];
            D(size(D,2)+1).string =['     ', componentName];
        elseif icatb_findstr(componentName, 'sub')
            subNumber = str2num(regexprep(componentName, '.*sub(\d+).*', '$1'));
            D(size(D,2)+1).string = ['Subject ', num2str(subNumber),  ' Component Results for Session ', num2str(jj), ...
                ' Stored as:'];
            D(size(D,2)+1).string = ['     ', componentName];
        end
    end
end

% Initial text in parameter info
newText = cellstr(str2mat(D.string));

newText = textwrap(handle_scroll, newText);
set(handle_scroll, 'string', newText);

% Apply conditions for dialog box differently for different platforms
if ispc
    set(handle_scroll, 'enable', 'inactive');
else
    set(handle_scroll, 'enable', 'on');
end

% set min, max, values for the listbox
set(handle_scroll, 'min', 0, 'max', 2);

% make no selection
set(handle_scroll, 'value', []);

set(handles, 'pointer', 'arrow');