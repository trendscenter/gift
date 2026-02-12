function icatb_setup_conn_ica(param_file)
%% Setup Connectivity ICA
%
% Inputs are:
% inputFile - Input file for batch analysis.
%
% Outputs:
% ICA Parameter file.

if ~exist('param_file', 'var')
    param_file = icatb_selectEntry('typeEntity', 'file', 'title', 'Select Parameter File', 'filter', '*ica*param*.mat');
    if isempty(param_file)
        error('Parameter file is not selected for analysis');
    end
end

load(param_file);

[outputDir, fN, extn] = fileparts(param_file);
if (isempty(outputDir))
    outputDir = pwd;
end

cd(outputDir);

param_file = fullfile(outputDir, [fN, extn]);

sesInfo.userInput.pwd = outputDir;
sesInfo.outputDir = outputDir;

sesInfo.userInput.param_file = param_file;

% Define input parameters
inputText = icatb_define_parameters(sesInfo);

[modalityType] = icatb_get_modality;

% Tags to be plotted in main figure window
tagsTobePlotted = {'numComp', 'numOfPC1', 'algorithm', 'which_analysis', 'backReconType', 'scaleType', 'parallel_info'};


handle_visibility = 'on';

drawnow;


%%%%%%%%% Draw figure and set up the controls %%%%%%%%

figureTag = 'setup_ICAGUI';

figHandle = findobj('tag', figureTag);

if ~isempty(figHandle)
    delete(figHandle);
end

figTitle = 'Setup Connectivity Domain ICA';

% Setup figure for GUI
[InputHandle] = icatb_getGraphics(figTitle, 'normal', figureTag, handle_visibility);

% store the directory information to sesInfo
sesInfo.userInput.modality = modalityType;

handles_data.sesInfo = sesInfo;

% set figure data
set(InputHandle, 'userdata', handles_data, 'Menubar', 'none');



algoIndex = strmatch('algoType', char(inputText.tag), 'exact');

inputText(algoIndex).tag = 'algorithm';


[dd, ia, ib] = intersect(cellstr(char(inputText.tag)), tagsTobePlotted);
inputText = inputText(ia);


% countTags_menu = 0;
%
% % total tags in menu
% tagsVec = zeros(1, length(inputText) - length(tagsTobePlotted));
%
% % store the tags that need to be plotted in a menu
% for ii = 1:length(inputText)
%     strIndex = strmatch(inputText(ii).tag, tagsTobePlotted, 'exact');
%     if isempty(strIndex)
%         countTags_menu = countTags_menu + 1;
%         tagsInMenu{countTags_menu} = inputText(ii).tag;
%         tagsVec(countTags_menu) = ii;
%     end
% end

% plot the input parameters to the parameter figure
icatb_plotInputPara('input_prefix', inputText, 'controls_to_plot', tagsTobePlotted, 'handles', ...
    InputHandle);

%%%%%%%%% End for drawing the figure and setting up the controls %%%%

for ii = 1:length(tagsTobePlotted)
    checkH = findobj(InputHandle, 'tag', tagsTobePlotted{ii});
    set(checkH, 'enable', 'on');
end


%% Set callback for cancel button
set(findobj(InputHandle, 'tag', 'Cancel'), 'callback', {@closeCallback, InputHandle});


%% Done callback
okHandle = findobj(InputHandle, 'tag', 'Done');
set(okHandle, 'callback', {@applyCallback, InputHandle});

%% Which Analysis Callback
WhichAnalysisH = findobj(InputHandle, 'tag', 'which_analysis');
set(WhichAnalysisH, 'callback', {@WhichAnalysisCallback, InputHandle});

%% ICA Algorithm callback
set(findobj(InputHandle, 'tag', 'algorithm'), 'callback', {@icaTypeCallback, InputHandle});

set(findobj(InputHandle, 'tag', 'parallel_info'), 'callback', {@parallelCallback, InputHandle});


%%%%%%%%%% End for setting up the function callbacks %%%%%%%%%%%%%%%%%

show_params(InputHandle);





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
    
    
elseif (selObject == 4)
    % cross isi
    
    % open input dialog box
    prompt = {'Enter no. of ICA Runs:'};
    dlg_title = 'ICA Runs';
    num_lines = 1;
    
    if (isfield(handles_data.sesInfo.userInput, 'cross_isi_opts'))
        num_ica_runs = handles_data.sesInfo.userInput.cross_isi_opts.num_ica_runs;
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
    
    handles_data.sesInfo.userInput.cross_isi_opts.num_ica_runs = num_ica_runs;
    
    
end

set(handles, 'userdata', handles_data);


function show_params(handles)

handles_data = get(handles, 'userdata');
sesInfo = handles_data.sesInfo;

%% Show selected params

numCompH = findobj(handles, 'tag', 'numComp');
numVoxels = length(sesInfo.userInput.mask_ind);
numComp = min([20, numVoxels]);
try
    numComp = sesInfo.userInput.numComp;
catch
end
set(numCompH, 'string', num2str(numComp));

set(numCompH, 'enable', 'on');

algoH = findobj(handles, 'tag', 'algorithm');
algorithm = 1;

try
    algorithm = sesInfo.userInput.algorithm;
catch
end

set(algoH, 'value', algorithm);

whichAnalysisH = findobj(handles, 'tag', 'which_analysis');
which_analysis = 1;
try
    which_analysis = sesInfo.userInput.which_analysis;
catch
end

set(whichAnalysisH, 'value', which_analysis);


backReconH = findobj(handles, 'tag', 'backReconType');
backReconType = 'gica';
backReconOptions = cellstr(get(backReconH, 'string'));
try
    backReconType = sesInfo.userInput.backReconType;
catch
end

if ~isnumeric(backReconType)
    backReconVal = strmatch(lower(backReconType), lower(backReconOptions), 'exact');
else
    backReconVal = backReconType;
end

if (isempty(backReconVal))
    backReconVal = 1;
end

set(backReconH, 'value', backReconVal);

try
    if isfield(sesInfo.userInput, 'scaleType')
        scaleOptions = {'No Scaling', 'Z-scores'};
        try
            matchIndex = strmatch(lower(sesInfo.userInput.scaleType), lower(cellstr(scaleOptions)), 'exact');
            if (isempty(matchIndex))
                matchIndex = 1;
            end
        catch
            matchIndex = 1;
        end
        scaleH = findobj(handles, 'tag', 'scaleType');
        set(scaleH, 'value', matchIndex);
    end
catch
end

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


function closeCallback(handleObj, event_data, handles)
% closes the figure window

% Close the current figure
try
    delete(handles);
catch
    %rethrow(lasterror);
end

function applyCallback(hObject, event_data, handles)
%% Apply Callback

icatb_defaults;
global ICAOPTIONS_WINDOW_DISPLAY;

handles_data = get(handles, 'userdata');

set(handles, 'pointer', 'watch');

drawnow;

sesInfo = handles_data.sesInfo;


%%ce021126 pc_multiplierH = findobj(handles, 'tag', 'pc_multiplier');
% pc_multiplier = str2num(get(pc_multiplierH, 'string'));
% sesInfo.userInput.pc_multiplier = pc_multiplier;

h_numOfPC1 = findobj(handles, 'tag', 'numOfPC1');
numOfPC1 = str2num(get(h_numOfPC1, 'string'));

numCompH = findobj(handles, 'tag', 'numComp');
numComp = str2num(get(numCompH, 'string'));
sesInfo.userInput.numComp = numComp;

if numComp > numOfPC1
    error('numComp must be same or less than numOfPC1')
end

numVoxels = length(sesInfo.userInput.mask_ind);
numofDatasets = sesInfo.userInput.numOfSub*sesInfo.userInput.numOfSess;
numReductionSteps = 2;
numOfPC3 = 0;

%ce021126 if (pc_multiplier < 1)
%     pc_multiplier = 1;
% end

if (numofDatasets == 1)
    numReductionSteps = 1;
end

numOfPC2 = numComp;

% if (numReductionSteps > 1)
%     numOfPC1 = min([numVoxels, ceil(pc_multiplier*numComp)]);
%     numOfPC2 = numComp;
% end

sesInfo.userInput.numReductionSteps = numReductionSteps;

sesInfo.userInput.numOfPC1 = numOfPC1;
sesInfo.userInput.numOfPC2 = numOfPC2;
sesInfo.userInput.numOfPC3 = numOfPC3;

sesInfo.userInput.numOfGroups1 = numofDatasets;
sesInfo.userInput.numOfGroups2 = 1;
sesInfo.userInput.numOfGroups3 = 0;

diffTimePoints =  repmat(length(sesInfo.userInput.mask_ind), 1, numofDatasets);
sesInfo.userInput.diffTimePoints = diffTimePoints;

algoH = findobj(handles, 'tag', 'algorithm');
sesInfo.userInput.algorithm = get(algoH, 'value');

whichAnalysisH = findobj(handles, 'tag', 'which_analysis');
which_analysis = get(whichAnalysisH, 'value');
sesInfo.userInput.which_analysis = which_analysis;

backReconH = findobj(handles, 'tag', 'backReconType');
backReconVal = get(backReconH, 'value');
backReconOptions = cellstr(get(backReconH, 'string'));
backReconType = lower(backReconOptions{backReconVal});
sesInfo.userInput.backReconType = backReconType;

scaleTypeH = findobj(handles, 'tag', 'scaleType');
scaleOptions = {'No Scaling', 'Z-scores'};
scaleVal = get(scaleTypeH, 'value');
scaleType = scaleOptions{scaleVal};
sesInfo.userInput.scaleType = scaleType;

dataSize = [numComp, numVoxels];

sesInfo.userInput.ICA_Options = icatb_icaOptions(dataSize, sesInfo.userInput.algorithm, ICAOPTIONS_WINDOW_DISPLAY);

sesInfo.isInitialized = 0;

param_file = sesInfo.userInput.param_file;
disp(['Saving session information in file ', param_file, ' ...']);
save(param_file, 'sesInfo');
disp('Done');
fprintf('\n');

set(handles, 'pointer', 'arrow');

try
    delete(handles);
catch
end

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

