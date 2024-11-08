function icatb_setup_fnc_ica(inputFile)
%% Setup FNC ICA
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
tagsTobePlotted = {'prefix', 'files', 'estimate_components', 'numComp', 'algorithm', 'which_analysis', 'scaleType'};


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

figTitle = 'Setup FNC ICA';

% Setup figure for GUI
[InputHandle] = icatb_getGraphics(figTitle, 'normal', figureTag, handle_visibility);

% store the directory information to sesInfo
sesInfo.userInput.pwd = oldDir;
sesInfo.userInput.outputDir = oldDir;
sesInfo.userInput.modality = modalityType;

handles_data.sesInfo = sesInfo;

% set figure data
set(InputHandle, 'userdata', handles_data, 'Menubar', 'none');


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
%setupDefaultsData.tagsTobePlotted = tagsInMenu;  setupDefaultsData.inputText = inputText;
%setupDefaultsData.tagsVec = tagsVec;

%% Open figure window when SetupICA-Defaults menu is clicked
%set(setupDefaults, 'callback', {@setupDefaultsCallback, InputHandle}, 'userdata', setupDefaultsData);

%% Answer function callbacks
set(findobj(InputHandle, 'tag', 'prefix'), 'callback', {@editboxCallback, InputHandle});

%% Load the functional files (data callback)
set(findobj(InputHandle, 'tag', 'files'), 'callback', ...
    {@dataCallback, InputHandle});

%% Set callback for cancel button
set(findobj(InputHandle, 'tag', 'Cancel'), 'callback', {@closeCallback, InputHandle});

estimateCompHandle = findobj(InputHandle, 'tag', 'estimate_components');

if ~isempty(estimateCompHandle)
    % set callback for estimate components callback
    set(estimateCompHandle, 'callback', {@estimateCompCallback, InputHandle});
end


%% Done callback
okHandle = findobj(InputHandle, 'tag', 'Done');
set(okHandle, 'callback', {@applyCallback, InputHandle});

%% Which Analysis Callback
WhichAnalysisH = findobj(InputHandle, 'tag', 'which_analysis');
set(WhichAnalysisH, 'callback', {@WhichAnalysisCallback, InputHandle});

%% ICA Algorithm callback
set(findobj(InputHandle, 'tag', 'algorithm'), 'callback', {@icaTypeCallback, InputHandle});


%%%%%%%%%% End for setting up the function callbacks %%%%%%%%%%%%%%%%%

if exist('output_prefix', 'var')
    prefixH = findobj(InputHandle, 'tag', 'prefix');
    set(prefixH, 'string', output_prefix);
end


% Execute this function for GUI
%editboxCallback(findobj(InputHandle, 'tag', 'prefix'), [], InputHandle);





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


function  editboxCallback (hObject, event_data, handles)
%% Editbox callback
%
handles_data = get(handles, 'userdata');

prefix = get(hObject, 'string');
if (isempty(prefix))
    error('Enter output prefix');
end

handles_data.sesInfo.userInput.prefix = prefix;

sesInfo = handles_data.sesInfo;

set(handles, 'pointer', 'watch');
param_file = fullfile(sesInfo.userInput.pwd, [sesInfo.userInput.prefix, '_ica_parameter_info.mat']);
if exist(param_file, 'file')
    load(param_file);
    handles_data.sesInfo = sesInfo;
    set(handles, 'userdata', handles_data);
    if (~isempty(sesInfo.userInput.dataInfo))
        
        show_params(handles);
        
        
    end
end

set(handles, 'userdata', handles_data);

set(handles, 'pointer', 'arrow');


function dataCallback(hObject, event_data, handles)
%% Data callback

handles_data = get(handles, 'userdata');

try
    handles_data.sesInfo.userInput.prefix;
catch
    error('Output prefix is not entered');
end

sesInfo = handles_data.sesInfo;

set(handles, 'pointer', 'watch');

sesInfo = icatb_select_fnc_data(sesInfo);

handles_data.sesInfo = sesInfo;

set(handles, 'userdata', handles_data);

show_params(handles);

set(handles, 'pointer', 'arrow');

set(hObject, 'foregroundcolor', [0, 1, 0]);


function estimateCompCallback(hObject, event_data, handles)
%% Estimate components

set(handles, 'pointer', 'watch');
drawnow;
handles_data = get(handles, 'userdata');
sesInfo = handles_data.sesInfo;

numSessions = length(sesInfo.userInput.dataInfo);
dataInfo = sesInfo.userInput.dataInfo;
for nSess = 1:numSessions
    sessName = dataInfo(nSess).name;
    files =  cellstr(dataInfo(nSess).files);
    if (nSess == 1)
        input_data_file_patterns = cell(size(files, 1), size(files, 2));
    end
    input_data_file_patterns(:, nSess) = cellstr(files);
end

contrast_vector = '';
try
    contrast_vector = sesInfo.userInput.contrast_vector;
catch
end

sesInfo.userInput.contrast_vector = contrast_vector;
inputData.contrast_vector = contrast_vector;

fnc_variable_mat_file = 'fnc_corrs_all';
try
    fnc_variable_mat_file = sesInfo.userInput.fnc_variable_mat_file;
catch
end

%% Load data
inputData.contrast_vector = contrast_vector;
inputData.fnc_variable_mat_file = fnc_variable_mat_file;
inputData.input_data_file_patterns = input_data_file_patterns;

fnc_file_name = fullfile(sesInfo.userInput.pwd, [sesInfo.userInput.prefix, '_fnc_data.mat']);

if exist(fnc_file_name, 'file')
    load(fnc_file_name, 'fnc_matrix');
else
    fnc_matrix = icatb_fnc_input(inputData);
    save(fnc_file_name, 'fnc_matrix');
end

fnc_matrix = icatb_mat2vec(fnc_matrix);

sesInfo.userInput.doEstimation = 1;


[numOfPC1, mdlVec] = order_selection(fnc_matrix);
sesInfo.userInput.est_opts.mdl = mdlVec;
sesInfo.userInput.est_opts.numOfPC1 = numOfPC1;
sesInfo.userInput.numComp = numOfPC1;


plotX.y = mdlVec;
plotX.x = (1:length(mdlVec));
plotX.title = 'Plot of MDL where minimum is the estimated dimensionality.';
msgString = ['The estimated independent components is found to be ', num2str(numOfPC1), ' using the MDL criteria.'];

helpButton = icatb_dialogBox('title', 'Estimated Components', 'textBody', msgString, 'textType', 'large', 'plotbutton', plotX);

handles_data.sesInfo = sesInfo;
set(handles, 'userdata', handles_data);

numCompH = findobj(handles, 'tag', 'numComp');
set(numCompH, 'string', num2str(numOfPC1));

set(handles, 'pointer', 'arrow');

function [comp_est, mdl, aic, kic] = order_selection(data, fwhm)

disp('Estimating dimension ...');

%% Remove mean
data = detrend(data, 0);

%% Arrange data based on correlation threshold
[V1, D1] = icatb_svd(data);

lam = diag(D1);

lam = sort(lam);

lam = lam(end:-1:1);

lam = lam(:)';

N = (size(data, 2));

if (exist('fwhm', 'var') && ~isempty(fwhm))
    N = N  / prod(fwhm);
    N = ceil(N);
end

%% Make eigen spectrum adjustment
tol = max(size(lam)) * eps(max(lam));

if (lam(end) < tol)
    lam(end) = [];
end

%% Correction on the ill-conditioned results (when tdim is large, some
% least significant eigenvalues become small negative numbers)
lam(real(lam) <= tol) = tol;

p = length(lam);
aic = zeros(1, p - 1);
kic = zeros(1, p - 1);
mdl = zeros(1, p - 1);
for k = 1:p-1
    LH = log(prod(lam(k+1:end).^(1/(p-k)) )/mean(lam(k+1:end)));
    mlh = 0.5*N*(p-k)*LH;
    df = 1 + 0.5*k*(2*p-k+1);
    aic(k) =  -2*mlh + 2*df;
    kic(k) =  -2*mlh + 3*df;
    mdl(k) =  -mlh + 0.5*df*log(N);
end

% Find the first local minimum of each ITC
itc = zeros(3, length(mdl));
itc(1,:) = aic;
itc(2,:) = kic;
itc(3,:) = mdl;

%% Use only mdl
dlap = squeeze(itc(end, 2:end)-itc(end, 1:end-1));
a = find(dlap > 0);
if isempty(a)
    comp_est = length(squeeze(itc(end, :)));
else
    comp_est = a(1);
end

disp(['Estimated components is found to be ', num2str(comp_est)]);


if ((comp_est < 2) || all(diff(mdl) < 0))
    comp_est = min([6, min(size(data))]);
    warning('FNC:DimensionalityEstimation', 'Estimated components is 1 or MDL function is monotonically decreasing. Using %d instead', comp_est);
end

disp('Done');

function show_params(handles)

handles_data = get(handles, 'userdata');
sesInfo = handles_data.sesInfo;

%% Show selected params
filesH = findobj(handles, 'tag', 'files');
set(filesH, 'foregroundcolor', [0, 1, 0]);
estimateH =  findobj(handles, 'tag', 'estimate_components');
set(estimateH, 'enable', 'on');

numCompH = findobj(handles, 'tag', 'numComp');
numComp = 2;
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

handles_data = get(handles, 'userdata');

set(handles, 'pointer', 'watch');

drawnow;

sesInfo = handles_data.sesInfo;

if (~isfield(sesInfo.userInput, 'dataInfo'))
    error('Data is not selected');
end

numCompH = findobj(handles, 'tag', 'numComp');
numComp = str2num(get(numCompH, 'string'));
sesInfo.userInput.numComp = numComp;

algoH = findobj(handles, 'tag', 'algorithm');
sesInfo.userInput.algorithm = get(algoH, 'value');

whichAnalysisH = findobj(handles, 'tag', 'which_analysis');
which_analysis = get(whichAnalysisH, 'value');
sesInfo.userInput.which_analysis = which_analysis;

scaleTypeH = findobj(handles, 'tag', 'scaleType');
scaleOptions = {'No Scaling', 'Z-scores'};
scaleVal = get(scaleTypeH, 'value');
scaleType = scaleOptions{scaleVal};
sesInfo.userInput.scaleType = scaleType;


reference_file = icatb_selectEntry('typeEntity', 'file', 'title', 'Select reference nifti file containing spatial maps', 'filter', '*.nii');
if ~isempty(reference_file)
    disp('Reference file is not selected for analysis. Connectogram cannot be plotted after ICA analysis.');
    
    
    reference_file_name = reference_file;
    if ~(exist(reference_file_name, 'file'))
        error([reference_file_name, ' doesn''t exist']);
    end
    [pathstr, fn, extn] = fileparts(reference_file_name);
    labels_file_name = fullfile(pathstr, [fn, '.txt']);
    if ~(exist(labels_file_name, 'file'))
        error([labels_file_name, ' doesn''t exist']);
    end
    
    
    sesInfo.userInput.reference_file = reference_file_name;
    compLabels = getCompLabels(labels_file_name);
    networkOpts = cell(length(compLabels), 2);
    for n = 1:length(compLabels)
        networkOpts{n, 1} = compLabels(n).name;
        networkOpts{n, 2} = compLabels(n).value;
    end
    sesInfo.userInput.network_summary_opts = networkOpts;
end

sesInfo.isInitialized = 0;

param_file = fullfile(sesInfo.userInput.pwd, [sesInfo.userInput.prefix, '_ica_parameter_info.mat']);
disp(['Saving session information in file ', param_file, ' ...']);
save(param_file, 'sesInfo');
disp('Done');
fprintf('\n');

set(handles, 'pointer', 'arrow');

try
    delete(handles);
catch
end

function compLabels = getCompLabels(txtFile)

fid = fopen(txtFile, 'r');
if (fid == -1)
    error(['File ', txtFile, ' cannot be opened for reading']);
end
try
    dd = textscan(fid, '%s', 'delimiter', '\t\n,', 'multipleDelimsAsOne', 1, 'whitespace', ' ');
    val = dd{1};
catch
    val = [];
end
fclose(fid);
val = val(icatb_good_cells(val));
chk = cellfun('isempty', regexp(val, '^\d+$'));

inds = find(chk == 1);

compLabels = repmat(struct('name', '', 'value', []), 1, length(inds));
for nI = 1:length(inds)
    compLabels(nI).name = val{inds(nI)};
    if (nI == length(inds))
        endT = length(val);
    else
        endT = inds(nI + 1) - 1;
    end
    dd = str2num(char(val{inds(nI) + 1:endT}));
    compLabels(nI).value = dd(:)';
end