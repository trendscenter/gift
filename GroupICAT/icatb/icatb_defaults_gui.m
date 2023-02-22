function icatb_defaults_gui
%% Open Defaults GUI
%

rehash path;

input_parameters = default_options;
[output_parameters, val] = icatb_OptionsWindow(input_parameters.inputParameters, input_parameters.defaults, 'on', 'title', 'Defaults');
if (~val)
    return;
end
results = output_parameters.results;
fnames = fieldnames(results);
defs = textread(which('icatb_defaults.m'), '%s', 'delimiter', '\n');

chkDeployed = 0;
try
    chkDeployed = isdeployed;
catch
end

if (chkDeployed)
    defs = regexprep(defs, 'function.*icatb_defaults', '');
end

defs(icatb_good_cells(defs) == 0) = {''};

for nF = 1:length(fnames)
    tmp = results.(fnames{nF});
    if (isnumeric(tmp))
        tmp = sprintf('%s = %s;', fnames{nF}, mat2str(tmp, 4));
    else
        tmp = sprintf('%s = ''%s'';', fnames{nF}, tmp);
    end
    
    inds = find(icatb_good_cells(regexp(defs, [fnames{nF}, '.*=.*;'])) == 1);
    
    if (isempty(inds))
        defs{end + 1} = tmp;
    else
        defs(inds) = regexprep(defs(inds), [fnames{nF}, '.*=.*;'], tmp);
    end
end

outDir = icatb_selectEntry('typeSelection', 'single', 'typeEntity', 'directory', 'title', 'Select a directory to save defaults...');

if (isempty(outDir))
    error('Directory is not selected');
end

fileN = fullfile(outDir, 'icatb_defaults.m');
try
    dlmwrite(fileN, char(defs), '');
    addpath(outDir);
    disp('Defaults file added to path dynamically. Save the path using pathtool for future sessions.');
    h = helpdlg('Defaults file added to path dynamically. Save the path using pathtool for future sessions.', 'About Defaults');
    waitfor(h);
    disp('');
    rehash path;
catch
    icatb_errorDialog(lasterr);
end

function input_parameters = default_options
%% Setup options that will be shown in the defaults GUI.
%
% input_parameters contains the information about user interface controls.
%

icatb_defaults;

global BG_COLOR;
global BG2_COLOR;
global FONT_COLOR;
global BUTTON_COLOR; % BUTTON background color
global BUTTON_FONT_COLOR; % BUTTON foreground color
global UI_FONTNAME;
global UI_FONTUNITS;
global UI_FS;
global FUNCTIONAL_DATA_FILTER;
global ZIP_IMAGE_FILES;
global CENTER_IMAGES;
global DEFAULT_MASK_OPTION;
global REMOVE_CONSTANT_VOXELS;
global DEFAULT_MASK_SBM_MULTIPLIER;
global SPM_STATS_WRITE_TAL;
global SPM_STATS_TTEST2_EXPLICIT_MASK;
global SPM_STATS_TTEST_THRESHOLD;
global SPM_STATS_AVG_RUNS;
global NPOINT_FFT_GROUP_COMPARISON;
global NUM_BINS_GROUP_COMPARISON;
global DEFAULT_TR_SPECTRAL_GROUP_COMPARE;
global DETRENDNUMBER;
global MAX_AVAILABLE_RAM;
global CONSERVE_DISK_SPACE;
global WRITE_ANALYSIS_STEPS_IN_DIRS;

%%%%%%%%% Input Parameters Structure
numParameters = 1;
inputParameters(numParameters).listString = 'Colors';
optionNumber = 1;
% Option 1 of parameter 1
options(optionNumber).promptString = 'Figure Background Color';
options(optionNumber).answerString = mat2str(BG_COLOR, 4);
options(optionNumber).uiType = 'edit'; options(optionNumber).value = 1;
options(optionNumber).tag = 'BG_COLOR'; options(optionNumber).answerType = 'numeric';
options(optionNumber).flag = 'delete';
options(optionNumber).callback = {@colorsCallback};
options(optionNumber).uiPos = [];
options(optionNumber).contextmenu = struct('label', 'Select Color', 'callback', @openColorsDialog);

optionNumber = optionNumber + 1;
% Option 2 of parameter 1
options(optionNumber).promptString = 'UI Background Color';
options(optionNumber).answerString = mat2str(BG2_COLOR, 4);
options(optionNumber).uiType = 'edit'; options(optionNumber).value = 1;
options(optionNumber).tag = 'BG2_COLOR'; options(optionNumber).answerType = 'numeric';
options(optionNumber).flag = 'delete';
options(optionNumber).callback = {@colorsCallback};
options(optionNumber).uiPos = [];
options(optionNumber).contextmenu = struct('label', 'Select Color', 'callback', @openColorsDialog);

optionNumber = optionNumber + 1;
% Option 3 of parameter 1
options(optionNumber).promptString = 'UI Foreground Color';
options(optionNumber).answerString = mat2str(FONT_COLOR, 4);
options(optionNumber).uiType = 'edit'; options(optionNumber).value = 1;
options(optionNumber).tag = 'FONT_COLOR'; options(optionNumber).answerType = 'numeric';
options(optionNumber).flag = 'delete';
options(optionNumber).callback = {@colorsCallback};
options(optionNumber).uiPos = [];
options(optionNumber).contextmenu = struct('label', 'Select Color', 'callback', @openColorsDialog);


optionNumber = optionNumber + 1;
% Option 4 of parameter 1
options(optionNumber).promptString = 'Button Background Color';
options(optionNumber).answerString = mat2str(BUTTON_COLOR, 4);
options(optionNumber).uiType = 'edit'; options(optionNumber).value = 1;
options(optionNumber).tag = 'BUTTON_COLOR'; options(optionNumber).answerType = 'numeric';
options(optionNumber).flag = 'delete';
options(optionNumber).callback = {@colorsCallback};
options(optionNumber).uiPos = [];
options(optionNumber).contextmenu = struct('label', 'Select Color', 'callback', @openColorsDialog);

optionNumber = optionNumber + 1;
% Option 5 of parameter 1
options(optionNumber).promptString = 'Button Foreground Color';
options(optionNumber).answerString = mat2str(BUTTON_FONT_COLOR, 4);
options(optionNumber).uiType = 'edit'; options(optionNumber).value = 1;
options(optionNumber).tag = 'BUTTON_FONT_COLOR'; options(optionNumber).answerType = 'numeric';
options(optionNumber).flag = 'delete';
options(optionNumber).callback = {@colorsCallback};
options(optionNumber).uiPos = [];
options(optionNumber).contextmenu = struct('label', 'Select Color', 'callback', @openColorsDialog);

inputParameters(numParameters).options = options;
clear options;

numParameters = numParameters + 1;
inputParameters(numParameters).listString = 'Fonts';

optionNumber = 1;
% Option 1 of parameter 2
val = strmatch(lower(UI_FONTNAME), lower(listfonts), 'exact');
if (isempty(val))
    val = 1;
end
options(optionNumber).promptString = 'Font Name';
options(optionNumber).answerString = char(listfonts);
options(optionNumber).uiType = 'popup'; options(optionNumber).value = val;
options(optionNumber).tag = 'UI_FONTNAME'; options(optionNumber).answerType = 'string';
options(optionNumber).flag = 'delete';
options(optionNumber).uiPos = [];


optionNumber = optionNumber + 1;
fontUnits = char('points', 'normalized', 'inches', 'centimeters', 'pixels');
val = strmatch(lower(UI_FONTUNITS), lower(fontUnits), 'exact');
if (isempty(val))
    val = 1;
end
% Option 2 of parameter 2
options(optionNumber).promptString = 'Font Units';
options(optionNumber).answerString = fontUnits;
options(optionNumber).uiType = 'popup'; options(optionNumber).value = val;
options(optionNumber).tag = 'UI_FONTUNITS'; options(optionNumber).answerType = 'string';
options(optionNumber).flag = 'delete';
options(optionNumber).uiPos = [];

% Option 3 of parameter 2
optionNumber = optionNumber + 1;
options(optionNumber).promptString = 'Font Size';
options(optionNumber).answerString = mat2str(UI_FS, 4);
options(optionNumber).uiType = 'edit';
options(optionNumber).tag = 'UI_FS'; options(optionNumber).answerType = 'numeric';
options(optionNumber).flag = 'delete';
options(optionNumber).uiPos = [0.1, 0.05];

inputParameters(numParameters).options = options;
clear options;


numParameters = numParameters + 1;
inputParameters(numParameters).listString = 'Component Image Defaults';
optionNumber = 1;
opts = char('*.img', '*.nii');
% Option 1 of parameter 3
val = strmatch(lower(FUNCTIONAL_DATA_FILTER), lower(opts), 'exact');
if (isempty(val))
    val = 1;
end
options(optionNumber).promptString = 'Image Format';
options(optionNumber).answerString = opts;
options(optionNumber).uiType = 'popup'; options(optionNumber).value = val;
options(optionNumber).tag = 'FUNCTIONAL_DATA_FILTER'; options(optionNumber).answerType = 'string';
options(optionNumber).flag = 'delete';
options(optionNumber).uiPos = [0.1, 0.05];

optionNumber = optionNumber + 1;
opts = char('Yes', 'No');
% Option 1 of parameter 2
val = strmatch(lower(ZIP_IMAGE_FILES), lower(opts), 'exact');
if (isempty(val))
    val = 1;
end
options(optionNumber).promptString = 'Compress Images?';
options(optionNumber).answerString = opts;
options(optionNumber).uiType = 'popup'; options(optionNumber).value = val;
options(optionNumber).tag = 'ZIP_IMAGE_FILES'; options(optionNumber).answerType = 'string';
options(optionNumber).flag = 'delete';
options(optionNumber).uiPos = [0.1, 0.05];

optionNumber = optionNumber + 1;
opts = char('0', '1');
% Option 1 of parameter 3
val = strmatch(mat2str(CENTER_IMAGES), lower(opts), 'exact');
if (isempty(val))
    val = 1;
end
options(optionNumber).promptString = 'Center Images?';
options(optionNumber).answerString = opts;
options(optionNumber).uiType = 'popup'; options(optionNumber).value = val;
options(optionNumber).tag = 'CENTER_IMAGES'; options(optionNumber).answerType = 'numeric';
options(optionNumber).flag = 'delete';
options(optionNumber).uiPos = [0.1, 0.05];

inputParameters(numParameters).options = options;
clear options;

% Option 1 of parameter 4
opts = char('first_file', 'all_files');
val = strmatch(lower(DEFAULT_MASK_OPTION), lower(opts), 'exact');
if (isempty(val))
    val = 1;
end
numParameters = numParameters + 1;
inputParameters(numParameters).listString = 'Mask';
optionNumber = 1;
options(optionNumber).promptString = 'Default Mask Option';
options(optionNumber).answerString = opts;
options(optionNumber).uiType = 'popup'; options(optionNumber).value = val;
options(optionNumber).tag = 'DEFAULT_MASK_OPTION'; options(optionNumber).answerType = 'string';
options(optionNumber).flag = 'delete';
options(optionNumber).uiPos = [];

opts = char('0', '1');
% Option 2 of parameter 4
val = strmatch(mat2str(REMOVE_CONSTANT_VOXELS), lower(opts), 'exact');
if (isempty(val))
    val = 1;
end

optionNumber = optionNumber + 1;
options(optionNumber).promptString = 'Remove Constant Voxels?';
options(optionNumber).answerString = opts;
options(optionNumber).uiType = 'popup'; options(optionNumber).value = val;
options(optionNumber).tag = 'REMOVE_CONSTANT_VOXELS'; options(optionNumber).answerType = 'numeric';
options(optionNumber).flag = 'delete';
options(optionNumber).uiPos = [0.1, 0.05];

optionNumber = optionNumber + 1;
options(optionNumber).promptString = 'Default Mask Multiplier(SBM)';
options(optionNumber).answerString = mat2str(DEFAULT_MASK_SBM_MULTIPLIER);
options(optionNumber).uiType = 'edit'; options(optionNumber).value = 1;
options(optionNumber).tag = 'DEFAULT_MASK_SBM_MULTIPLIER'; options(optionNumber).answerType = 'numeric';
options(optionNumber).flag = 'delete';
options(optionNumber).uiPos = [0.1, 0.05];


inputParameters(numParameters).options = options;
clear options;

opts = char('0', '1', '2');
% Option 1 of parameter 5
val = strmatch(mat2str(SPM_STATS_WRITE_TAL), lower(opts), 'exact');
if (isempty(val))
    val = 1;
end
numParameters = numParameters + 1;
inputParameters(numParameters).listString = 'SPM Stats';
optionNumber = 1;
options(optionNumber).promptString = 'SPM Stats and Write Talairach?';
options(optionNumber).answerString = opts;
options(optionNumber).uiType = 'popup'; options(optionNumber).value = val;
options(optionNumber).tag = 'SPM_STATS_WRITE_TAL'; options(optionNumber).answerType = 'numeric';
options(optionNumber).flag = 'delete';
options(optionNumber).uiPos = [0.1, 0.05];

opts = char('0', '1');
% Option 2 of parameter 5
val = strmatch(mat2str(SPM_STATS_TTEST2_EXPLICIT_MASK), lower(opts), 'exact');
if (isempty(val))
    val = 1;
end
optionNumber = optionNumber + 1;
options(optionNumber).promptString = 'Explicit Mask?';
options(optionNumber).answerString = opts;
options(optionNumber).uiType = 'popup'; options(optionNumber).value = val;
options(optionNumber).tag = 'SPM_STATS_TTEST2_EXPLICIT_MASK'; options(optionNumber).answerType = 'numeric';
options(optionNumber).flag = 'delete';
options(optionNumber).uiPos = [0.1, 0.05];

optionNumber = optionNumber + 1;
options(optionNumber).promptString = 'T-test Threshold';
options(optionNumber).answerString = mat2str(SPM_STATS_TTEST_THRESHOLD);
options(optionNumber).uiType = 'edit'; options(optionNumber).value = 1;
options(optionNumber).tag = 'SPM_STATS_TTEST_THRESHOLD'; options(optionNumber).answerType = 'numeric';
options(optionNumber).flag = 'delete';
options(optionNumber).uiPos = [0.1, 0.05];

opts = char('0', '1');
% Option 2 of parameter 4
val = strmatch(mat2str(SPM_STATS_AVG_RUNS), lower(opts), 'exact');
if (isempty(val))
    val = 1;
end
optionNumber = optionNumber + 1;
options(optionNumber).promptString = 'Average Runs/sessions?';
options(optionNumber).answerString = opts;
options(optionNumber).uiType = 'popup'; options(optionNumber).value = val;
options(optionNumber).tag = 'SPM_STATS_AVG_RUNS'; options(optionNumber).answerType = 'numeric';
options(optionNumber).flag = 'delete';
options(optionNumber).uiPos = [0.1, 0.05];

inputParameters(numParameters).options = options;
clear options;

numParameters = numParameters + 1;
inputParameters(numParameters).listString = 'FFT';
optionNumber = 1;
options(optionNumber).promptString = 'N Point FFT';
options(optionNumber).answerString = mat2str(NPOINT_FFT_GROUP_COMPARISON);
options(optionNumber).uiType = 'edit'; options(optionNumber).value = 1;
options(optionNumber).tag = 'NPOINT_FFT_GROUP_COMPARISON'; options(optionNumber).answerType = 'numeric';
options(optionNumber).flag = 'delete';
options(optionNumber).uiPos = [0.1, 0.05];

optionNumber = optionNumber + 1;
options(optionNumber).promptString = 'No. Of Bins';
options(optionNumber).answerString = mat2str(NUM_BINS_GROUP_COMPARISON);
options(optionNumber).uiType = 'edit'; options(optionNumber).value = 1;
options(optionNumber).tag = 'NUM_BINS_GROUP_COMPARISON'; options(optionNumber).answerType = 'numeric';
options(optionNumber).flag = 'delete';
options(optionNumber).uiPos = [0.1, 0.05];

optionNumber = optionNumber + 1;
options(optionNumber).promptString = 'TR in sec';
options(optionNumber).answerString = mat2str(DEFAULT_TR_SPECTRAL_GROUP_COMPARE);
options(optionNumber).uiType = 'edit'; options(optionNumber).value = 1;
options(optionNumber).tag = 'DEFAULT_TR_SPECTRAL_GROUP_COMPARE'; options(optionNumber).answerType = 'numeric';
options(optionNumber).flag = 'delete';
options(optionNumber).uiPos = [0.1, 0.05];

inputParameters(numParameters).options = options;
clear options;

opts = char('0', '1', '2', '3');
% Option 2 of parameter 4
val = strmatch(mat2str(DETRENDNUMBER), lower(opts), 'exact');
if (isempty(val))
    val = 1;
end
numParameters = numParameters + 1;
inputParameters(numParameters).listString = 'Other';
optionNumber = 1;
options(optionNumber).promptString = 'Detrend Level';
options(optionNumber).answerString = opts;
options(optionNumber).uiType = 'popup'; options(optionNumber).value = val;
options(optionNumber).tag = 'DETRENDNUMBER'; options(optionNumber).answerType = 'numeric';
options(optionNumber).flag = 'delete';
options(optionNumber).uiPos = [0.1, 0.05];

optionNumber = optionNumber + 1;
options(optionNumber).promptString = 'Max RAM Available in GB?';
options(optionNumber).answerString = mat2str(MAX_AVAILABLE_RAM);
options(optionNumber).uiType = 'edit'; options(optionNumber).value = 1;
options(optionNumber).tag = 'MAX_AVAILABLE_RAM'; options(optionNumber).answerType = 'numeric';
options(optionNumber).flag = 'delete';
options(optionNumber).uiPos = [0.1, 0.05];

opts = char('0', '1');
% Option 2 of parameter 4
val = strmatch(mat2str(WRITE_ANALYSIS_STEPS_IN_DIRS), lower(opts), 'exact');
if (isempty(val))
    val = 1;
end
optionNumber = optionNumber + 1;
options(optionNumber).promptString = 'Write Analysis Steps in Directories?';
options(optionNumber).answerString = opts;
options(optionNumber).uiType = 'popup'; options(optionNumber).value = val;
options(optionNumber).tag = 'WRITE_ANALYSIS_STEPS_IN_DIRS'; options(optionNumber).answerType = 'numeric';
options(optionNumber).flag = 'delete';
options(optionNumber).uiPos = [0.1, 0.05];


opts = char('0', '1', '2');
% Option 2 of parameter 4
val = strmatch(mat2str(CONSERVE_DISK_SPACE), lower(opts), 'exact');
if (isempty(val))
    val = 1;
end
optionNumber = optionNumber + 1;
options(optionNumber).promptString = 'Conserve Disk Space?';
options(optionNumber).answerString = opts;
options(optionNumber).uiType = 'popup'; options(optionNumber).value = val;
options(optionNumber).tag = 'CONSERVE_DISK_SPACE'; options(optionNumber).answerType = 'numeric';
options(optionNumber).flag = 'delete';
options(optionNumber).uiPos = [0.1, 0.05];

inputParameters(numParameters).options = options;
clear options;

% store the input parameters in two fields
input_parameters.inputParameters = inputParameters;
input_parameters.defaults = inputParameters;


function colorsCallback(hObject, event_data, handles)
%% Set colors callback
%

tag = get(hObject, 'tag');
tag = [tag, 'ContextMenu'];

if strcmpi(get(get(hObject, 'parent'), 'SelectionType'), 'alternate')
    if (~isempty(findobj(handles, 'tag', tag)))
        hcmenu = uicontextmenu;
        uimenu(hcmenu, 'Label', 'Select Colors', 'Callback', {@openColorsDialog, hObject}, 'tag', tag);
    end
end

function openColorsDialog(hObject, event_data, handles)
%% Open uisetcolor
%

tag = get(hObject, 'tag');
tag = strrep(tag, 'ContextMenu', '');

hd = findobj(get(get(hObject, 'parent'), 'parent'), 'tag', tag);
ud = get(hd, 'string');
try
    ud = uisetcolor(str2num(ud));
    set(hd, 'string', mat2str(ud, 4));
catch
end