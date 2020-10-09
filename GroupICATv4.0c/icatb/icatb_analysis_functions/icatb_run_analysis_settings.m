function run_analysis_settings = icatb_run_analysis_settings(sesInfo)
%% Run analysis settings
%

icatb_defaults;
global BG_COLOR;
global FONT_COLOR;
global HELP_FONT_COLOR;

checkH = findobj('tag', 'run_analysis_window');

if (ishandle(checkH))
    try
        delete(checkH);
    catch
    end
end

handleFig = icatb_getGraphics('Run Analysis', 'normal', 'run_analysis_window', 'off');
set(handleFig, 'menubar', 'none');

%% Offsets
yOffset = 0.05;
xOffset = 0.05;


%% Plot push buttons
buttonHeight = 0.05;
buttonWidth = 0.12;

helpPos = [0.75 - 0.5*buttonWidth, yOffset + 0.5*buttonHeight, buttonWidth, buttonHeight];
helpH = icatb_uicontrol('parent', handleFig, 'units', 'normalized', 'style', 'pushbutton', 'position', ...,
    helpPos, 'string', '?', 'tag', 'help_button', 'foregroundcolor', [1, 1, 0], 'horizontalalignment', 'center', 'fontweight', 'bold', 'callback', {@helpCallback, handleFig});

donePos = [0.5 - 0.5*buttonWidth, yOffset + 0.5*buttonHeight, buttonWidth, buttonHeight];
doneH = icatb_uicontrol('parent', handleFig, 'units', 'normalized', 'style', 'pushbutton', 'position', ...,
    donePos, 'string', 'Done', 'tag', 'done_button', 'horizontalalignment', 'center', 'callback', {@doneCallback, handleFig});

cancelPos = [0.25 - 0.5*buttonWidth, yOffset + 0.5*buttonHeight, buttonWidth, buttonHeight];
cancelH = icatb_uicontrol('parent', handleFig, 'units', 'normalized', 'style', 'pushbutton', 'position', ...,
    cancelPos, 'string', 'Cancel', 'tag', 'cancel_button', 'horizontalalignment', 'center', 'callback', 'delete(gcbf)');

remSpace = (1 - cancelPos(2) - cancelPos(4) - 3*yOffset);

opts = icatb_get_analysis_settings;

opts = cellstr(str2mat(opts));

showPerfsettings = 'on';

if (sesInfo.userInput.numOfSub*sesInfo.userInput.numOfSess > 1)
    
    max_mem = zeros(length(opts), 1);
    for nOpts = 1:length(opts)
        max_mem(nOpts) = icatb_get_analysis_settings(sesInfo, opts{nOpts});
    end
    
    opts = strcat(opts, repmat({' ('}, length(opts), 1), cellstr(num2str(max_mem, '%0.4f')), repmat({' GB)'}, length(opts), 1));
    
end


icaAlgo = icatb_icaAlgorithm; % available ICA algorithms

algoVal = sesInfo.userInput.algorithm; % algorithm index

% selected ICA algorithm
algorithmName = deblank(icaAlgo(algoVal, :));

if strcmpi(algorithmName, 'moo-icar')
    algorithmName = 'gig-ica';
end

if ((sesInfo.userInput.numOfSub*sesInfo.userInput.numOfSess == 1) || strcmpi(algorithmName, 'gig-ica') || strcmpi(algorithmName, 'constrained ica (spatial)'))
    showPerfsettings = 'off';
end

frameWidth = (1 - 2*xOffset);
frameHeight = 0.5*remSpace;

%% Frame comprising radiobuttons
framePos = [xOffset, cancelPos(2) + cancelPos(4) + yOffset, frameWidth, frameHeight];
frameH = icatb_uicontrol('parent', handleFig, 'units', 'normalized', 'style', 'frame', 'position', framePos, 'tag', 'frame_radiobuttons', ...
    'backgroundcolor', BG_COLOR, 'foregroundcolor', FONT_COLOR, 'visible', showPerfsettings, 'enable', showPerfsettings);

text1Pos = framePos + xOffset;
textHeight = 0.05;
text1Pos(3) = framePos(3) - 2*xOffset;
text1Pos(2) = framePos(2) + framePos(4) - textHeight;
text1Pos(4) = textHeight;

text1H = icatb_uicontrol('parent', handleFig, 'units', 'normalized', 'String', 'Group PCA Performance Settings', 'horizontalalignment', 'center', 'style', 'text', ...
    'position', text1Pos, 'backgroundcolor', BG_COLOR, 'foregroundcolor', FONT_COLOR, 'visible', showPerfsettings, 'enable', showPerfsettings);
icatb_wrapStaticText(text1H);
tmpPos = get(text1H, 'extent');
tmpPos = tmpPos;
text1Pos(2) = framePos(2) + framePos(4) - tmpPos(4);
text1Pos(3) = tmpPos(3);
text1Pos(1) = framePos(1) + (framePos(1) + framePos(3) - text1Pos(3))/2;
set(text1H, 'position', text1Pos);

radioHeight = (frameHeight - 2*yOffset)/3;
radioWidth = (framePos(1) + framePos(3) - 3*xOffset);

%% Plot Radio buttons
radioPos = [framePos(1) + xOffset, framePos(2) + yOffset, radioWidth, radioHeight];

radio3H = icatb_uicontrol('parent', handleFig, 'units', 'normalized', 'position', radioPos, 'style', 'radiobutton', 'string', opts{end}, 'tag', 'analysisType3', ...
    'min', 0, 'max', 1, 'callback', {@radioButtonCallback, handleFig, '^analysisType'}, 'value', 1, 'visible', showPerfsettings, 'enable', showPerfsettings);

radioPos(2) = radioPos(2) + radioHeight;

radio2H = icatb_uicontrol('parent', handleFig, 'units', 'normalized', 'position', radioPos, 'style', 'radiobutton', 'string', opts{end-1}, 'tag', 'analysisType2', ...
    'min', 0, 'max', 1, 'callback', {@radioButtonCallback, handleFig, '^analysisType'}, 'value', 0, 'visible', showPerfsettings, 'enable', showPerfsettings);

radioPos(2) = radioPos(2) + radioHeight;
radio1H = icatb_uicontrol('parent', handleFig, 'units', 'normalized', 'position', radioPos, 'style', 'radiobutton', 'string', opts{end-2}, 'tag', 'analysisType1', ...
    'min', 0, 'max', 1, 'callback', {@radioButtonCallback, handleFig, '^analysisType'}, 'value', 0, 'visible', showPerfsettings, 'enable', showPerfsettings);

%% Draw text box and Listboxes
textWidth = 0.35;
textHeight = 0.05;

listboxOrigin = framePos(2) + framePos(4) + yOffset;
listHeight = (1 - yOffset - listboxOrigin);
listWidth = 1 - textWidth - 3*xOffset;
listPos = [2*xOffset + textWidth, listboxOrigin, listWidth, listHeight];

textOrigin = framePos(2) + framePos(4) + yOffset + 0.5*remSpace;
text2Pos = [xOffset, listboxOrigin + 0.5*listHeight - 0.5*textHeight, textWidth, textHeight];

text2H = icatb_uicontrol('parent', handleFig, 'units', 'normalized', 'style', 'text', 'position', text2Pos, 'string', 'Select Analysis Step/Steps');
icatb_wrapStaticText(text2H);
text2Pos = get(text2H, 'position');
text2Pos(2) = listPos(2) + 0.5*listPos(4) - 0.5*text2Pos(4);
set(text2H, 'position', text2Pos);


analysisStr = {'All***', 'Resume', 'Initialize Parameters', 'Group Data Reduction', 'Calculate ICA/IVA', ...
    'Back Reconstruct', 'Calibrate Components', 'Group Stats'};
userdata = {'all', 'resume', 'parameter_initialization', 'group_pca', 'calculate_ica', 'back_reconstruct', 'scale_components', 'group_stats'};

conserve_disk_space = 0;
if (isfield(sesInfo.userInput, 'conserve_disk_space'))
    conserve_disk_space = sesInfo.userInput.conserve_disk_space;
end

if (conserve_disk_space == 1)
    analysisStr([6, 8]) = [];
    userdata([6, 8]) = [];
end

if (strcmpi(algorithmName, 'gig-ica') || strcmpi(algorithmName, 'constrained ica (spatial)'))
    
    chkSpInds = strcmpi(analysisStr, 'back reconstruct');
    analysisStr(chkSpInds) = [];
    userdata(chkSpInds) = [];
    
    chkSpInds = strcmpi(analysisStr, 'group data reduction');
    analysisStr(chkSpInds) = [];
    userdata(chkSpInds) = [];
    
end

listH = icatb_uicontrol('parent', handleFig, 'units', 'normalized', 'style', 'listbox', 'position', listPos, 'string', analysisStr, 'max', 2, 'min', 0, 'value', 1, ...
    'userdata', userdata, 'tag', 'analysis_steps');

set(handleFig, 'visible', 'on');

try
    waitfor(handleFig);
catch
end

drawnow;

%% Return output
if (isappdata(0, 'run_analysis_settings_data'))
    run_analysis_settings = getappdata(0, 'run_analysis_settings_data');
    rmappdata(0, 'run_analysis_settings_data');
else
    error('Run analysis step/steps are not selected');
end


function radioButtonCallback(hObject, event_data, handles, matchstr)
%% Radio button callback

radioHandles = findobj(handles, 'style', 'radiobutton');

allTags = get(radioHandles, 'tag');

check = regexp(allTags, matchstr);
check = icatb_good_cells(check);

allTags = allTags(check);
radioHandles = radioHandles(check);

isRadio = strcmp(allTags, get(hObject, 'tag'));

radioHandles(isRadio) = [];

set(radioHandles, 'value', 0);


function doneCallback(hObject, event_data, handles)
%% Done callback

opts = icatb_get_analysis_settings;
radioVals = [get(findobj(handles, 'tag', 'analysisType1'), 'value'), get(findobj(handles, 'tag', 'analysisType2'), 'value'), get(findobj(handles, 'tag', 'analysisType3'), 'value')];

if (all(radioVals == 0))
    error('Please select the appropriate group PCA performance settings');
end

perfType = lower(opts{radioVals == 1});

listH = findobj(handles, 'tag', 'analysis_steps');
listVal = get(listH, 'value');
userdata = get(listH, 'userdata');

if (isempty(listVal))
    error('Analysis step is not selected');
end

delete(handles);

drawnow;

stepsToRun = userdata(listVal);

%% Set Application data
run_analysis_settings_data.perfType = perfType;
run_analysis_settings_data.stepsToRun = lower(stepsToRun);
setappdata(0, 'run_analysis_settings_data', run_analysis_settings_data);


function helpCallback(hObject, event_data, handles)
%% Open help dialog box
%

msgString = str2mat('I. Select Analysis Step/Steps - There are 8 options available: ', ...
    '', ...
    '  a. All*** - Run all steps including Parameter Initialization, Data Reduction, Calculate ICA, Back Reconstruction, Calibrate Components and Group Stats.', ...
    '', ...
    '  b. Resume - Resume an incomplete analysis.', ...
    '', ...
    '  c. Parameter Initialization - Parameters initialization step parses user input and does error check before running group ICA.', ...
    '', ...
    '  d. Data Reduction - Subjects data is reduced to a few no. of components and passed to group data reduction step for multi-subject analysis.', ...
    '', ...
    '  e. Calculate ICA - ICA is run on the final reduced data from the data reduction step. ', ...
    '', ...
    '  f. Back Reconstruction - Individual subject components are reconstructed using the ICA information', ...
    '', ...
    '  g. Calibrate Components - Component maps and time courses could be scaled to either percent signal change or z-scores.', ...
    '', ...
    '  h. Group Stats - Group statistics are done on the individual subject components ', ...
    '', ...
    '', ...
    'II. Group PCA Performance Settings - There are three options available. The best match for each option is dependent on the MAX_AVAILABLE_MEM variable in icatb_defaults.', ...
    '', ...
    '  a. Maximize Performance - Data from subjects are blindly stacked. Covariance based PCA or expectation maximization is used.', ...
    '', ...
    '  b. Less Memory Usage - Typically data from subjects will not be stacked. For large data-sets covariance matrix size might be larger than stacked data-sets. In such cases expectation maximization is used.', ...
    '', ...
    '  c. User Specified settings - User settings will be selected.');

h = icatb_dialogBox('title', 'Run Analysis Settings', 'textBody', msgString, 'textType', 'large');

try
    waitfor(h);
catch
end