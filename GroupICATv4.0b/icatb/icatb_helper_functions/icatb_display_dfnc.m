function icatb_display_dfnc(param_file)
%% Display Mancovan results
%

icatb_defaults;
global DFNC_ZOI_THRESH;
global DFNC_ZOI_TMAP_THRESH;
global FONT_COLOR;
global UI_FONTNAME;


%% Select dFNC file
if (~exist('param_file', 'var'))
    param_file = icatb_selectEntry('title', 'Select dFNC Parameter File', 'typeEntity', 'file', 'typeSelection', 'single', 'filter', '*dfnc.mat');
    drawnow;
    if (isempty(param_file))
        error('dFNC parameter file is not selected');
    end
end

load (param_file);

outputDir = fileparts(param_file);

if (isempty(outputDir))
    outputDir = pwd;
end

cd (outputDir);

outputDir = pwd;

dfncInfo.outputDir = outputDir;

post_process_file = fullfile(outputDir,  [dfncInfo.prefix, '_post_process.mat']);

%% Plot Correlations
load(post_process_file);

drawnow;

%% Delete a previous figure of display GUI
checkDispGUI = findobj('tag', 'display_dfnc');

if ~isempty(checkDispGUI)
    for ii = 1:length(checkDispGUI)
        delete(checkDispGUI(ii));
    end
end

% display figure
graphicsHandle = icatb_getGraphics('Display dFNC', 'displayGUI', 'display_dfnc', 'off');

set(graphicsHandle, 'CloseRequestFcn', @figCloseCallback);

% set graphics handle menu none
set(graphicsHandle, 'menubar', 'none');

htmlMenu = uimenu('parent', graphicsHandle, 'label', 'HTML Report', 'callback', ...
    {@htmlCallback, graphicsHandle});

optionsmenu = uimenu('parent', graphicsHandle, 'label', 'Options');
uimenu(optionsmenu,  'label', 'Connectogram (Cluster states)', 'callback', {@plotConnectoGram, graphicsHandle});

% offsets
xOffset = 0.03; yOffset = 0.065;

% title color
titleColor = [0 0.9 0.9];
% fonts
titleFont = 13;
axes('Parent', graphicsHandle, 'position', [0 0 1 1], 'visible', 'off');
xPos = 0.5; yPos = 0.97;
text(xPos, yPos, 'Display dFNC Results', 'color', titleColor, 'FontAngle', 'italic', 'fontweight', 'bold', ...
    'fontsize', titleFont, 'HorizontalAlignment', 'center', 'FontName', UI_FONTNAME);

% plot display button
buttonWidth = 0.2; buttonHeight = 0.05;
displayButtonPos = [0.75 yOffset buttonWidth buttonHeight];


displayButtonH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'pushbutton', ...
    'position', displayButtonPos, 'string', 'Display', 'tag', 'display_button', 'callback', ...
    {@displayButtonCallback, graphicsHandle});

popupHeight = 0.25;
popupWidth = 0.3;

promptWidth = 0.52;
promptHeight = 0.05;


promptPos =  [xOffset, yPos - 0.5*popupHeight - yOffset - 0.5*promptHeight, promptWidth, promptHeight];
textH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'text', 'position', promptPos, 'string', 'Select results to display', ...
    'HorizontalAlignment', 'center');
% horizontal alignment - center, vertical alignment - middle
icatb_wrapStaticText(textH);
promptPos = get(textH, 'position');

listPos = [promptPos(1) + promptPos(3) + xOffset, yPos - yOffset - popupHeight, popupWidth, popupHeight];
listH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'listbox', 'position', listPos, 'string', ...
    char('FNC Oscillations', 'Cluster stats', 'Meta State Analysis'), 'tag', 'select_list', 'min', 0, 'max', 1);

set(graphicsHandle, 'visible', 'on');
set(graphicsHandle, 'userdata', dfncInfo);


function displayButtonCallback(hObject, event_data, handles)
%% Display button callback
%

set(handles, 'pointer', 'watch');

dfncInfo = get(handles, 'userdata');
listH = findobj(handles, 'tag', 'select_list');
strs = cellstr(get(listH, 'string'));
val = get(listH, 'value');
selected_str = deblank(strs{val});

icatb_dfnc_results(dfncInfo, selected_str);

set(handles, 'pointer', 'arrow');


function htmlCallback(hObject, event_data, handles)
%% HTML Report Callback
%

set(handles, 'pointer', 'watch');

dfncInfo = get(handles, 'userdata');

icatb_dfnc_results_html(dfncInfo);

set(handles, 'pointer', 'arrow');



function plotConnectoGram(hObject, event_data, handles)
%% Plot connectogram
%

set(handles, 'pointer', 'watch');

dfncInfo = get(handles, 'userdata');

% Nifti file names
file_names = dfncInfo.userInput.compFiles;

% component network names and numbers
comp_network_names = cell(length(dfncInfo.userInput.comp), 2);
for nC = 1:size(comp_network_names, 1)
    comp_network_names{nC, 1} = dfncInfo.userInput.comp(nC).name;
    comp_network_names{nC, 2} = dfncInfo.userInput.comp(nC).value;
end

drawnow;

outputDir = dfncInfo.outputDir;
post_process_file = fullfile(outputDir,  [dfncInfo.prefix, '_post_process.mat']);
load(post_process_file);

num_clusters = dfncInfo.postprocess.num_clusters;


numParameters = 1;

inputText(numParameters).promptString = ['Enter state number/numbers'];
inputText(numParameters).uiType = 'edit';
inputText(numParameters).answerString = num2str(1:num_clusters);
inputText(numParameters).dataType = 'numeric';
inputText(numParameters).tag = 'states';
inputText(numParameters).enable = 'on';


numParameters = numParameters + 1;

inputText(numParameters).promptString = ['Enter connectivity threshold'];
inputText(numParameters).uiType = 'edit';
inputText(numParameters).answerString = num2str(-Inf);
inputText(numParameters).dataType = 'numeric';
inputText(numParameters).tag = 'conn_threshold';
inputText(numParameters).enable = 'on';

answer = icatb_inputDialog('inputtext', inputText, 'Title', 'Select connectogram options', 'handle_visibility',  'on', 'windowStyle', 'modal');

drawnow;

if (isempty(answer))
    error('Parameters not selected');
end

states = answer{1};
conn_threshold = answer{2};

for cluster_state_num = states
    
    % state #1 (connectivity matrix)
    C = icatb_vec2mat(clusterInfo.Call(cluster_state_num,:));
    
    disp(['State #', num2str(cluster_state_num)]);
    
    try
        % Slice view
        icatb_plot_connectogram([], comp_network_names, 'C', C, 'convert_to_zscores', 'yes', 'image_file_names', ...
            file_names, 'colorbar_label', ['State #', num2str(cluster_state_num), ' Corr'], 'slice_plane', 'sagittal', 'conn_threshold', conn_threshold);
    catch
        disp(lasterr);
    end
    
end


set(handles, 'pointer', 'arrow');


function figCloseCallback(hObject, event_data, handles)
% figure close callback

delete(hObject);
