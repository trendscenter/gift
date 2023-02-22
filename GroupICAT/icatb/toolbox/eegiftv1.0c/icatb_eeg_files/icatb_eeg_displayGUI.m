function icatb_eeg_displayGUI(paramFile)
% Display GUI for displaying EEG components

% Load defaults
icatb_defaults;

global PARAMETER_INFO_MAT_FILE;
global UI_FONTNAME;

% load parameters file
filterP = ['*', PARAMETER_INFO_MAT_FILE, '*.mat'];
if ~exist('paramFile', 'var')
    [paramFile] = icatb_selectEntry('typeEntity', 'file', 'title', 'Select Parameter File', ...
        'filter', filterP);
end

drawnow;

[pathstr, fileName] = fileparts(paramFile);

cd(pathstr); % Cd to parameter file location
load(paramFile);

if ~exist('sesInfo', 'var')
    error('Please select the ICA parameter file');
end

if ~sesInfo.isInitialized
    error('Please run the analysis to display the results');
end


if isfield(sesInfo, 'modality')
    if ~strcmpi(sesInfo.modality, 'eeg')
        error('Use GIFT toolbox to display ICA results');
    end
end

displayGUIMsg = 'Opening EEG Display GUI. Please wait ...';
disp(displayGUIMsg);
helpHandle = helpdlg(displayGUIMsg, 'Opening EEG Display GUI');

zipContents.zipFiles = {};
zipContents.files_in_zip(1).name = {};
if isfield(sesInfo, 'zipContents')
    zipContents = sesInfo.zipContents;
end

[subjectICAFiles, calibratedMATFiles] = getSubjectMatFileNames(sesInfo);

% Get results from sesInfo file
dispParameters.icaOutputFiles = sesInfo.icaOutputFiles;
dispParameters.subjectICAFiles = subjectICAFiles;
dispParameters.calibratedMATFiles = calibratedMATFiles;
dispParameters.numOfSess = sesInfo.numOfSess;
dispParameters.numOfSub = sesInfo.numOfSub;
dispParameters.numOfComp = sesInfo.numComp;
dispParameters.inputFiles = sesInfo.inputFiles;
dispParameters.HInfo = sesInfo.HInfo;
dispParameters.outputDir = pathstr; % save the output directory
dispParameters.paramFile = paramFile; % parameter file
dispParameters.inputPrefix = sesInfo.userInput.prefix;
dispParameters.zipContents = zipContents;
dispParameters.mask_ind = sesInfo.mask_ind;
dispParameters.diffTimePoints = sesInfo.diffTimePoints;

if length(find(dispParameters.diffTimePoints == dispParameters.diffTimePoints(1))) ~= ...
        length(dispParameters.diffTimePoints)
    error('Cannot display EEG components as the number of images must be the same between subjects');
end

clear sesInfo;

% store the information about the structural file
%[dispParameters] = icatb_getStructuralFile(dispParameters);

% delete a previous figure of display GUI
checkDispGUI = findobj('tag', 'display_gui_eeg');

if ~isempty(checkDispGUI)
    for ii = 1:length(checkDispGUI)
        delete(checkDispGUI(ii));
    end
end

% display figure
graphicsHandle = icatb_getGraphics('EEG Display GUI', 'displayGUI', 'display_gui_eeg', 'off');

% set graphics handle menu none
set(graphicsHandle, 'menubar', 'none');

% plot display defaults
display_defaultsMenu = uimenu('parent', graphicsHandle, 'label', 'Display Defaults', 'callback', ...
    {@display_defaults_callback, graphicsHandle});


% EEG Display GUI Options Menu
dispGUIOptionsMenu = uimenu('parent', graphicsHandle, 'label', 'Options');

% Channel Locations Menu
channelLocsH = uimenu(dispGUIOptionsMenu, 'label', 'Select EEG channel locs file', 'callback', ...
    {@locsFileCallback, graphicsHandle});

% Time Scale Menu
timeScaleH = uimenu(dispGUIOptionsMenu, 'label', 'Select time scale', 'callback', ...
    {@timeScaleCallback, graphicsHandle});

% Help Menu
helpMenu = uimenu('parent', graphicsHandle, 'label', 'EEGIFT-Help');
htmlHelpMenu = uimenu(helpMenu, 'label', 'Display GUI', 'callback', 'icatb_openHTMLHelpFile(''icatb_eeg_displayGUI.htm'');');

% offsets
xOffset = 0.05; yOffset = 0.05;

%%%%%%%%%%%%% Draw Title here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% title color
titleColor = [0 0.9 0.9];
% fonts
titleFont = 13;
axes('Parent', graphicsHandle, 'position', [0 0 1 1], 'visible', 'off');
xPos = 0.5; yPos = 0.97;
text(xPos, yPos, 'EEG Display GUI', 'color', titleColor, 'FontAngle', 'italic', 'fontweight', 'bold', ...
    'fontsize', titleFont, 'HorizontalAlignment', 'center', 'FontName', UI_FONTNAME);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot controls in two frames
% Upper frame: Images
% Lower frame: UIcontrols

% plot display button
buttonWidth = 0.2; buttonHeight = 0.05;
displayButtonPos = [0.75 yOffset buttonWidth buttonHeight];


displayButtonH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'pushbutton', ...
    'position', displayButtonPos, 'string', 'Display', 'tag', 'display_button', 'callback', ...
    {@displayButtonCallback, graphicsHandle});


% Plot sort components followed by listboxes for viewing set and
% component number

% plot sort components here
textHeight = 0.05; textWidth = 0.35;
control_width = 0.35; control_height = 0.05;

% text origin
textYOrigin = 0.9 - 0.5*textHeight;

textPos = [xOffset, textYOrigin, textWidth, textHeight];

% plot text
textH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'text', 'position', ...
    textPos, 'string', 'Select display method', 'HorizontalAlignment', 'center');

% horizontal alignment - center, vertical alignment - middle
align(textH,'center','middle');

popupPos = [textPos(1) + textPos(3) + 2*xOffset, textPos(2), ...
    control_width, control_height];

clear options;

popupH = icatb_getUIPopUp(graphicsHandle, str2mat('Component', 'Subject'), popupPos, '', 'on', 'display_method');

% Set callback for display type
set(popupH, 'callback', {@displayMethodCallback, graphicsHandle});


% plot component listbox and viewing set listbox
textPos(2) = textPos(2) - textHeight - yOffset;

textWidth = 0.25; textHeight = 0.05;
% draw listbox
listWidth = 0.5; listHeight = 0.3;

textPos = [textPos(1) + 0.5*listWidth - 0.5*textWidth textPos(2) textWidth textHeight];

% plot text
textH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'text', 'position', ...
    textPos, 'string', 'Viewing Set', 'HorizontalAlignment', 'center');

% horizontal alignment - center, vertical alignment - middle
align(textH,'center','middle');

outputFiles = chkOutputFiles(dispParameters.icaOutputFiles, dispParameters.outputDir, subjectICAFiles, calibratedMATFiles);
counter = 1;
numOfSets = length(outputFiles);
% loop over number of sets
for jj = 1 : numOfSets
    % loop over sessions
    for kk = 1 : length(outputFiles(jj).ses)
        str = deblank(outputFiles(jj).ses(kk).name(1, :));
        [pathstr fileName] = fileparts(str);
        underScoreIndex = icatb_findstr(fileName, '_');
        str = fileName(1:underScoreIndex(end) - 1);
        str = [num2str(jj),'-', num2str(kk), ' ', str];
        options(counter).str = str;
        counter = counter + 1;
    end
end

dispParameters.icaOutputFiles = outputFiles;

listboxPos = [xOffset textPos(2) - listHeight - yOffset  listWidth listHeight];
% list box position
listH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'listbox', 'position', ...
    listboxPos, 'HorizontalAlignment', 'center', 'string', str2mat(options.str), 'value', 1, 'tag', ...
    'viewing_set', 'min', 0, 'max', 1);


listWidth = 0.25;

listboxPos = [1 - xOffset - listWidth, listboxPos(2), listWidth, listHeight];
textPos = [listboxPos(1) + 0.5*listWidth - 0.5*textWidth textPos(2) textWidth textHeight];

text2H = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'text', 'position', ...
    textPos, 'string', 'Component No:', 'HorizontalAlignment', 'center');

% horizontal alignment - center, vertical alignment - middle
align(text2H, 'center', 'middle');

% number of components
numComp = dispParameters.numOfComp;

clear options;

% return component index
for ii = 1:numComp
    options(ii).str = icatb_returnFileIndex(ii);
end

% list box position
list2H = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'listbox', 'position', ...
    listboxPos, 'string', str2mat(options.str), 'HorizontalAlignment', 'center', ...
    'tag', 'component_number', 'min', 0, 'max', 2);

clear options
% store the data to handles data
%handles_data.dispParameters = dispParameters;

[dispParameters, inputText] = icatb_eeg_displayGUI_defaults('init', [], dispParameters, 'off');
%dispParameters.visualizationmethods = 'component explorer';
handles_data.dispParameters = dispParameters; handles_data.inputText = inputText;

% set the figure data
set(graphicsHandle, 'userdata', handles_data);

clear handles_data;

displayMethodCallback(popupH, [], graphicsHandle);

try
    delete(helpHandle);
catch
end

% Make the figure visible
set(graphicsHandle, 'visible', 'on');

function display_defaults_callback(hObject, event_data, handles)
% get the display defaults

handles_data = get(handles, 'userdata');

[parameters, inputText] = icatb_eeg_displayGUI_defaults('no-init', handles_data.inputText, handles_data.dispParameters);

handles_data.dispParameters = parameters;
handles_data.inputText = inputText;

set(handles, 'userdata', handles_data);


function displayButtonCallback(hObject, event_data, handles)
% display button callback

try

    % get the figure data
    handles_data = get(handles, 'userdata');
    dispParameters = handles_data.dispParameters;
    % load parameter file
    load(dispParameters.paramFile);
    % indices in the mask
    dispParameters.mask_ind = sesInfo.mask_ind;
    zipContents.zipFiles = {};
    zipContents.files_in_zip(1).name = {};
    if isfield(sesInfo, 'zipContents')
        zipContents = sesInfo.zipContents;
    end
    % zip contents
    dispParameters.zipContents = zipContents;


    if ~isfield(sesInfo, 'locsFile')
        if ~isfield(sesInfo.userInput, 'locsFile')
            sesInfo = selectLocsFile(dispParameters.paramFile);
        else
            sesInfo.locsFile = sesInfo.userInput.locsFile;
        end
    end

    locsFile = sesInfo.locsFile;

    timePoints = dispParameters.diffTimePoints;

    if (sesInfo.scaleType == 0)
        ylabelStr = 'Arbitrary Units';
    elseif (sesInfo.scaleType == 1)
        ylabelStr = 'Data Units';
    else
        ylabelStr = 'Z-scores';
    end

    EEGTimeAxis = [];

    if isfield(sesInfo.userInput, 'EEGTimeAxis')
        sesInfo.EEGTimeAxis = sesInfo.userInput.EEGTimeAxis;
    end

    if isfield(sesInfo, 'EEGTimeAxis')
        EEGTimeAxis = sesInfo.EEGTimeAxis;
        if length(EEGTimeAxis) ~= dispParameters.HInfo.DIM(1)
            error('Error:TimeAxis', 'Length of time scale (%s) does not match the data (%s)', ...
                num2str(length(EEGTimeAxis)), num2str(dispParameters.HInfo.DIM(1)));
        end
    end

    clear sesInfo;

    try
        EEGlocs = icatb_eeg_pop_chanedit([], 'load', {locsFile});
    catch
        error('Error:EEGlocFile', 'Selected file %s is not a EEG location file or doesn''t exist', locsFile);
    end


    if timePoints(1) ~= length(EEGlocs)
        error('Error:ChannelLocs', ['Please check the EEG channel locs file as no. of channel locations (%s)', ...
            '\n does not match the data (%s)'], num2str(length(EEGlocs)), num2str(timePoints(1)));
    end

    if ~isempty(EEGTimeAxis)
        xlabelStr = 'Time (ms)';
    else
        xlabelStr = 'Time Points';
    end


    % find viewing set listbox
    viewListH = findobj(handles, 'tag', 'viewing_set');
    viewingsetStr = get(viewListH, 'string');
    viewingsetVal = get(viewListH, 'value');

    % viewing set selected
    viewingset = deblank(viewingsetStr(viewingsetVal, :));

    % component listbox
    compListH = findobj(handles, 'tag', 'component_number');
    % get the component numbers
    component_numbers = get(compListH, 'value');

    if isempty(component_numbers)
        component_numbers = 1;
    end

    dispParameters.componentnumber = component_numbers;

    % store viewing set
    dispParameters.viewingset = viewingset;
    returnValue = dispParameters.returnValue;
    convertToZ = dispParameters.convertToZ;
    threshValue = dispParameters.thresholdvalue;
    images_per_figure = dispParameters.imagesperfigure;

    icatb_defaults;
    global EEG_TOPOPLOT_COLORMAP;

    displayMethod = dispParameters.displayMethod;

    structDIM = dispParameters.HInfo.DIM(1:3);

    structuralData = ones(structDIM);

    % Number of components
    numComp = dispParameters.numOfComp;

    VOX = [1 1 1];

    set(handles, 'pointer', 'watch');

    switch lower(displayMethod)
        case 'component'
            %icatb_defaults;
            global TMAP_AN3_FILE;

            % check the viewing set
            if ~isfield(dispParameters, 'viewingset')
                error('Viewing set is not selected');
            else
                if isempty(dispParameters.viewingset)
                    error('Viewing set is not selected');
                end
            end

            % load components
            dashIndex = icatb_findstr(dispParameters.viewingset, '-');
            spaceIndex = icatb_findstr(dispParameters.viewingset, ' ');
            a = dispParameters.viewingset(1 : dashIndex(1) - 1);
            b = dispParameters.viewingset(dashIndex(1) + 1 : spaceIndex(1));
            compFiles = deblank(dispParameters.icaOutputFiles(str2num(a)).ses(str2num(b)).name);

            chkSub = strmatch(compFiles, dispParameters.subjectICAFiles, 'exact');
            if (~isempty(chkSub))
                compFiles = dispParameters.calibratedMATFiles{chkSub};
            end

            %compFiles = icatb_eeg_checkOutPutFiles(P, dispParameters.inputPrefix);

            % component files
            compFiles = icatb_fullFile('directory', dispParameters.outputDir, 'files', compFiles);

            helpHandle = helpdlg('Loading components...', 'Loading components');

            % resize the image and return header Info
            [icasig, weights] = icatb_eeg_loadICAData('structData', structuralData, 'compFiles', compFiles, 'comp_numbers', component_numbers);

            try
                delete(helpHandle);
            catch
            end

            % Calculate mean ICA signal
            meanICASig = (mean(icasig, 3))';

            icasig = reshape(icasig, [length(component_numbers), prod(size(structuralData))]);

            % apply display parameters
            [icasig] = icatb_applyDispParameters(icasig, convertToZ, returnValue, threshValue, [], []);

            icasig = reshape(icasig, [length(component_numbers), size(structuralData)]);

            % Store the data in a data structure
            for nComp = 1:length(component_numbers)
                [fileIndex] = icatb_returnFileIndex(component_numbers(nComp));
                [funcImg, minICAIm, maxICAIm, minInterval, maxInterval] = icatb_make_composite(structuralData, icasig(nComp, :, :), ...
                    returnValue, 'axial', structDIM, VOX);
                plotData.comp(nComp).minICAIm = minICAIm;
                plotData.comp(nComp).maxICAIm = maxICAIm;
                plotData.comp(nComp).minInterval = minInterval;
                plotData.comp(nComp).maxInterval = maxInterval;

                axesNum = 1;
                % Topoplot
                plotData.comp(nComp).axes(axesNum).data = weights(:, nComp);
                plotData.comp(nComp).axes(axesNum).plotType = 'topoplot';
                plotData.comp(nComp).axes(axesNum).title = ['Comp ', fileIndex];

                axesNum = axesNum + 1;
                % Image
                plotData.comp(nComp).axes(axesNum).data = squeeze(funcImg);
                plotData.comp(nComp).axes(axesNum).plotType = 'image';
                plotData.comp(nComp).axes(axesNum).title = ['EEG Signal'];

                axesNum = axesNum + 1;
                % Timecourse
                plotData.comp(nComp).axes(axesNum).data =  meanICASig(:, nComp);
                plotData.comp(nComp).axes(axesNum).sem =  std(squeeze(icasig(nComp, :, :)), 0, 2) ./ sqrt(size(icasig, 3) - 1);
                plotData.comp(nComp).axes(axesNum).plotType = 'timecourse';
                plotData.comp(nComp).axes(axesNum).title = ['Mean Signal Over Trials'];

            end

            % Image colormap
            plotData.figLabel = dispParameters.viewingset;

        case 'subject'

            %icatb_defaults;
            global TMAP_AN3_FILE;

            numComp = 0;
            % get the total number of components
            for ii = 1:length(dispParameters.icaOutputFiles)
                numComp = numComp + length(dispParameters.icaOutputFiles(ii).ses);
            end


            compNum = dispParameters.componentnumber(end);
            counter = 0;
            outputDir = dispParameters.outputDir;

            % resize the image and return header Info
            helpHandle = helpdlg('Loading components...', 'Loading components');
            % loop over remaining files
            for kk = 1:length(dispParameters.icaOutputFiles)
                for ii = 1:length(dispParameters.icaOutputFiles(kk).ses)
                    counter = counter + 1;
                    P = deblank(dispParameters.icaOutputFiles(kk).ses(ii).name);
                    PFullFiles = P;
                    %PFullFiles = icatb_eeg_checkOutPutFiles(P, dispParameters.inputPrefix);
                    chkSub = strmatch(PFullFiles, dispParameters.subjectICAFiles, 'exact');
                    if (~isempty(chkSub))
                        PFullFiles = dispParameters.calibratedMATFiles{chkSub};
                    end
                    PFullFiles = icatb_fullFile('directory', dispParameters.outputDir, 'files', PFullFiles);

                    [icasig, weights] = icatb_eeg_loadICAData('structData', structuralData, 'compFiles', PFullFiles, 'comp_numbers', compNum);

                    % Calculate mean ICA signal
                    tempMeanICASig = (mean(icasig, 3))';

                    icasig = reshape(icasig, [1, prod(size(structuralData))]);

                    % apply display parameters
                    [icasig] = icatb_applyDispParameters(icasig, convertToZ, returnValue, threshValue, [], []);

                    icasig = reshape(icasig, [1, size(structuralData)]);

                    [funcImg, minICAIm, maxICAIm, minInterval, maxInterval] = icatb_make_composite(structuralData, icasig, ...
                        returnValue, 'axial', structDIM, VOX);

                    plotData.comp(counter).minICAIm = minICAIm;
                    plotData.comp(counter).maxICAIm = maxICAIm;
                    plotData.comp(counter).minInterval = minInterval;
                    plotData.comp(counter).maxInterval = maxInterval;

                    axesNum = 1;
                    % Topoplot
                    plotData.comp(counter).axes(axesNum).data = weights;
                    plotData.comp(counter).axes(axesNum).plotType = 'topoplot';

                    axesNum = axesNum + 1;
                    % Image
                    plotData.comp(counter).axes(axesNum).data = squeeze(funcImg);
                    plotData.comp(counter).axes(axesNum).plotType = 'image';
                    plotData.comp(counter).axes(axesNum).title = ['EEG Signal'];

                    axesNum = axesNum + 1;
                    % Timecourse
                    plotData.comp(counter).axes(axesNum).data =  tempMeanICASig;
                    plotData.comp(counter).axes(axesNum).sem =  std(squeeze(icasig), 0, 2) ./ sqrt(size(icasig, 3) - 1);
                    plotData.comp(counter).axes(axesNum).plotType = 'timecourse';
                    plotData.comp(counter).axes(axesNum).title = ['Mean Signal Over Trials'];


                    % end for checking the data type structure or double
                    [pathstr, fileName, extn] = fileparts(deblank(P(1, :)));
                    % image extension
                    under_score_positions = icatb_findstr(fileName, '_');
                    fileName = [fileName(1:under_score_positions(end)), icatb_returnFileIndex(compNum)];
                    plotData.comp(counter).axes(1).title = fileName;
                    clear P1 funcImg
                end

            end

            if ishandle(helpHandle)
                delete(helpHandle);
            end
            plotData.figLabel = 'Figure ';

    end

    plotData.timeAxis = EEGTimeAxis;
    plotData.timeAxisXlabel = xlabelStr;
    plotData.timeAxisYlabel = ylabelStr;

    image_cmap = icatb_getColormap(1, returnValue, 1, 'other');
    nPoints = size(image_cmap, 1);
    image_cmap = image_cmap(1:ceil(nPoints/2), :);
    image_cmap = [image_cmap; ones(ceil(nPoints/2), 3)];
    topo_cmap = EEG_TOPOPLOT_COLORMAP;

    icatb_eeg_display_comp('plot_data', plotData, 'title_color', 'c', 'color_map', image_cmap, 'time_course_color', 'm', 'number_per_figure', ...
        images_per_figure, 'EEGlocs', EEGlocs, 'topo_cmap', topo_cmap);

    set(handles, 'pointer', 'arrow');

catch
    set(handles, 'pointer', 'arrow');
    icatb_displayErrorMsg;
end


function displayMethodCallback(hObject, event_data, handles)
% Display method callback


handles_data = get(handles, 'userdata');
dispParameters = handles_data.dispParameters;
value = get(hObject, 'value');
options = get(hObject, 'string');

% Display method
dispParameters.displayMethod = deblank(lower(options(value, :)));

numComp = dispParameters.numOfComp;

% Component numbers
compNumberH = findobj(handles, 'tag', 'component_number');

% Viewing set
viewingSetH = findobj(handles, 'tag', 'viewing_set');

selComp = get(compNumberH, 'value');

if strcmpi(dispParameters.displayMethod, 'component')
    set(compNumberH, 'min', 0, 'max', 2);
    set(compNumberH, 'value', [1:numComp]);
    set(viewingSetH, 'enable', 'on');
else
    set(compNumberH, 'value', selComp(1));
    set(compNumberH, 'min', 0, 'max', 1);
    set(viewingSetH, 'enable', 'off');
end

handles_data.dispParameters = dispParameters;

set(handles, 'userdata', handles_data);


function [sesInfo] = selectLocsFile(paramFile)
% Select EEG location file

load(paramFile);

locsFile = icatb_selectEntry('title', 'Select EEG channel location file for topoplot', 'typeEntity', 'file', 'filter', '*.locs', ...
    'typeSelection', 'single');

if isempty(locsFile)
    error('EEG channel location file must be selected for topoplot.');
end

% Update fields
sesInfo.userInput.locsFile = locsFile;
sesInfo.locsFile = locsFile;

drawnow;

icatb_save(paramFile, 'sesInfo');


function locsFileCallback(hObject, event_data, handles)
% Select EEG locs file

handles_data = get(handles, 'userdata');
dispParameters = handles_data.dispParameters;

paramFile = dispParameters.paramFile;

selectLocsFile(paramFile);

drawnow;

function timeScaleCallback(hObject, event_data, handles)
% Select time scale file

handles_data = get(handles, 'userdata');
dispParameters = handles_data.dispParameters;

paramFile = dispParameters.paramFile;

load(paramFile);

timeScaleFile = icatb_selectEntry('title', 'Select time scale for the EEG signal', 'typeEntity', 'file', 'filter', '*.set;*.asc;*.dat', ...
    'typeSelection', 'single');

drawnow;

if isempty(timeScaleFile)
    error('Time scale file for plotting EEG signal is not selected.');
end

sesInfo.userInput.timeScaleFile = timeScaleFile;

sesInfo.timeScaleFile = timeScaleFile;

[pp, fName, extn] = fileparts(timeScaleFile);

disp(['Loading time scale file ', timeScaleFile, ' ...']);
set(handles, 'pointer', 'watch');
try
    EEG = icatb_eeg_pop_loadset(timeScaleFile); % load .set/.dat file
    EEGTimeAxis = EEG.times;
catch
    data = load(timeScaleFile, '-ascii');
    if size(data, 1) == 1
        data = data';
    end
    EEGTimeAxis = data(:, 1);
end

set(handles, 'pointer', 'arrow');
EEGTimeAxis = EEGTimeAxis(:);

sesInfo.userInput.EEGTimeAxis = EEGTimeAxis;
sesInfo.EEGTimeAxis = EEGTimeAxis;

icatb_save(paramFile, 'sesInfo');


function [subjectICAFiles, calibratedMATFiles] = getSubjectMatFileNames(sesInfo)

subjectICAFiles = icatb_parseOutputFiles('icaOutputFiles', sesInfo.icaOutputFiles, 'numOfSub', sesInfo.numOfSub, 'numOfSess', sesInfo.numOfSess, 'flagTimePoints', sesInfo.flagTimePoints);

subFiles = cell(sesInfo.numOfSub*sesInfo.numOfSess, 1);
calibratedMATFiles = subFiles;
count = 0;
for nSub = 1:sesInfo.numOfSub
    for nSess = 1:sesInfo.numOfSess
        count = count + 1;
        subFiles{count} = deblank(subjectICAFiles(nSub).ses(nSess).name(1, :));
        calibratedMATFiles{count} = [sesInfo.calibrate_components_mat_file, num2str(nSub), '-', num2str(nSess), '.mat'];
    end
end

clear subjectICAFiles;
subjectICAFiles = subFiles;

function status = chkFile(file)

[pathstr, fN, extn] = fileparts(deblank(file));

fN2 = regexprep(fN, '(_\d.*)$', '_');

status = 0;
if (exist(fullfile(pathstr, [fN, extn]), 'file'))
    status = 1;
end

function outputFiles = chkOutputFiles(outputFiles, outDir, subjectICAFiles, calibratedMATFiles)
%% Truncate output files
%

setsToInclude = [];
for jj = 1:length(outputFiles)
    file = deblank(outputFiles(jj).ses(1).name(1, :));
    chkIndex = strmatch(file, subjectICAFiles, 'exact');
    compName = file;
    if (~isempty(chkIndex))
        compName = calibratedMATFiles{chkIndex};
    end
    compName = fullfile(outDir, compName);
    status = chkFile(compName);
    if (~status)
        continue;
    end
    setsToInclude = [setsToInclude, jj];
end

outputFiles = outputFiles(setsToInclude);

