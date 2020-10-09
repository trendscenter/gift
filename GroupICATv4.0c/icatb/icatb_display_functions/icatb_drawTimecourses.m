function icatb_drawTimecourses(parameters)
% Expanded view of the time course
% uses function callbacks and menus to display the options and the
% utilities

% load defaults
icatb_defaults;

% Screen Color Defaults
global BG_COLOR;
global BG2_COLOR;
global BUTTON_COLOR;
global BUTTON_FONT_COLOR;
global FG_COLOR;
global AXES_COLOR;
global FONT_COLOR;
global UI_FONTNAME;
global UI_FONTUNITS;
global UI_FS;

% set figure data the parameters
figureData = parameters;
outputDir = parameters.filesOutputDir;
clear parameters;

% open time course figure
timecourse_handle = icatb_getGraphics(figureData.titleFig, 'timecourse', 'expaned_view_timecourse');

% get all the options to include in options menu
input_parameters = icatb_displayOptions;

if (strcmpi(icatb_get_modality, 'fmri'))
    
    % Plot Options menu
    optionsMenu = uimenu(timecourse_handle, 'label', 'Options');
    timecourseOptionsMenu = uimenu(optionsMenu, 'label', 'Timecourse Options', 'callback', ...
        {@optionsWindowCallback, timecourse_handle}, 'userdata', input_parameters);
    clear input_parameters;
    
    figureData.optionsMenu = timecourseOptionsMenu;
    
    % Plot Utilities menu
    utilitiesMenu = uimenu(timecourse_handle, 'label', 'Utilities');
    % FFT
    fftMenu = uimenu(utilitiesMenu, 'label', 'Power Spectrum', 'callback', {@fftCallback, timecourse_handle});
    figureData.fftMenu = fftMenu;
    
    
    % show split time course menu
    if figureData.num_DataSets > 1
        splitTimecourseMenu = uimenu(utilitiesMenu, 'label', 'Split Timecourses', 'callback', ...
            {@splitTimecoursesCallback, timecourse_handle}); % SPLIT TIME COURSES
        figureData.splitTimecourseMenu = splitTimecourseMenu;
    end
    
    % Plot Help menu
    helpMenu = uimenu(timecourse_handle, 'label', 'GIFT-Help');
    % Time courses help
    timecoursesHelpMenu = uimenu(helpMenu, 'label', 'Expanded View of Timecourse', 'callback', ...
        'icatb_openHTMLHelpFile(''icatb_expanded_view_tc.htm'');');
    
    % Initialise vars
    sortingCriteria = [];
    sortingType = [];
    
    % get the sorting criteria and sorting type
    if isfield(figureData, 'sortParameters')
        sortingCriteria = lower(figureData.sortParameters.sortingCriteria);
        sortingType = lower(figureData.sortParameters.sortingType);
    end
    
    if strcmp(sortingCriteria, 'multiple regression') & ~strcmp(sortingType, 'spatial')
        % adjust ICA only for multiple regression
        adjustICAMenu = uimenu(optionsMenu, 'label', 'Adjust ICA', 'callback', ...
            {@adjustICACallback, timecourse_handle}); % ADJUST ICA
        % Regular method to compute event average
        eventAvgMenu = uimenu(utilitiesMenu, 'label', 'Event Average (Regular Method)', 'callback', ...
            {@eventRelatedCallback, timecourse_handle}, 'userdata', figureData.sortParameters.refInfo, 'tag', 'event_avg_regular'); % EVENT AVG
        % Deconvolution method to compute event average
        deconv_eventAvgMenu = uimenu(utilitiesMenu, 'label', 'Event Average (Deconvolution Method)', 'callback', ...
            {@eventRelatedCallback, timecourse_handle}, 'userdata', figureData.sortParameters.refInfo, 'tag', 'event_avg_deconv'); % Deconv method for EVENT AVG
        figureData.adjustICAMenu = adjustICAMenu;
        figureData.eventAvgMenu = eventAvgMenu;
        figureData.deconv_eventAvgMenu = deconv_eventAvgMenu;
    elseif strcmp(sortingCriteria, 'correlation') & ~strcmp(sortingType, 'spatial')
        % Regular method to compute event average
        eventAvgMenu = uimenu(utilitiesMenu, 'label', 'Event Average (Regular Method)', 'callback', ...
            {@eventRelatedCallback, timecourse_handle}, 'userdata', figureData.sortParameters.refInfo, 'tag', 'event_avg_regular'); % EVENT AVG
        % Deconvolution method to compute event average
        deconv_eventAvgMenu = uimenu(utilitiesMenu, 'label', 'Event Average (Deconvolution Method)', 'callback', ...
            {@eventRelatedCallback, timecourse_handle}, 'userdata', figureData.sortParameters.refInfo, 'tag', 'event_avg_deconv'); % Deconv method for EVENT AVG
        figureData.eventAvgMenu = eventAvgMenu;
        figureData.deconv_eventAvgMenu = deconv_eventAvgMenu;
    end
    
    
    
    % % Save figure data menu
    % saveTimecourseDataMenu = uimenu(optionsMenu, 'label', 'Save Figure Parameters', 'callback', ...
    %     {@saveTimecourseDataCallback, timecourse_handle}, 'userdata', outputDir);
    %
    % % Save figure data menu
    % loadTimecourseDataMenu = uimenu(optionsMenu, 'label', 'Load Figure Parameters', 'callback', ...
    %     {@loadTimecourseDataCallback, timecourse_handle}, 'userdata', outputDir);
    
    % Procedure: Initial plot will be an exact replica of the draw time courses
    % figure but plots in expanded view. Options can be selected using ICA
    % Options window. The results can be viewed using plot button.
    
    % draw axis
    axesPosition = [0.05 .1 .8 .8];
    axesTag = 'plot';
    % draw plot button
    plotWidth = 0.1; plotHeight = 0.1;
    pos = [axesPosition(1) + axesPosition(3) + 0.02 0.5 - 0.5*plotHeight plotWidth plotHeight];
    plotH = icatb_uicontrol('parent', timecourse_handle, 'units', 'normalized', 'style', 'pushbutton', 'string', ...
        'Reset', 'position', pos, 'callback', {@plotTC_expandedMode, timecourse_handle}, 'tag', 'plot_button');
    figureData.axesPosition = axesPosition;
    figureData.axesTag = axesTag;
    
    set(timecourse_handle, 'userdata', figureData);
    
    plotTC_expandedMode(plotH, [], timecourse_handle); % execute the plot time course callback
    
    
else
    axesPosition = [0.05 .1 .8 .8];
    axisH = axes('units', 'normalized', 'position', axesPosition);
    plot(figureData.undetrendICA, 'm', 'parent', axisH);
    xlabel('Subjects', 'parent', axisH);
    ylabel('ICA Loadings', 'parent', axisH);
    title(figureData.titleFig, 'parent', axisH);
    axis(axisH, 'tight');
    set(axisH, 'YColor', FONT_COLOR, 'XColor', FONT_COLOR);
    
end

function plotTC_expandedMode(handleObj, event_data, handles)
% plot time courses in expanded mode

% get the options from the Time course options menu
% pass the handle visibility off
% get the results here.

% use the undetrended timecourse
% apply smoothing defaults
% remove the mean (detrend(.., 0))
% use the sorting criteria as selected in sort components GUI

% adjust ICA uses the modeling parameters to remove the nuisance parameters
% and plots in the same axis

% Utilities menu use the options from the plot of the timecourse:
% 1ike split timecourse, event average

set(handles, 'pointer', 'watch');

if ~exist('flag_calculate_stats', 'var')
    flag_calculate_stats = 'yes';
end

% set the fields from timecourse options to figure data
figureData = get(handles, 'userdata'); % get the figure data

%%%%%%%%%%%% Options from options window %%%%%%%%%%%%%%%%%%%%
% get the options from the Time course options menu
% pass the handle visibility off
% get the results here.

optionsWindowCallback(figureData.optionsMenu, [], handles, 'off'); % options callback
input_parameters = get(figureData.optionsMenu, 'userdata'); % input parameters
timeCourseOptions = input_parameters.results; % get the time course options
clear input_parameters;
%%%%%%%%%%%% End of Options from options window %%%%%%%%%%%%%%%%%%%%
timecourseFields = fieldnames(timeCourseOptions); % get all the fields

for ii = 1:length(timecourseFields)
    figureData = setfield(figureData, timecourseFields{ii}, getfield(timeCourseOptions, timecourseFields{ii}));
end

clear timecourseFields;
clear timeCourseOptions;

% get all the options from the options window
smoothPara = figureData.icaSmoothPara;
smoothValue = figureData.icaSmoothValue;
modelNormalize = figureData.modelNormalize;
%modelZeroMean = figureData.modelZeroMean;
modelLinearDetrend = figureData.modelLinearDetrend;
icaDetrend = figureData.icaDetrend;
icaFlip = figureData.icaFlip;
icaSmoothPara = figureData.icaSmoothPara;
icaSmoothValue = figureData.icaSmoothValue;
eventWindowSize = figureData.eventWindowSize;
eventInterpFactor = figureData.eventInterpFactor;

% data sets information
num_DataSets = figureData.num_DataSets;
numSubjects = figureData.numSubjects;
numSessions = figureData.numSessions;
diffTimePoints = [];
if isfield(figureData, 'sortParameters')
    % different time points variable
    if isfield(figureData.sortParameters, 'diffTimePoints')
        diffTimePoints = figureData.sortParameters.diffTimePoints;
    else
        diffTimePoints = figureData.diffTimePoints;
    end
end

% Depending on the title name
% Select the appropriate sorting criteria
title_fig = get(handles, 'name');
compareTitle = lower(title_fig);
sortingCriteria = [];
sortingType = [];

if isfield(figureData, 'sortParameters')
    sortingCriteria = lower(figureData.sortParameters.sortingCriteria);
    sortingType = lower(figureData.sortParameters.sortingType);
end

% load defaults
icatb_defaults;

% Screen Color Defaults
global BG_COLOR;
global BG2_COLOR;
global FG_COLOR;
global AXES_COLOR;
global FONT_COLOR;
global LEGENDCOLOR;

currentAxes = get(handles, 'currentaxes');

delete(currentAxes);

axisH = axes('parent', handles, 'position', figureData.axesPosition, 'tag', figureData.axesTag);

%%%%%%%%%%% Time course calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use original ICA Timecourse and detrend ICA by removing mean
undetrendICA = figureData.undetrendICA;

if isempty(sortingType) | strcmp(sortingType, 'spatial')
    % remove the mean using the specified information in the options window
    timecourse = icatb_detrend(undetrendICA, num_DataSets, diffTimePoints, icaDetrend);
else
    timecourse = icatb_detrend(undetrendICA, num_DataSets, diffTimePoints, 0); % remove the mean
end

%timecourse = figureData.icaTimecourse;

% apply smoothing values
timecourse = icatb_smoothTimecourse(timecourse, diffTimePoints, num_DataSets, 1, ...
    icaSmoothPara, icaSmoothValue);

% ica timecourse calculations
if icaFlip == 1
    timecourse = timecourse*-1;
end
%%%%%%%%%%% End for Time course calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%% MODEL TC Calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model Timecourse calculations
if ~isempty(figureData.modelTimecourse)
    % get the model time course
    modelTimecourse = figureData.modelTimecourse;
    
    if modelLinearDetrend == 1
        % remove the linear trend
        modelTimecourse = icatb_detrend(modelTimecourse, num_DataSets, diffTimePoints, 1);
    end
    
    %     if modelZeroMean == 1
    %         % remove the mean
    %        modelTimecourse = icatb_detrend(modelTimecourse, num_DataSets, diffTimePoints, 0);
    %     end
    % normalize the model time course
    if modelNormalize == 1
        modelTimecourse = icatb_normalizeTimecourse(modelTimecourse, timecourse, diffTimePoints, num_DataSets);
    end
end

%%%%%%%%%%%%%%%%%% End for Model Time course Calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot appropriate sorting criteria title
if strcmp(sortingCriteria, 'correlation') & ~strcmp(sortingType, 'spatial')
    [sortingVal, timecourse] = icatb_correlateFunctions(modelTimecourse, timecourse, num_DataSets, diffTimePoints, icaDetrend);
    titlestr = ['Correlation = ', num2str(sortingVal)];
elseif strcmp(sortingCriteria, 'kurtosis') & ~strcmp(sortingType, 'spatial')
    [sortingVal, timecourse] = icatb_kurtosis(timecourse, num_DataSets, diffTimePoints, icaDetrend);
    titlestr = ['Kurtosis = ', num2str(sortingVal)];
    % Plot the linefit to model timecourse and the ICA timecourse
elseif strcmp(sortingCriteria, 'multiple regression') & ~strcmp(sortingType, 'spatial')
    [sortingVal, bRegress, ModelIndices, otherIndices, linearRegress, removeTrend, timecourse] = ...
        icatb_multipleRegression(modelTimecourse, timecourse, size(modelTimecourse, 2), num_DataSets, ...
        diffTimePoints, icaDetrend);
    % store the data for adjusting timecourses
    adjustData = struct('bRegress', bRegress, 'ModelIndices', ModelIndices, 'otherIndices', otherIndices, ...
        'linearRegress', linearRegress, 'removeTrend', removeTrend, 'icaLine', 'icaline');
    
    set(figureData.adjustICAMenu, 'userdata', adjustData);
    
    titlestr = ['R-Square Statistic = ', num2str(sortingVal)];
else
    titlestr = title_fig;
end

titlestr = [titlestr, ' (Temporal Detrend value = '  num2str(icaDetrend), ')'];

figureData.timecourse = timecourse; % set the field time course to figure data
if num_DataSets > 1
    parameters.numSubjects = numSubjects; % number of subjects
    parameters.numSessions = numSessions; % number of sessions
    parameters.titleFig = 'Time courses window';
    %parameters.timecourse = timecourse; % time course
    if exist('modelTimecourse', 'var')
        parameters.modelTimecourse = modelTimecourse; % model time course
    else
        parameters.modelTimecourse = [];
    end
    set(figureData.splitTimecourseMenu, 'userdata', parameters); % split time courses menu
end

storeHandlesCount = 0; legendString = {};

% check the model time courses
if ~isempty(figureData.modelTimecourse)
    referenceString = 'Reference Function';
    % plot the reference functions and append the legend string
    % Plot all the model timecourses
    [nrows, ncols] = size(modelTimecourse);
    for i = 1:ncols
        modelPlotHandle(i) = plot(modelTimecourse(:, i), 'y:', 'parent', axisH);
        legendString{length(legendString) + 1} = referenceString;
        hold on;
    end
    % store handles to be used later for plotting legend
    storeHandlesCount = storeHandlesCount + 1;
    %legendString{length(legendString) + 1} = referenceString; % store the reference string
    storeHandles(storeHandlesCount) = modelPlotHandle(end); % store only the end handle
    hold on;
    if strcmp(sortingCriteria, 'multiple regression') & ~strcmp(sortingType, 'spatial')
        % Plot ICA Timecourses
        icaTimecourseLine = plot(timecourse, 'm', 'tag', adjustData.icaLine, 'parent', axisH);
        hold on;
        storeHandlesCount = storeHandlesCount + 1;
        storeHandles(storeHandlesCount) = icaTimecourseLine;
        legendString{length(legendString) + 1} = 'ICA Timecourse';
        % % plot regressor
        lineRegressH = plot(linearRegress, '-.c', 'parent', axisH);
        hold off;
        storeHandlesCount = storeHandlesCount + 1;
        storeHandles(storeHandlesCount) = lineRegressH;
        legendString{length(legendString) + 1} = 'Linefit to ICA';
    else
        % Plot ICA Time course
        icaTimecourseLine = plot(timecourse, 'm', 'parent', axisH);
        storeHandlesCount = storeHandlesCount + 1;
        storeHandles(storeHandlesCount) = icaTimecourseLine;
        legendString{length(legendString) + 1} = 'ICA Timecourse';
        hold off;
    end
else
    icaTimecourseLine = plot(timecourse, 'm', 'parent', axisH);
    storeHandlesCount = storeHandlesCount + 1;
    storeHandles(storeHandlesCount) = icaTimecourseLine;
    legendString{length(legendString) + 1} = 'ICA Timecourse';
    hold off;
end

icatb_legend(legendString{:});

% Plot a dotted line separating the subjects
% Change later to incorporate different timepoints for each subject
if num_DataSets > 1
    lengthTotalTimecourse = length(timecourse);
    %lengthIndividualTimecourse = lengthTotalTimecourse/num_DataSets;
    axisYlim = get(axisH, 'ylim');
    ySpacing = (axisYlim(2) - axisYlim(1))/4;
    plotYY = (axisYlim(1):ySpacing:axisYlim(2));
    for numdataSets = 2:num_DataSets
        hold on;
        %lengthIndividualTimecourse = sum(figureData.refInfo.diffTimePoints(1:numdataSets - 1));
        %valX = lengthIndividualTimecourse*(numdataSets - 1);
        valX =  sum(figureData.sortParameters.diffTimePoints(1:numdataSets - 1));
        plotXX(1:length(plotYY)) = valX;
        plot(plotXX, plotYY, '--g', 'parent', axisH);
    end
end

axis(axisH, 'tight');

set(axisH, 'YColor', FONT_COLOR, 'XColor', FONT_COLOR);

set(handles, 'userdata', figureData);

clear figureData;

hold off;

title(titlestr, 'parent', axisH);
set(handles, 'pointer', 'arrow');

%%%%%%%%%%%%%%%%%% Function callbacks %%%%%%%%%%%%%%%%
function fftCallback(hObject, event_data, handles)
% powerSpectrum

psdFile = which('psd.m');

if ~isempty(psdFile)
    
    figureData = get(handles, 'userdata');
    
    icatb_defaults;
    % Screen Color Defaults
    global BG_COLOR;
    global BG2_COLOR;
    global FG_COLOR;
    global AXES_COLOR;
    global FONT_COLOR;
    global SMOOTHINGVALUE;
    global  SMOOTHPARA;
    
    fftHandle = icatb_getGraphics('Power Spectrum', 'normal', 'Power Spectrum Figure');
    timecourse = figureData.timecourse;
    fftH = axes('parent', fftHandle, 'units', 'normalized', 'position', [0.1 0.1 0.8 0.8]);
    % calculate power spectrum
    psd(timecourse, 500);
    axis(fftH, 'tight');
    set(fftH, 'YColor', FONT_COLOR, 'XColor', FONT_COLOR);
    % set the color to magenta
    set(get(gca, 'children'), 'color', [1 0 1]);
    clear figureData;
else
    disp('Need signal processing toolbox to view power spectrum');
    return;
end
% end for checking psd

% Splits time courses in the main window if there are more than one time
% course
function splitTimecoursesCallback(handleObj, evd, handles)

% get the user data associated with the push button
parameters = get(handleObj, 'userdata');
handles_data = get(handles, 'userdata');
% Get the variables here
numSubjects = parameters.numSubjects;
numSessions = parameters.numSessions;
titleFig = parameters.titleFig;
%timecourse = parameters.timecourse;
timecourse = handles_data.timecourse;
modelTimecourse = parameters.modelTimecourse;
diffTimePoints = [];
% get different time points
if isfield(handles_data, 'sortParameters')
    diffTimePoints = handles_data.sortParameters.diffTimePoints;
end
% convert all vars to structure format
[timeCourseStruct, modelTimeCourseStruct, xAxis] = icatb_convertTo_struct('timecourse', timecourse, 'model', ...
    modelTimecourse, 'numSubjects', numSubjects, 'numSessions', numSessions, 'difftimepoints', diffTimePoints);

% show multiple plots in a figure
icatb_showTimeCourses('TC', timeCourseStruct, 'x_axis', xAxis, 'numsubjects', numSubjects, ...
    'numsessions', numSessions, 'titlefig', 'Splitting Time courses', 'meandisplay', 'yes', 'SEM', ...
    'yes', 'model', modelTimeCourseStruct, 'timecourse_color', 'm', 'model_timecourse_color', 'y:');

function eventRelatedCallback(hObject, evd, handles)

% Computes the event related averages. Gets the ICA time course under
% consideration. If the units is seconds then selected onsets will be
% divided by TR.
% get the event related vector (correct Onset???). Interpolate the time course by a factor
% which is factor to multiply * TR.

% Get the related sequence where the events occur

% get the user data
figureData = get(handles, 'userdata');

eventAvgObjTag = get(hObject, 'tag');

if strcmpi(eventAvgObjTag, 'event_avg_regular')
    figTitle = 'Event Average (Regular Method)';
    fprintf('\n');
    disp('Using regular method to compute event average ...');
    fprintf('\n');
else
    figTitle = 'Event Average (Deconvolution Method)';
    fprintf('\n');
    disp('Using deconvolution method to compute event average ...');
    fprintf('\n');
end

% get the options from the options menu handle
optionsWindowCallback(figureData.optionsMenu, [], handles, 'off'); % options callback
input_parameters = get(figureData.optionsMenu, 'userdata'); % input parameters
timeCourseOptions = input_parameters.results; % get the time course options
clear input_parameters;
%%%%%%%%%%%% End of Options from options window %%%%%%%%%%%%%%%%%%%%

timecourseFields = fieldnames(timeCourseOptions); % get all the fields
% set the figure data with the corresponding fields
for ii = 1:length(timecourseFields)
    figureData = setfield(figureData, timecourseFields{ii}, getfield(timeCourseOptions, timecourseFields{ii}));
end

icatb_defaults;

global FONT_COLOR;

% Window size
windowSize = figureData.eventWindowSize; % in seconds
interpFactor = figureData.eventInterpFactor; % interpolation factor

eventData = struct('eventInterpFactor', interpFactor, 'eventWindowSize', windowSize);

icaTimecourse = figureData.timecourse; % time course

% Get the information about the data sets here
% number of subjects and sessions
num_DataSets = figureData.num_DataSets;
numSubjects = figureData.numSubjects;
numSessions = figureData.numSessions;
numOfSub = figureData.numOfSub;
numOfSess = figureData.numOfSess;
actualDataSets =  numOfSub*numOfSess;
detrendNumber = figureData.icaDetrend;

% get the subject string to be displayed as a string in listbox
if numSubjects > 1 & numSessions > 1
    count = 0;
    for ii = 1:numSubjects
        for jj = 1:numSessions
            count = count + 1;
            subject_string(count).name = ['Subject ', num2str(ii), ' Session ', num2str(jj)];
        end
    end
elseif numSubjects > 1 & numSessions == 1
    for ii = 1:numSubjects
        subject_string(ii).name = ['Subject ', num2str(ii)];
    end
elseif numSubjects == 1 & numSessions > 1
    for ii = 1:numSessions
        subject_string(ii).name = ['Session ', num2str(ii)];
    end
else
    subject_string.name = ['Subject''s Session'];
end

diffTimePoints = [];
% get different time points
if isfield(figureData, 'sortParameters')
    diffTimePoints = figureData.sortParameters.diffTimePoints;
    if isfield(figureData.sortParameters, 'session_number')
        data_sessionNumber = figureData.sortParameters.session_number;
    end
end

% check the spm mat flag
if isfield(figureData, 'spmMatFlag')
    spmMatFlag = figureData.spmMatFlag;
else
    spmMatFlag = 'not_specified';
end

clear figureData;
% convert timecourses to required structure format: incorporate different
% time points
[timeCourses] = icatb_convertTo_struct('timecourse', icaTimecourse, 'numSubjects', numSubjects, 'numSessions', ...
    numSessions, 'diffTimePoints', diffTimePoints);

clear icaTimecourse;

% get the figure data
refInfo = get(hObject, 'userdata');

%%%%%%%%%%%%%%%%% REGRESSOR GUI %%%%%%%%%%%%%%%%%%%
[selectedRegressors] = icatb_returnRegressor(refInfo, num_DataSets, numSubjects, ...
    numSessions, spmMatFlag, subject_string);

countDataSet = 0;
% Loop over number of subjects and sessions
for numSub = 1:numSubjects
    % check if the size of the SPM file is greater than one
    if size(refInfo.SPMFile, 2) > 1
        spmFileNumber = (numSub - 1)*numSessions + 1;
        % get the current spm file
        SPMfileName = refInfo.SPMFile(spmFileNumber).name;
        spmData = refInfo.spmAppData(spmFileNumber).data;
    else
        % get the current spm file
        SPMfileName = refInfo.SPMFile.name;
        % check the special case for session specific regressors
        if length(refInfo.spmAppData) > 1
            spmData = refInfo.spmAppData(1).data;
        else
            spmData = refInfo.spmAppData.data;
        end
    end
    
    % SPM information
    spmVar = spmData.SPM;
    
    flag_refFunc = spmData.flag_refFunc;
    
    clear spmData;
    
    % get the SPM mat file
    TimeResponse = spmVar.xY.RT;
    
    % Original Interpolation factor
    origInterpFactor = eventData.eventInterpFactor;
    
    % Interpolation factor
    interpFactor = ceil(origInterpFactor*TimeResponse);
    
    % loop over the sessions
    for numSess = 1:numSessions
        countDataSet = countDataSet + 1;
        % pull the session number from the selected regressors
        currentRegressor = deblank(selectedRegressors(countDataSet, :)); % get the current regressor
        % get the onsets and data session number using the regressor
        % name
        
        if strcmp(lower(flag_refFunc), 'session_specific')
            if numSessions == numOfSess
                data_sessionNumber = numSess;
            end
            [temp_data_sessionNumber, onset_number] = icatb_get_onsets(currentRegressor, spmVar);
        else
            [data_sessionNumber, onset_number] = icatb_get_onsets(currentRegressor, spmVar);
        end
        
        if numSubjects*numSessions == 1
            % displaying the reference function
            disp(['Selected reference function to compute event related average for subject''s session is ', currentRegressor]);
        elseif numSubjects == 1 & numSessions > 1
            % displaying the reference function
            disp(['Selected reference function to compute event related average for session  ', num2str(numSess), ' is ', currentRegressor]);
        elseif numSessions == 1 & numSubjects > 1
            % displaying the reference function
            disp(['Selected reference function to compute event related average for subject  ', num2str(numSub), ' is ', currentRegressor]);
        else
            % displaying the reference function
            disp(['Selected reference function to compute event related average for subject ', num2str(numSub), ...
                ' session ', num2str(numSess), ' is ', currentRegressor]);
        end
        
        % % determine the selected onset
        % % or current onset number
        if data_sessionNumber > length(spmVar.Sess)
            data_sessionNumber = length(spmVar.Sess);
        end
        
        %%%%%% Modifying selected onset
        if strcmp(lower(spmVar.xBF.UNITS), 'scans')
            selectedOnset = round(spmVar.Sess(data_sessionNumber).U(onset_number).ons);
        else
            selectedOnset = round((spmVar.Sess(data_sessionNumber).U(onset_number).ons) / TimeResponse);
        end
        
        timeCourses.sub(numSub).sess(numSess).tc = icatb_detrend(timeCourses.sub(numSub).sess(numSess).tc, 1, ...
            length(timeCourses.sub(numSub).sess(numSess).tc), detrendNumber);
        
        if strcmpi(eventAvgObjTag, 'event_avg_regular')
            % Calculate event average
            eventAvg = icatb_calculate_eventAvg(timeCourses.sub(numSub).sess(numSess).tc, origInterpFactor, TimeResponse, eventData.eventWindowSize, selectedOnset);
            %timeAxis = (1:length(eventAvg)) / origInterpFactor; % time axis
            timeAxis = icatb_interp((1:eventData.eventWindowSize), origInterpFactor);
        else
            eventAvg = icatb_calculate_eventAvg_deconv(timeCourses.sub(numSub).sess(numSess).tc, TimeResponse, eventData.eventWindowSize, selectedOnset);
            timeAxis = (1:length(eventAvg));
        end
        
        timeCourseStruct.sub(numSub).sess(numSess).tc = eventAvg; % event average
        xAxis.sub(numSub).sess(numSess).tc = timeAxis; % time axis
    end
    if exist('spmVar', 'var')
        clear spmVar;
    end
end
% end for loop over number of subjects and sessions


% number of data sets is more than 1
if num_DataSets == 1
    % if the number of data sets is 1 else show the plots in a figure
    % browse
    
    plotH = icatb_getGraphics(figTitle, 'normal', 'figure');
    
    plot(timeAxis, eventAvg, 'm');
    
    title(figTitle, 'color', FONT_COLOR, 'HorizontalAlignment', 'center');
    
    xlabel(['Time in seconds where ', ' window size = ', num2str( eventData.eventWindowSize)]);
    
    set(gca, 'YColor', FONT_COLOR, 'XColor', FONT_COLOR, 'position', [0.1 0.1 0.8 0.8]);
    
    axis(gca, 'tight');
    
else
    icatb_showTimeCourses('TC', timeCourseStruct, 'x_axis', xAxis, 'numsubjects', numSubjects, 'numsessions', numSessions, ...
        'titlefig', figTitle, 'meandisplay', 'yes', 'SEM', 'yes', 'timecourse_color', 'm');
end

% context callback
function adjustICACallback(hObject, evd, handles)
% context menu callback
% open figure window prompting the user to select the reference function
% remove the nuisance parameters and the other model coefficients from the
% ica timecourse
% plot the new ica timecourse

% get the userdata
adjustData = get(hObject, 'userdata');

% store the data for adjusting timecourses
bRegress = adjustData.bRegress; ModelIndices = adjustData.ModelIndices;
otherIndices = adjustData.otherIndices; linearRegress = adjustData.linearRegress;
removeTrend = adjustData.removeTrend;
icaLine = adjustData.icaLine;

% figure data
figData = get(handles, 'userdata');

timecourse = figData.timecourse;
detrendNumber = figData.icaDetrend;

% check the spm mat flag
if isfield(figData, 'spmMatFlag')
    spmMatFlag = figData.spmMatFlag;
else
    spmMatFlag = 'not_specified';
end

% data information
numSubjects = figData.numSubjects;
numSessions = figData.numSessions;
num_DataSets = figData.num_DataSets;
modelTimecourse = figData.modelTimecourse;

% reference function information
if ~isfield(figData, 'sortParameters')
    error('sortParameters field is missing in figData');
end
diffTimePoints = figData.sortParameters.diffTimePoints;
refInfo = figData.sortParameters.refInfo;
num_Regress = size(refInfo.modelIndex, 2);

% spm file names
SPMFile = refInfo.SPMFile;

if length(SPMFile) > 1
    % spmNames
    spmNames = str2mat(SPMFile.name);
else
    spmNames = SPMFile.name;
end

% % details of the help push button
% helpDetails.title = 'Adjust ICA time course';
%
% D(1).string = ['For adjusting ICA time-course, only one reference function is used. Please click on the left top listbox to select the data set and the reference functions for that data set will be displayed in the left bottom listbox.'];
% D(size(D,2)+1).string = '';
% D(size(D,2)+1).string = ['Select the reference function in the left bottom listbox and the selected string will appear in the right top listbox. These selections can be made as long as Ok button is pressed.'];
% str = {D.string};

% string for help button
%helpDetails.str = str;

% get the subject string to be displayed as a string in listbox
if numSubjects > 1 & numSessions > 1
    count = 0;
    for ii = 1:numSubjects
        for jj = 1:numSessions
            count = count + 1;
            subject_string(count).name = ['Subject ', num2str(ii), ' Session ', num2str(jj)];
        end
    end
elseif numSubjects > 1 & numSessions == 1
    for ii = 1:numSubjects
        subject_string(ii).name = ['Subject ', num2str(ii)];
    end
elseif numSubjects == 1 & numSessions > 1
    for ii = 1:numSessions
        subject_string(ii).name = ['Session ', num2str(ii)];
    end
else
    subject_string.name = ['Subject''s Session'];
end

if (size(refInfo.modelIndex, 2) == 1)
    icatb_dialogBox('title', 'Adjust ICA', 'textBody', 'Cannot adjust ICA timecourse as only one regressor is selected during sorting.', ...
        'textType', 'large');
    return;
end

%%%%%%%%%%%%%%%%% REGRESSOR GUI %%%%%%%%%%%%%%%%%%%

[selectedRegressors, selectedIndices] = icatb_returnRegressor(refInfo, num_DataSets, numSubjects, ...
    numSessions, spmMatFlag, subject_string, 'adjust_ica');

% Remove the nuisance parameters and the other than selected model indices
rowStart = 1;
for nn = 1:num_DataSets
    rowEnd = sum(diffTimePoints(1:nn));
    % Regressor matrix
    [X] = icatb_modelX(modelTimecourse(rowStart:rowEnd, :), diffTimePoints(nn), detrendNumber);
    % Total terms
    TotalTerms = size(X, 2);
    regressIndices = [1:TotalTerms];
    % Not selected indices
    not_selected_indices = find(regressIndices ~= selectedIndices(nn));
    % Other indices
    otherIndxs = TotalTerms*(nn - 1) + (num_Regress + 1 : TotalTerms);
    % Model indices
    modelIndxs = TotalTerms*(nn - 1) + (1:num_Regress);
    allIndxs = [modelIndxs, otherIndxs];
    notSelIndices = allIndxs(not_selected_indices);
    timecourse(rowStart:rowEnd) = timecourse(rowStart:rowEnd) - X(:, not_selected_indices)*bRegress(notSelIndices);
    rowStart = rowEnd + 1;
end

icaTimecourseLine = findobj(handles, 'tag', icaLine);

axisH = get(handles, 'currentaxes');

newStr = ['ICA time course (ICA tc) adjusted by the selected reference function. Click plot to view the old ICA tc.'];
title(newStr, 'parent', axisH);

% plot the adjusted time course
% Only plot the adjusted time course and don't calculate the ICA time
% course
set(icaTimecourseLine, 'Ydata', timecourse);

axis(axisH, 'tight');

figData.timecourse = timecourse;

set(handles, 'userdata', figData);

function optionsWindowCallback(hObject, event_data, handles, figVisibility)
% hObject - menu string
% handles - time course figure

if ~exist('figVisibility', 'var')
    figVisibility = 'on';
end
input_parameters = get(hObject, 'userdata'); % get the menu user data
defaults = input_parameters.defaults; % get the defaults
inputParameters = input_parameters.inputParameters; % get the input parameters
clear input_parameters;
% plot the figure window to display the options
[output_parameters, valueButton] = icatb_OptionsWindow(inputParameters, defaults, figVisibility);
input_parameters.inputParameters = output_parameters.inputParameters;
input_parameters.defaults = defaults;
input_parameters.results = output_parameters.results; % get the results
set(hObject, 'userdata', input_parameters); % set the object user data (updated options)
% Show the title after the options are updated by clicking on the ok button
if valueButton == 1
    plotH = findobj(handles, 'tag', 'plot_button');
    plotTC_expandedMode(plotH, [], handles); % execute the plot time course callback
    %     % set the title at the plot window
    %     axisH = get(handles, 'currentaxes');
    %     newStr = ['Click plot to view the time courses with the updated options.'];
    %     title(newStr, 'parent', axisH);
    
    % call the plot function here
    
end

function saveTimecourseDataCallback(hObject, event_data, handles)
% save figure data

outputDir = get(hObject, 'userdata');
if isempty(outputDir)
    outputDir = pwd;
end
figureData = get(handles, 'userdata');

% open input dialog box
prompt = {'Enter Valid File Name:'};
dlg_title = 'Save Figure Parameters as';
num_lines = 1;
def = {'TimecourseParameters'};
% save the file with the file name specified
fileName = icatb_inputdlg2(prompt, dlg_title, num_lines, def);

if ~isempty(fileName)
    % make a full file
    fileName = [fileName{1}, '.mat'];
    fileName = fullfile(outputDir, fileName);
    icatb_save(fileName, 'figureData');
    disp(['Figure parameters are saved in file: ', fileName]);
end

function loadTimecourseDataCallback(hObject, event_data, handles)
% load figure data

try
    outputDir = get(hObject, 'userdata');
    
    global FIGURE_DATA
    
    % open file selection box
    [P] = icatb_selectEntry('typeEntity', 'file', 'title', 'Select Timecourse Parameters File', 'filter', '*.mat', ...
        'startPath', outputDir);
    if ~isempty(P)
        load(P);
        if ~exist('figureData', 'var')
            error('Not a valid timecourse parameters file');
        end
        FIGURE_DATA = figureData;
        clear figureData;
        D(1).string = 'Global FIGURE_DATA variable is created.';
        D(size(D, 2) + 1).string = '';
        D(size(D, 2) + 1).string = 'Type the following statements at the command prompt:';
        D(size(D, 2) + 1).string = '';
        D(size(D, 2) + 1).string = 'global FIGURE_DATA; FIGURE_DATA';
        disp(str2mat(D.string));
    end
catch
    icatb_errorDialog(lasterr, 'Figure Parameters');
    icatb_displayErrorMsg;
end