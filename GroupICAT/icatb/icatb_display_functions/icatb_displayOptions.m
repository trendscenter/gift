% Expanded view of the time course Plot.
% using function callbacks to handle the plotting more efficiently
% USE MENUS to display the options and the utilities

% MENU 1:
% Options menu with list like Time course, Adjust ICA
% Time course Options string contains:
% a.) Model time course options
% Include: Normalize and Zero mean 
% b.) ICA time course options
% Include: Flip ICA time course, DETRENDNUMBER, SMOOTHING FACTOR
% c.) Event Related Average Options
% Include: Window Size, Interpolation Factor

% after applying these options adjust the time course

function input_parameters = icatb_displayOptions

icatb_defaults;

global DETRENDNUMBER;
% smoothing defaults
global SMOOTHINGVALUE;
global SMOOTHPARA;
% event average defaults
global EVENTAVG_WINDOW_SIZE;
global EVENTAVG_INTERP_FACTOR;

% List includes FFT, Split Time courses, Event Related Average, ...

%%%%%%%%% Input Parameters Structure
numParameters = 1;
inputParameters(numParameters).listString = 'Model Timecourse';
optionNumber = 1;
% Option 1 of parameter 1
options(optionNumber).answerString = 'Normalize Model Timecourse'; options(optionNumber).promptString = '';
options(optionNumber).uiType = 'checkbox'; options(optionNumber).value = 1;
options(optionNumber).tag = 'modelNormalize'; options(optionNumber).answerType = 'numeric';
options(optionNumber).flag = 'delete';
options(optionNumber).uiPos = [];

% optionNumber = optionNumber + 1;
% % Option 2 of parameter 1
% options(optionNumber).answerString = 'Zero Mean Model Timecourse'; options(optionNumber).promptString = '';
% options(optionNumber).uiType = 'checkbox'; options(optionNumber).value = 0;
% options(optionNumber).tag = 'modelZeroMean'; options(optionNumber).answerType = 'numeric';
% options(optionNumber).flag = 'delete';
% inputParameters(numParameters).options = options; % will be used in plotting the controls in a frame

optionNumber = optionNumber + 1;
% Option 2 of parameter 1
options(optionNumber).answerString = 'Linear Detrend Model Timecourse'; options(optionNumber).promptString = '';
options(optionNumber).uiType = 'checkbox'; options(optionNumber).value = 0;
options(optionNumber).tag = 'modelLinearDetrend'; options(optionNumber).answerType = 'numeric';
options(optionNumber).flag = 'delete';
options(optionNumber).uiPos = [];
inputParameters(numParameters).options = options; % will be used in plotting the controls in a frame
clear options;

numParameters = numParameters + 1; % second parameter
inputParameters(numParameters).listString = 'ICA Timecourse';
optionNumber = 1;
% Option 1 of parameter 2
options(optionNumber).answerString = str2mat('0', '1', '2', '3'); options(optionNumber).promptString = 'Detrend Number';
matchedIndex = strmatch(num2str(DETRENDNUMBER), options(optionNumber).answerString, 'exact');
options(optionNumber).uiType = 'popup'; options(optionNumber).value = matchedIndex;
options(optionNumber).tag = 'icaDetrend'; options(optionNumber).answerType = 'numeric';
options(optionNumber).flag = 'delete';
options(optionNumber).uiPos = [0.16 0.045];

optionNumber = optionNumber + 1;
% Option 2 of parameter 2
options(optionNumber).answerString = 'Flip ICA Timecourse'; options(optionNumber).promptString = '';
options(optionNumber).uiType = 'checkbox'; options(optionNumber).value = 0;
options(optionNumber).tag = 'icaFlip'; options(optionNumber).answerType = 'numeric';
options(optionNumber).flag = 'delete';
options(optionNumber).uiPos = [];

optionNumber = optionNumber + 1;
% Option 2 of parameter 2
options(optionNumber).answerString = str2mat('Yes', 'No');
matchedIndex = strmatch(lower(SMOOTHPARA), lower(options(optionNumber).answerString), 'exact');
options(optionNumber).promptString = 'Smoothing Parameter';
options(optionNumber).uiType = 'popup'; options(optionNumber).value = matchedIndex;
options(optionNumber).tag = 'icaSmoothPara'; options(optionNumber).answerType = 'string';
options(optionNumber).flag = 'delete';
options(optionNumber).uiPos = [0.16 0.045];

optionNumber = optionNumber + 1;
% Option 3 of parameter 2
options(optionNumber).answerString = num2str(SMOOTHINGVALUE); options(optionNumber).promptString = 'Smoothing Value';
options(optionNumber).uiType = 'edit'; options(optionNumber).value = 0;
options(optionNumber).tag = 'icaSmoothValue'; options(optionNumber).answerType = 'numeric';
options(optionNumber).flag = 'delete';
options(optionNumber).uiPos = [0.16 0.045];

inputParameters(numParameters).options = options; % will be used in plotting the controls in a frame
clear options;

numParameters = numParameters + 1; % Third parameter
inputParameters(numParameters).listString = 'Event Average';
optionNumber = 1;
% Option 1 of parameter 3
options(optionNumber).answerString = num2str(EVENTAVG_WINDOW_SIZE); options(optionNumber).promptString = 'Window Size in Seconds';
options(optionNumber).uiType = 'edit'; options(optionNumber).value = 0;
options(optionNumber).tag = 'eventWindowSize'; options(optionNumber).answerType = 'numeric';
options(optionNumber).flag = 'delete';
options(optionNumber).uiPos = [0.16 0.045];

optionNumber = optionNumber + 1;
% Option 2 of parameter 3
options(optionNumber).answerString = num2str(EVENTAVG_INTERP_FACTOR); options(optionNumber).promptString = 'Interpolation Factor';
options(optionNumber).uiType = 'edit'; options(optionNumber).value = 0;
options(optionNumber).tag = 'eventInterpFactor'; options(optionNumber).answerType = 'numeric';
options(optionNumber).flag = 'delete';
options(optionNumber).uiPos = [0.16 0.045];

inputParameters(numParameters).options = options; % will be used in plotting the controls in a frame
clear options;

% store the input parameters in two fields
input_parameters.inputParameters = inputParameters;
input_parameters.defaults = inputParameters;


