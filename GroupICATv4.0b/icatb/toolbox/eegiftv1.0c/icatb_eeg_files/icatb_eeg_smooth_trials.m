function icatb_eeg_smooth_trials(erpFiles, baseLine, smoothTrials, moving_average_window, sorting_trials, outFileName, outputDir)
% Smooth trials using specified moving average and write out the smoothed
% images. This code is based on Tom Eichele's work. We use moving average
% code written by Carlos Adri�n Vargas Aguilera that was posted on Matlab file
% exchange.
%
% Inputs:
%
% 1. erpFiles - ERP files in .set format
% 2. baseLine - Vector of length 2 to correct baseline
% 3. smoothTrials - Option is provided to smooth trials using moving point
%                   average (moving_average_window).
% 4. moving_average_window - Type of moving average window.
% 5. sorting_trials - Options are 'No', 'Condition & Latency'.
% 6. outFileName - File name used when writing analyze or Nifti data.
% 7. outputDir - Output directory to save analyze or Nifti data.
%

icatb_defaults;
% EEG defaults
global EEG_RMBASE;

if ~exist('erpFiles', 'var')
    erpFiles = icatb_selectEntry('typeSelection', 'multiple', 'typeEntity', 'file', 'title', 'Select ERP files for subjects...', ...
        'filter', '*.set');
end

if isempty(erpFiles)
    error('ERP files are not selected for the analysis');
end


drawnow;

% When input parameters are not specified use a figure window to get the
% EEG parameter info.
if ~exist('baseLine', 'var')

    numParameters = 1;

    % define all the input parameters in a structure
    inputText(numParameters).promptString = 'Enter vector to correct baseline (ms)';
    inputText(numParameters).uiType = 'edit';
    inputText(numParameters).answerString = ['[', num2str(EEG_RMBASE), ']'];
    inputText(numParameters).dataType = 'numeric';
    inputText(numParameters).tag = 'baseLine';
    inputText(numParameters).enable = 'on';

    numParameters = numParameters + 1;

    % define all the input parameters in a structure
    inputText(numParameters).promptString = 'Do you want to smooth EEG trials?';
    inputText(numParameters).uiType = 'popup';
    inputText(numParameters).answerString = {'Yes', 'No'};
    inputText(numParameters).dataType = 'string';
    inputText(numParameters).tag = 'smoothTrials';
    inputText(numParameters).enable = 'on';

    numParameters = numParameters + 1;

    % define all the input parameters in a structure
    inputText(numParameters).promptString = 'Enter a number for moving average window';
    inputText(numParameters).uiType = 'edit';
    inputText(numParameters).answerString = '3';
    inputText(numParameters).dataType = 'numeric';
    inputText(numParameters).tag = 'moving_average_window';
    inputText(numParameters).enable = 'on';

    numParameters = numParameters + 1;

    % define all the input parameters in a structure
    inputText(numParameters).promptString = 'How do you want to sort trials?';
    inputText(numParameters).uiType = 'popup';
    inputText(numParameters).answerString = {'No', 'Condition & Latency'};
    inputText(numParameters).dataType = 'string';
    inputText(numParameters).tag = 'sorting_trials';
    inputText(numParameters).enable = 'on';


    numParameters = numParameters + 1;

    % define all the input parameters in a structure
    inputText(numParameters).promptString = 'Enter file name (not full path) for writing data';
    inputText(numParameters).uiType = 'edit';
    inputText(numParameters).answerString = 'erp';
    inputText(numParameters).dataType = 'string';
    inputText(numParameters).tag = 'outFileName';
    inputText(numParameters).enable = 'on';


    figureTag = 'eeg_parameters';

    % EEG parameters GUI
    inputHandle = icatb_getGraphics('EEG Parameters', 'normal', figureTag, 'on');

    % set figure data
    set(inputHandle, 'menubar', 'none');

    % Help Menu
    helpMenu = uimenu('parent', inputHandle, 'label', 'EEGIFT-Help');
    htmlHelpMenu = uimenu(helpMenu, 'label', 'Import Data', 'callback', 'icatb_openHTMLHelpFile(''icatb_eeg_import_data.htm'');');

    % plot the input parameters to the parameter figure
    %icatb_plotInputPara('input_prefix', inputText, 'controls_to_plot', cellstr(str2mat(inputText.tag)), 'handles', inputHandle);
    inputHandle = icatb_plot_controls_fig(inputText, figureTag, 'on', 'Done', 'Cancel', cellstr(str2mat(inputText.tag)), inputHandle);

    handles_data.inputText = inputText;
    set(inputHandle, 'userdata', handles_data);

    % Specify object callbacks
    set(findobj(inputHandle, 'tag', 'smoothTrials'), 'callback', {@smoothTrialsCallback, inputHandle});

    set(findobj(inputHandle, 'tag', 'Cancel'), 'callback', 'delete(gcbf)');

    set(findobj(inputHandle, 'tag', 'Done'), 'callback', {@doneCallback, inputHandle});

    try
        set(inputHandle, 'visible', 'on');
        waitfor(inputHandle);
    catch
        if ishandle(inputHandle)
            delete(inputHandle);
        end
        %error('Figure window to enter eeg parameters was quit');
    end

    % Get input from user
    if isappdata(0, 'eeg_para_app_data')
        eeg_para_app_data = getappdata(0, 'eeg_para_app_data');
        rmappdata(0, 'eeg_para_app_data');
    else
        error('Need EEG input parameters for smoothing trials');
    end


    %%%% Generate vars using field names from structure %%%%%%
    fldNames = fieldnames(eeg_para_app_data);

    for nn = 1:length(fldNames)
        tempVar = getfield(eeg_para_app_data, fldNames{nn});
        if isnumeric(tempVar)
            tempVar = num2str(tempVar);
            strToEvaluate = [fldNames{nn}, ' = [', tempVar, '];'];
        else
            strToEvaluate = [fldNames{nn}, ' = ', '''', tempVar, '''', ';'];
        end
        eval(strToEvaluate);
    end
    %%%% End for generating vars from fieldnames %%%%%%

end

drawnow;

if ~exist('moving_average_window', 'var')
    moving_average_window = 3;
end

if ~exist('smoothTrials', 'var')
    smoothTrials = 'no';
end

if ~exist('sorting_trials', 'var')
    sorting_trials = 'no';
end

if ~exist('outFileName', 'var')
    outFileName = 'erp';
end

if (length(baseLine) ~= 2)
    error('Vector to correct baseline must be numeric and of length 2');
end

if strcmpi(smoothTrials, 'yes')
    if moving_average_window <= 0
        error('Enter a positive integer for moving average window');
    end
    moving_average_window = ceil(moving_average_window);
end

[status, message] = icatb_errorCheck(outFileName, 'output_prefix', 'outFileName');

if (status == 0)
    error(message);
end

if ~exist('outputDir', 'var')
    outputDir = icatb_selectEntry('typeSelection', 'single', 'typeEntity', 'directory', 'title', 'Select a directory to save data');
end


% Number of subjects
numSubjects = size(erpFiles, 1);


if isempty(outputDir)
    error('Output directory for writing the images is not selected');
end

%%% End for getting vars %%%%

cd(outputDir);


% save smooth trials parameters
EEGParameters = struct('erpFiles', erpFiles, 'baseLine', baseLine, 'smoothTrials', smoothTrials, 'moving_average_window', moving_average_window, 'sorting_trials', sorting_trials, ...
    'outFileName', outFileName);

eegParaFile = fullfile(outputDir, [outFileName, '_smooth_trials_parameters.mat']);
icatb_save(eegParaFile, 'EEGParameters');

disp(['Smoothing trials information is stored in ', eegParaFile]);

fprintf('\n');

drawnow;

% Loop over subjects
for nSub = 1:numSubjects

    fprintf('\n');

    currentErpFile = deblank(erpFiles(nSub, :));
    disp(['Subject ', currentErpFile]);

    [pathstr, fileN, extn] = fileparts(currentErpFile);

    % Use EEGLAB function to load set file
    EEG = icatb_eeg_pop_loadset(currentErpFile); % load .set/.dat file

    % Correct baseline
    EEG = icatb_eeg_pop_rmbase(EEG, baseLine);

    % No. of trials
    n_trials = size(EEG.data, 3);
    % No. of events
    events = (1:n_trials);

    if isfield(EEG, 'event')
        % Get events and latency
        if isfield(EEG.event, 'urevent') && (isfield(EEG, 'urevent') && ~isempty(EEG.urevent))
            % Update number of trials
            n_trials = length(EEG.urevent);
            % Update number of events
            events = [EEG.event.urevent];
            if ~isnumeric(EEG.urevent(1).type)
                eventType = cellstr(str2mat(EEG.urevent.type));
            else
                eventType = [EEG.urevent.type];
                eventType = eventType(:);
            end
            latency = [EEG.urevent.latency];
        else
            if isfield(EEG.event, 'type')
                if ~isnumeric(EEG.event(1).type)
                    eventType = cellstr(str2mat(EEG.event.type));
                else
                    eventType = [EEG.event.type];
                    eventType = eventType(:);
                end
                latency = [EEG.event.latency];
            end
        end
        % End for getting events and latency
    end

    % No. of channels and samples
    n_chans = size(EEG.data, 1);
    n_samples = size(EEG.data, 2);

    % Initialise data
    data = zeros(n_chans, n_samples, n_trials);

    % Check missing trials
    indices = ones(n_trials, 1);
    indices(events) = 0;
    indices = find(indices ~= 0);
    % End for checking missing trials

    meanEEGData = mean(double(EEG.data), 3);

    % Pad missing trials with grand mean
    for nM = 1:length(indices)
        data(:, :, indices(nM)) =  meanEEGData;
    end
    % End for padding missing trials


    clear meanEEGData;

    data(:, :, events) = double(EEG.data);

    clear EEG;

    %%%%%%% Sort Trials %%%%%%
    if exist('eventType', 'var') && exist('latency', 'var')
        data = sort_trials(data, latency, eventType, sorting_trials);
    end

    %%% End for sorting trials %%%%%

    if strcmpi(smoothTrials, 'yes')

        disp(['Smoothing trials using ', num2str(moving_average_window), ' point moving average ...']);

        %%%% Smooth across single trials with 3-point moving average %%%%
        for nChans = 1:n_chans
            disp(['Smoothing subject ', num2str(nSub), ' channel ', num2str(nChans), ' trials ...']);

            if nChans == 1
                statusH = icatb_statusBar('init', n_chans, 'Smoothing Trials...', '', '');
            end

            % Loop over samples
            for nSamp = 1:n_samples
                data(nChans, nSamp, :) = icatb_moving_average(data(nChans, nSamp, :), moving_average_window);
            end
            % End loop over samples

            statusH = icatb_statusBar('increment', 1);

        end
        % finish status bar
        icatb_statusBar('finished');
        %%% End for smoothing across single trials %%%%

    end

    % Create directory
    if ~exist(fullfile(outputDir, fileN), 'dir')
        mkdir(outputDir, fileN);
    end


    % Output file
    outFile = fullfile(outputDir, fileN, [outFileName, '.mat']);

    disp(['Saving channels in file ', outFile, ' ...']);

    % Data must be of dimensions time by trials by electrodes
    data = permute(data, [2 3 1]);

    icatb_save(outFile, 'data');

    clear data;

    disp('Done writing channels');

    fprintf('\n');

end
% End loop over subjects


disp('Done');


function doneCallback(hObject, event_data, handles)
% Done callback

handles_data = get(handles, 'userdata');
inputText = handles_data.inputText;

eeg_para_app_data = struct;
for nI = 1:length(inputText)
    currentObj = findobj(handles, 'tag', inputText(nI).tag);
    getStr = get(currentObj, 'string');
    getVal = get(currentObj, 'value');
    if strcmpi(inputText(nI).uiType, 'popup')
        currentStr = deblank(getStr{getVal});
    else
        currentStr = getStr;
    end

    if strcmpi(inputText(nI).dataType, 'numeric')
        try
            currentStr = str2num(currentStr);
        catch
            error('Error:CheckInput', 'Number entered for %s\n is not a valid number', inputText(nI).promptString);
        end
    end

    eeg_para_app_data = setfield(eeg_para_app_data, inputText(nI).tag, currentStr);
end

% Set application data
setappdata(0, 'eeg_para_app_data', eeg_para_app_data);

delete(handles);


function smoothTrialsCallback(hObject, event_data, handles)
% Avg trials callback

movingAverageH = findobj(handles, 'tag', 'moving_average_window');

getStr = get(hObject, 'string');
getVal = get(hObject, 'value');

if strcmpi(deblank(getStr{getVal}), 'yes')
    set(movingAverageH, 'enable', 'on');
else
    set(movingAverageH, 'enable', 'off');
end

function data = sort_trials(data, latency, eventType, sorting_trials)
% Sort trials callback


% Sort by condition and latency
if strcmpi(sorting_trials, 'condition & latency')

    disp(['Sorting trials by ', sorting_trials, ' ...']);
    % Get unique event types
    uniqueCond = unique(eventType);
    % Sort by conditions
    uniqueCond = sort(uniqueCond);
    index = repmat(struct('ind', []), 1, length(uniqueCond));
    % Deal numbers
    if isnumeric(uniqueCond)
        % Loop over conditions
        for n = 1:length(uniqueCond)
            inds = find(eventType == uniqueCond(n));
            % Sort by latency
            [temp, sortInd] = sort(latency(inds));
            index(n).ind = inds(sortInd);
            index(n).ind = index(n).ind(:)';
        end
        % End loop over conditions
    else
        % Treat eventTypes as strings
        % Loop over conditions
        for n = 1:length(uniqueCond)
            inds = strmatch(uniqueCond{n}, eventType, 'exact');
            % Sort by latency
            [temp, sortInd] = sort(latency(inds));
            index(n).ind = inds(sortInd);
            index(n).ind = index(n).ind(:)';
        end
        % End loop over conditions
    end

    indices = [index.ind];

    if length(indices) ~= size(data, 3)
        error('Error:SortTrials', 'No. of trials (%s) does not match the number of events (%s)\n', num2str(size(data, 3)), num2str(length(indices)));
    end
    % Sorted trials
    data = data(:, :, indices);

    disp(['Done sorting trials by ', sorting_trials]);

end
% End for sorting by condition and latency