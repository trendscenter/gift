function icatb_eventAverage(paramFile, selSubjects, selSessions, compNumber, eventAvgMethod)
% Calculate event average
% Figure window will open to select subjects sessions, component and
% regressors

icatb_defaults;

% event average defaults
global EVENTAVG_WINDOW_SIZE;
global EVENTAVG_INTERP_FACTOR;
global DETRENDNUMBER;
global PARAMETER_INFO_MAT_FILE;
global METHOD_ENTERING_REGRESSORS;
global TXTFILE_REGRESSORS;
global SMOOTHINGVALUE;
global SMOOTHPARA;

% Get info from defaults
eventWindowSize = EVENTAVG_WINDOW_SIZE;
eventInterpFactor = EVENTAVG_INTERP_FACTOR;
detrendNumber = DETRENDNUMBER;
% type of entering regressors
analType = lower(METHOD_ENTERING_REGRESSORS);

% Load parameter file
filterStr = ['*', PARAMETER_INFO_MAT_FILE, '*.mat'];
if ~exist('paramFile', 'var')
    [paramFile] = icatb_selectEntry('typeEntity', 'file', 'title', 'Select Parameter File', 'filter', filterStr);
end
[pathstr, fileName] = fileparts(paramFile);
outputDir = pathstr;
cd(pathstr);
load(paramFile);

drawnow;

if ~exist('sesInfo', 'var')
    error('Please select the ICA parameter file');
end

if isfield(sesInfo, 'modality')
    if ~strcmpi(sesInfo.modality, 'fmri')
        error('You have selected EEG parameter file. Please select the parameter file using GIFT toolbox.');
    end
end

diffTimePoints = sesInfo.diffTimePoints;
flagTimePoints = sesInfo.flagTimePoints;

%get results from sesInfo file
dispParameters.numOfSess = sesInfo.numOfSess;
dispParameters.numOfSub = sesInfo.numOfSub;
dispParameters.numOfComp = sesInfo.numComp;
dispParameters.inputFiles = sesInfo.inputFiles;
dispParameters.mask_ind = sesInfo.mask_ind;
dispParameters.spmMatrices = sesInfo.userInput.designMatrix;
dispParameters.outputDir = pathstr; % save the output directory
dispParameters.paramFile = paramFile; % parameter file
dispParameters.inputPrefix = sesInfo.userInput.prefix;

try
    dispParameters.zipContents = sesInfo.zipContents;
catch
end

% include spm mat flag
if isfield(sesInfo.userInput, 'spmMatFlag')
    dispParameters.spmMatFlag = sesInfo.userInput.spmMatFlag;
else
    dispParameters.spmMatFlag = 'not_specified';
end

sortingTextFile = [];
handleVisibility = 'on';

if strcmpi(analType, 'batch')
    handleVisibility = 'off';
    % sorting text file (optional way to enter the regressors)
    if isfield(sesInfo.userInput, 'sortingTextFile')
        sortingTextFile = sesInfo.userInput.sortingTextFile;
    else
        sortingTextFile = TXTFILE_REGRESSORS;
    end
end

subjectICAFiles = icatb_parseOutputFiles('icaOutputFiles', sesInfo.icaOutputFiles, 'numOfSub', sesInfo.numOfSub, 'numOfSess', sesInfo.numOfSess);
dispParameters.calibrate_components_mat_file = sesInfo.calibrate_components_mat_file;

clear sesInfo;

dispParameters.flagTimePoints = flagTimePoints;
dispParameters.diffTimePoints = diffTimePoints;

% spm files
spmFiles = dispParameters.spmMatrices;

if isempty(spmFiles(1).name)
    figureData.dispParameters = dispParameters;
    figHandle = figure('visible', 'off', 'userdata', figureData);

    % set the spm matrices
    icatb_selectDesignMatrix(figHandle);
    % figure data
    figureData = get(figHandle, 'userdata');
    clear dispParameters;
    dispParameters = figureData.dispParameters;
    clear figureData;
    delete(figHandle);

end

% subject information
numOfSub = dispParameters.numOfSub;
numOfSess = dispParameters.numOfSess;
spmFiles = dispParameters.spmMatrices;
numComp = dispParameters.numOfComp;

% Results directory to put event averages
eventAvgDir = [dispParameters.inputPrefix, '_event_avg_results'];

spmFiles = str2mat(spmFiles.name);

if isempty(spmFiles)
    error('Spm design matrix is not selected');
end

% Select event average method
choiceStr = str2mat('Regular', 'Deconvolution');
if ~strcmpi(analType, 'batch')
    InputHandle = icatb_getGraphics('Event Average Method', 'normal', 'Event Average Method'); % figure handle
    % set menubar none
    set(InputHandle, 'menubar', 'none', 'windowstyle', 'modal');

    value = icatb_promptUI('popup', 'Which method do you want to use for calculating event average?', str2mat('Regular', 'Deconvolution'), 'numeric', InputHandle);
    delete(InputHandle);
else
    value = 1;
end

if ~exist('eventAvgMethod', 'var')
    % Event average method
    eventAvgMethod = deblank(choiceStr(value, :));
end


if ~exist('selSessions', 'var')

    %%%%%%%%%%%%%%%%%%%%
    %%%% Select Sessions
    %%%%%%%%%%%%%%%%%%%%
    selSessions = 1;
    if numOfSess > 1
        [selSessions] = icatb_listdlg('PromptString', 'Select sessions', 'SelectionMode', 'multiple', ...
            'ListString', str2mat(num2str((1:numOfSess)')), 'movegui', 'center', 'windowStyle', 'modal', 'title_fig', 'Select sessions');
        if isempty(selSessions)
            error('Sessions are not selected for the event average');
        end
    end

end

if ~exist('selSubjects', 'var')

    %%%%%%%%%%%%%%%%%%%%
    %%%% Select Subjects
    %%%%%%%%%%%%%%%%%%%%
    selSubjects = 1;
    if numOfSub > 1
        [selSubjects] = icatb_listdlg('PromptString', 'Select subjects', 'SelectionMode', 'multiple', ...
            'ListString', str2mat(num2str((1:numOfSub)')), 'movegui', 'center', 'windowStyle', 'modal', ...
            'title_fig', 'Select subjects');
        if isempty(selSubjects)
            error('Subjects are not selected for the event average');
        end
    end

end


disp('Loading spm.mat files ...');

if strcmpi(dispParameters.spmMatFlag, 'same_sub_same_sess')
    % same regressors for subjects and sessions
    % loop over each component time course and get the onsets of the regressor
    diffTimePoints = dispParameters.diffTimePoints(1:numOfSess);

    % check the regressor information
    if strcmpi(analType, 'batch')

        % define a new function that takes input from the text file and returns the corresponding variable
        % spmData  or refAppData
        refAppData.data = icatb_pullRegressors_file('txtfile', deblank(sortingTextFile), 'numOfSub', ...
            numOfSub, 'numOfSess', numOfSess, 'keywd', 'single_subject_session', 'spmmatflag', ...
            dispParameters.spmMatFlag, 'spmNames', deblank(spmFiles), 'countTimePoints', ...
            diffTimePoints(1), 'typeSelection', 'multiple', 'data_sessionNumber', 1, 'function_to_use', ...
            'icatb_loadspm_new', 'subjectNumber', []);
        refAppData.data.sel_refnames =  refAppData.data.selectedRegressors;
    else

        refAppData.data = icatb_loadSPM_new('spmName', deblank(spmFiles), 'countTimePoints', ...
            diffTimePoints(1), 'typeStr', 'multiple', 'data_sessionNumber', 1, ...
            'check_spm_design_matrix', 'yes', 'check_regressors', 'yes');

    end
    % end for regressor selection

elseif strcmpi(dispParameters.spmMatFlag, 'same_sub_diff_sess')
    % regressors different over sessions
    diffTimePoints = dispParameters.diffTimePoints(1:numOfSess);

    countS = 0;
    for ii = selSessions
        countS = countS + 1;
        % type for entering regressors
        if strcmpi(analType, 'batch')
            refAppData(countS).data = icatb_pullRegressors_file('txtfile', deblank(sortingTextFile), 'numOfSub', ...
                numOfSub, 'numOfSess', numOfSess, 'keywd', 'single_subject_session', 'spmmatflag', ...
                dispParameters.spmMatFlag, 'spmNames', deblank(spmFiles), 'countTimePoints', diffTimePoints, ...
                'typeSelection', 'multiple', 'data_sessionNumber', ii, 'function_to_use', 'icatb_loadspm_new', ...
                'subjectNumber', 1);
            refAppData(countS).data.sel_refnames =  refAppData(countS).data.selectedRegressors;
        else
            refAppData(countS).data = icatb_loadSPM_new('spmName', deblank(spmFiles), 'countTimePoints', ...
                diffTimePoints, 'typeStr', 'multiple', 'data_sessionNumber', ii, ...
                'check_spm_design_matrix', 'yes', 'check_regressors', 'yes');

        end
        % end for type for entering regressors
    end

elseif strcmpi(dispParameters.spmMatFlag, 'diff_sub_diff_sess')
    % regressors different over subjects and sessions
    count = 0;
    for ii = selSubjects
        % get the time points for the session
        diffTimePoints = dispParameters.diffTimePoints((ii - 1)*numOfSess + 1:ii*numOfSess);
        % loop over sessions
        for jj = selSessions
            count = count + 1;

            % type for entering regressors
            if strcmpi(analType, 'batch')
                refAppData(count).data = icatb_pullRegressors_file('txtfile', deblank(sortingTextFile), ...
                    'numOfSub', numOfSub, 'numOfSess', numOfSess, 'keywd', 'single_subject_session', 'spmmatflag', ...
                    dispParameters.spmMatFlag, 'spmNames', deblank(spmFiles(ii, :)), 'countTimePoints', diffTimePoints, ...
                    'typeSelection', 'multiple', 'data_sessionNumber', jj, 'function_to_use', 'icatb_loadspm_new', ...
                    'subjectNumber', ii);
                refAppData(count).data.sel_refnames =  refAppData(count).data.selectedRegressors;
            else

                refAppData(count).data = icatb_loadSPM_new('spmName', deblank(spmFiles(ii, :)), 'countTimePoints', ...
                    diffTimePoints, 'typeStr', 'multiple', 'data_sessionNumber', jj, ...
                    'check_spm_design_matrix', 'yes', 'check_regressors', 'yes');

            end
            % end for selecting regressors

        end
        % end loop over sessions
    end
    % end loop over subjects
end
% end for checking the state of the design matrix

inputPrefix = dispParameters.inputPrefix;
spmMatFlag = dispParameters.spmMatFlag;

%clear dispParameters;

% Form list string
if strcmpi(spmMatFlag, 'same_sub_diff_sess')
    listString = repmat(struct('name', ''), 1, length(selSessions));
    for ii = 1:length(selSessions)
        listString(ii).name = ['Sess ', num2str(selSessions(ii))];
    end
    temp = str2mat(listString.name); clear listString;
    listString = temp; clear temp;
elseif strcmpi(spmMatFlag, 'same_sub_same_sess')
    listString = 'Selected Data';
else
    listString = repmat(struct('name', ''), 1, length(selSubjects)*length(selSessions));
    listCount = 0;
    for ii = 1:length(selSubjects)
        for jj = 1:length(selSessions)
            listCount = listCount + 1;
            listString(listCount).name = ['Sub ', num2str(selSubjects(ii)), ' Sess ', ...
                num2str(selSessions(jj))];
        end
    end
    temp = str2mat(listString.name); clear listString;
    listString = temp; clear temp;
end
% end for forming list string

% Get the event information from the selected data-sets, regressors
eventAppData = icatb_getInfo_EventAvg(listString, refAppData, numComp, handleVisibility);

drawnow;

if exist('compNumber', 'var')
    eventAppData.component_number = compNumber;
    clear compNumber;
end

% Selected components
component_numbers = [eventAppData.component_number];

% Make a directory
if (exist(fullfile(outputDir, eventAvgDir), 'dir') ~= 7)
    mkdir(outputDir, eventAvgDir);
end

numPlots = 1;
if length(selSubjects)*length(selSessions) > 1
    numPlots = length(selSubjects)*length(selSessions) + 1;
end

% Loop over selected components
for compNumber = component_numbers
    % Loop over selected regressors
    for nRegress = 1:length(eventAppData.eventData(1).regressorInfo)
        fprintf('\n');
        count = 0;
        % loop over selected subjects and sessions for the particular component
        for ii = selSubjects
            countS = 0;
            for jj = selSessions
                count = count + 1;
                countS = countS + 1;
                % check the spm flag
                if strcmpi(spmMatFlag, 'same_sub_same_sess')
                    selectedOnset = eventAppData.eventData(1).regressorInfo(nRegress).selectedOnset;
                    TR = eventAppData.eventData(1).TR;
                    selectedRegressor = eventAppData.eventData(1).regressorInfo(nRegress).selectedRegressor;
                elseif strcmpi(spmMatFlag, 'same_sub_diff_sess')
                    selectedOnset = eventAppData.eventData(countS).regressorInfo(nRegress).selectedOnset;
                    TR = eventAppData.eventData(countS).TR;
                    selectedRegressor =  eventAppData.eventData(countS).regressorInfo(nRegress).selectedRegressor;
                else
                    selectedOnset = eventAppData.eventData(count).regressorInfo(nRegress).selectedOnset;
                    TR = eventAppData.eventData(count).TR;
                    selectedRegressor = eventAppData.eventData(count).regressorInfo(nRegress).selectedRegressor;
                end
                % end for checking the spm flag

                disp(['Calculating Event Average for Subject ', num2str(ii), ' Session ', ...
                    num2str(jj), ' component ', num2str(compNumber)]);

                disp(['Selected regressor for event average is ', selectedRegressor]);

                % Load calibrated component files
                %                 calibratedCompFile = [calibrateMATFileName, num2str(ii), '-', num2str(jj), '.mat'];
                %                 load(fullfile(outputDir, calibratedCompFile), 'tc');

                tempTC = icatb_loadComp(dispParameters, compNumber, 'subjects', ii, 'sessions', jj, 'vars_to_load', 'tc', 'subject_ica_files', subjectICAFiles);

                % smooth ica time course
                if strcmpi(SMOOTHPARA, 'yes')
                    tempTC = icatb_gauss_smooth1D(tempTC, SMOOTHINGVALUE);
                end

                tempTC = icatb_detrend(tempTC, 1, length(tempTC), detrendNumber);

                if strcmpi(eventAvgMethod, 'regular')
                    % Calculate event average
                    eventAvg = icatb_calculate_eventAvg(tempTC, eventInterpFactor, TR, eventWindowSize, selectedOnset);
                    %timeAxis = (1:length(eventAvg)) / eventInterpFactor;
                    timeAxis = icatb_interp((1:eventWindowSize), eventInterpFactor);
                else
                    % Deconvolution method
                    eventAvg = icatb_calculate_eventAvg_deconv(tempTC, TR, eventWindowSize, selectedOnset);
                    timeAxis = (1:length(eventAvg));
                end

                if count == 1
                    eventTc = zeros(length(eventAvg), numPlots);
                    titleStr = cell(numPlots, 1);
                end

                titleStr{count} = ['Event average of sub ', num2str(ii), ' sess ', num2str(jj), ...
                    ' comp ', num2str(compNumber), ' (', selectedRegressor, ')'];
                eventTc(:, count) = eventAvg(:);
            end
            % end loop over sessions
        end
        % end loop over subjects

        if count > 1
            count = count + 1;
            eventTc(:, count) = mean(eventTc(:, 1:count - 1), 2);
            titleStr{count} = ['Mean over selected datasets for comp ', num2str(compNumber), ' (', selectedRegressor, ')'];
        end

        % Saving the event averages
        eventMatFile = fullfile(outputDir, eventAvgDir, [inputPrefix, '_ev_avg_', lower(eventAvgMethod), '_comp_', num2str(compNumber), ...
            '_regress_', num2str(nRegress), '.mat']);

        disp(['Saving event average information in file: ', eventMatFile]);

        icatb_save(eventMatFile, 'eventTc', 'timeAxis', 'titleStr', 'selectedRegressor', 'eventAvgMethod');

        clear eventAvg eventTc;
        fprintf('\n');
    end
    % End loop over selected regressors

end
% End loop over selected components

disp('Done');