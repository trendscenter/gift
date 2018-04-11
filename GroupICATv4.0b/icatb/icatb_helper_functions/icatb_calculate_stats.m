function icatb_calculate_stats
% GUI for calculates stats like mean, std and t-map for
% data-sets selected

% load defaults
icatb_defaults;
global PARAMETER_INFO_MAT_FILE;
global MEAN_ALL_AN3_FILE;
global MEAN_AN3_FILE;
global TMAP_AN3_FILE;
global STD_AN3_FILE;
global COMPONENT_NAMING;
global SESSION_POSTFIX;
global TIMECOURSE_NAMING;
global FUNCTIONAL_DATA_FILTER;


imageExtns = FUNCTIONAL_DATA_FILTER(2:end);

% filter for parameter file
filterP = ['*', PARAMETER_INFO_MAT_FILE, '*.mat'];

% select the parameter file
param_file = icatb_selectEntry('typeEntity', 'file', 'typeSelection', 'single', 'filter', ...
    filterP, 'title', 'Select a valid parameter file');

drawnow;

% output directory
[outputDir, fileName, extn] = fileparts(param_file);

% change directory
cd(outputDir);

% load the parameter file
load(param_file);

if ~exist('sesInfo', 'var')
    error(['Selected file: ', param_file, ' is not a valid parameter file']);
end

% check if the sesInfo is initialized
if sesInfo.isInitialized == 0
    error('stats will be calculated only after you run the analysis.');
end

inputPrefix = sesInfo.userInput.prefix;
meanAllNaming = MEAN_ALL_AN3_FILE;
meanNaming = MEAN_AN3_FILE;
stdNaming = STD_AN3_FILE;
tmapNaming = TMAP_AN3_FILE;

writeInDirs = 0;
try
    writeInDirs = sesInfo.write_analysis_steps_in_dirs;
catch
end

if (writeInDirs == 1)
    meanNaming = fullfile('_group_stats_files', regexprep(meanNaming, '^(_)', ''));
    stdNaming = fullfile('_group_stats_files', regexprep(stdNaming, '^(_)', ''));
    tmapNaming = fullfile('_group_stats_files', regexprep(tmapNaming, '^(_)', ''));
    meanAllNaming = fullfile('_group_stats_files', regexprep(meanAllNaming, '^(_)', ''));
end

% Number of subjects and sessions
numOfSub = sesInfo.numOfSub;
numOfSess = sesInfo.numOfSess;

mask_ind = sesInfo.mask_ind;

% Output directory
sesInfo.outputDir = outputDir;

% Number of data sets is one
if numOfSub*numOfSess == 1
    disp('Stats cannot be calculated as only one data set is present.');
    return;
end


% get the count for time points
if isfield(sesInfo, 'diffTimePoints')
    countFiles = sesInfo.diffTimePoints; % different time points
else
    % get the count of the files
    [countFiles] = icatb_get_countTimePoints(sesInfo.userInput.files);
    sesInfo.userInput.diffTimePoints = countFiles;
    sesInfo.diffTimePoints = sesInfo.userInput.diffTimePoints;
    checkTimePoints = find(countFiles ~= countFiles(1));

    if ~isempty(checkTimePoints)
        sesInfo.flagTimePoints = 'different_time_points';
    else
        sesInfo.flagTimePoints = 'same_time_points';
    end

    icatb_save(param_file, 'sesInfo');
end
% end for checking the time points


sesInfo.subjectICAFiles = icatb_parseOutputFiles('icaOutputFiles', sesInfo.icaOutputFiles, 'numOfSub', sesInfo.numOfSub, 'numOfSess', sesInfo.numOfSess);

[modalityType, dataTitle, compSetFields] = icatb_get_modality;


if isfield(sesInfo, 'modality')
    if ~strcmpi(sesInfo.modality, modalityType)
        if strcmpi(sesInfo.modality, 'fmri')
            error('You have selected the fMRI parameter file. Use GIFT toolbox to calculate stats.');
        else
            error('You have selected the EEG parameter file. Use EEGIFT toolbox to calculate stats.');
        end
    end
end

sesInfo.compSetFields = compSetFields;

if strcmpi(modalityType, 'fmri')
    if ~strcmpi(sesInfo.flagTimePoints, 'same_time_points')
        disp(['Calculating stats for only component signals as weights have different dimensions...']);
        fprintf('\n');
    end
end

statsStr = str2mat('Mean', 'Standard Deviation', 'T Map');

% get the value from the popup box
InputHandle = icatb_getGraphics('Selecting data-sets', 'normal', 'data-sets');

% set menubar none
set(InputHandle, 'menubar', 'none');

% Select subjects and sessions
if (numOfSub == 1 & numOfSess > 1)

    flagStr = 'over_data_sets';

    % subject and session numbers
    subjectNumber = 1;
    sessionNumber = [1:numOfSess];

    % question for calculating
    questionStr = 'What do you want to calculate?';
    choiceStr = statsStr;
    value2 = icatb_promptUI('popup', questionStr, choiceStr, 'numeric', InputHandle);

    if isempty(value2)
        error(['What do you want to calculate is not selected?']);
    end

    if strcmpi(deblank(choiceStr(value2, :)), 'mean')
        funcEval = 'mean';
    elseif strcmpi(deblank(choiceStr(value2, :)), 'standard deviation')
        funcEval = 'std';
    else
        funcEval = 'tmap';
    end
    delete(InputHandle);

elseif (numOfSub > 1 & numOfSess == 1)

    flagStr = 'over_data_sets';

    % subject and session numbers
    subjectNumber = [1:numOfSub];
    sessionNumber = 1;

    % question for calculating
    questionStr = 'What do you want to calculate?';
    choiceStr = statsStr;
    value2 = icatb_promptUI('popup', questionStr, choiceStr, 'numeric', InputHandle);

    if isempty(value2)
        error(['What do you want to calculate is not selected?']);
    end

    if strcmpi(deblank(choiceStr(value2, :)), 'mean')
        funcEval = 'mean';
    elseif strcmpi(deblank(choiceStr(value2, :)), 'standard deviation')
        funcEval = 'std';
    else
        funcEval = 'tmap';
    end
    delete(InputHandle);

else

    questionStr = 'How do you want to select the data-sets?';
    choiceStr = str2mat('Select all datasets', 'Select all sessions for certain subjects', 'Select certain sessions for all subjects');
    value1 = icatb_promptUI('popup', questionStr, choiceStr, 'numeric', InputHandle);

    questionStr = 'What do you want to calculate?';
    choiceStr = statsStr;
    value2 = icatb_promptUI('popup', questionStr, choiceStr, 'numeric', InputHandle);

    delete(InputHandle);

    if isempty(value1)
        error('Data-sets are not specified');
    end

    if isempty(value2)
        error(['What do you want to calculate is not selected?']);
    end

    if value1 == 1
        % all data sets
        flagStr = 'over_data_sets';
        subjectNumber = [1:numOfSub];
        sessionNumber = [1:numOfSess];
        %
    elseif value1 == 2
        % over sessions
        flagStr = 'over_sessions';
        sessionNumber = [1:numOfSess];
        for ii = 1:numOfSub
            subStr(ii).name = ['Subject ', num2str(ii)];
        end
        % select subjects
        subjectNumber = icatb_listdlg('PromptString', 'Select Subject/Subjects', 'SelectionMode', 'multiple', ...
            'ListString', str2mat(subStr.name), 'movegui', 'center', 'windowStyle', 'modal', ...
            'title_fig', 'Select subjects');
        if isempty(subjectNumber)
            error('Subject/subjects are not selected');
        end
    else
        % over subjects
        flagStr = 'over_subjects';
        subjectNumber = [1:numOfSub];
        % select sessions
        for ii = 1:numOfSess
            sessStr(ii).name = ['Session ', num2str(ii)];
        end
        % select sessions
        sessionNumber = icatb_listdlg('PromptString', 'Select Session/Sessions', 'SelectionMode', 'multiple', ...
            'ListString', str2mat(sessStr.name), 'movegui', 'center', 'windowStyle', 'modal', ...
            'title_fig', 'Select sessions');
        if isempty(sessionNumber)
            error('Session/sessions are not selected');
        end
    end
    % end for checking the type of data-sets


    if strcmpi(deblank(choiceStr(value2, :)), 'mean')
        funcEval = 'mean';
    elseif strcmpi(deblank(choiceStr(value2, :)), 'standard deviation')
        funcEval = 'std';
    else
        funcEval = 'tmap';
    end

end
% end for checking data-sets

drawnow;

numComp = sesInfo.numComp;

helpHandle = helpdlg('Calculating Stats. Please wait ...', 'Calculating Stats');

if strcmpi(modalityType, 'fmri')

    % read the volume of the first functional image
    [V, HInfo] = icatb_returnHInfo(deblank(sesInfo.inputFiles(1).name(1, :)));

else

    HInfo = [];

end

switch lower(flagStr)

    case 'over_data_sets'

        disp(['Calculating mean over data-sets ...']);

        % calculate mean
        [meanIm, meanA] = calculate_mean(subjectNumber, sessionNumber, sesInfo);

        % calculate std
        if strcmpi(funcEval, 'std') | strcmpi(funcEval, 'tmap')
            disp(['Calculating std over data-sets ...']);
            [stdIm, stdA] = calculate_std(subjectNumber, sessionNumber, sesInfo, meanIm, meanA);
        end
        % end for calculating std

        % calculate tmap
        if strcmpi(funcEval, 'tmap')
            disp(['Calculating tmap over data-sets ...']);
            [tmapIm] = calculate_tmap(meanIm, stdIm, length(subjectNumber)*length(sessionNumber));
        end
        % end for calculating tmap

        % write the data
        if strcmpi(funcEval, 'mean')

            fileNaming = [inputPrefix, meanAllNaming, COMPONENT_NAMING, SESSION_POSTFIX, '_all_'];
            disp(['Writing mean components as ', fileNaming]);
            fileNaming = [fileNaming, imageExtns];
            % write data
            icatb_saveICAData(fileNaming, meanIm, meanA, mask_ind, numComp, HInfo, 'real', [], outputDir);

        elseif strcmpi(funcEval, 'std')

            fileNaming = [inputPrefix, stdNaming, COMPONENT_NAMING, SESSION_POSTFIX, '_all_'];
            % write data
            disp(['Writing std components as ', fileNaming]);
            fileNaming = [fileNaming, imageExtns];
            icatb_saveICAData(fileNaming, stdIm, stdA, mask_ind, numComp, HInfo, 'real', [], outputDir);

        elseif strcmpi(funcEval, 'tmap')

            fileNaming = [inputPrefix, tmapNaming, COMPONENT_NAMING, SESSION_POSTFIX, '_all_'];
            % write data
            disp(['Writing tmaps components ', fileNaming]);
            fileNaming = [fileNaming, imageExtns];
            icatb_saveICAData(fileNaming, tmapIm, meanA, mask_ind, numComp, HInfo, 'real', [], outputDir);

        end
        % end for writing the data



    case 'over_sessions'

        % loop over selected subjects
        for ii = 1:length(subjectNumber)

            disp(['Calculating mean over sessions for subject ', num2str(subjectNumber(ii))]);

            % calculate mean
            [meanIm, meanA] = calculate_mean(subjectNumber(ii), sessionNumber, sesInfo);

            % calculate std
            if strcmpi(funcEval, 'std') | strcmpi(funcEval, 'tmap')
                disp(['Calculating std over sessions for subject ', num2str(subjectNumber(ii))]);
                [stdIm, stdA] = calculate_std(subjectNumber(ii), sessionNumber, sesInfo, meanIm, meanA);
            end
            % end for calculating std

            % calculate tmap
            if strcmpi(funcEval, 'tmap')
                disp(['Calculating tmaps over sessions for subject ', num2str(subjectNumber(ii))]);
                [tmapIm] = calculate_tmap(meanIm, stdIm, length(sessionNumber));
            end
            % end for calculating tmap

            % write the data
            if strcmpi(funcEval, 'mean')

                fileNaming = [inputPrefix, meanNaming, '_sub_', icatb_returnFileIndex(subjectNumber(ii)), ...
                    COMPONENT_NAMING];
                % write data
                disp(['Writing mean components for subject ', num2str(subjectNumber(ii)), ' as ', fileNaming]);
                fileNaming = [fileNaming, imageExtns];
                icatb_saveICAData(fileNaming, meanIm, meanA, mask_ind, numComp, HInfo, 'real', [], outputDir);

            elseif strcmpi(funcEval, 'std')

                fileNaming = [inputPrefix, stdNaming, '_sub_', icatb_returnFileIndex(subjectNumber(ii)), ...
                    COMPONENT_NAMING];
                % write data
                disp(['Writing std components for subject ', num2str(subjectNumber(ii)), ' as ', fileNaming]);
                fileNaming = [fileNaming, imageExtns];
                icatb_saveICAData(fileNaming, stdIm, stdA, mask_ind, numComp, HInfo, 'real', [], outputDir);

            elseif strcmpi(funcEval, 'tmap')

                fileNaming = [inputPrefix, tmapNaming, '_sub_', icatb_returnFileIndex(subjectNumber(ii)), ...
                    COMPONENT_NAMING];
                % write data
                disp(['Writing tmaps for subject ', num2str(subjectNumber(ii)), ' as ', fileNaming]);
                fileNaming = [fileNaming, imageExtns];
                icatb_saveICAData(fileNaming, tmapIm, meanA, mask_ind, numComp, HInfo, 'real', [], outputDir);

            end
            % end for writing the data

        end
        % end loop over subjects

    case 'over_subjects'

        % loop over selected sessions
        for ii = 1:length(sessionNumber)

            disp(['Calculating mean over subjects for session ', num2str(sessionNumber(ii))]);

            % calculate mean
            [meanIm, meanA] = calculate_mean(subjectNumber, sessionNumber(ii), sesInfo);

            % calculate std
            if strcmpi(funcEval, 'std') | strcmpi(funcEval, 'tmap')
                disp(['Calculating std over subjects for session ', num2str(sessionNumber(ii))]);
                [stdIm, stdA] = calculate_std(subjectNumber, sessionNumber(ii), sesInfo, meanIm, meanA);
            end
            % end for calculating std

            % calculate tmap
            if strcmpi(funcEval, 'tmap')
                disp(['Calculating tmaps over subjects for session ', num2str(sessionNumber(ii))]);
                [tmapIm] = calculate_tmap(meanIm, stdIm, length(subjectNumber));
            end
            % end for calculating tmap


            % write the data
            if strcmpi(funcEval, 'mean')

                fileNaming = [inputPrefix, meanNaming, '_sess_', icatb_returnFileIndex(sessionNumber(ii)), ...
                    COMPONENT_NAMING];
                % write data
                disp(['Writing mean components over subjects for session ', num2str(sessionNumber(ii)), ' as ', ...
                    fileNaming]);
                fileNaming = [fileNaming, imageExtns];
                icatb_saveICAData(fileNaming, meanIm, meanA, mask_ind, numComp, HInfo, 'real', [], outputDir);

            elseif strcmpi(funcEval, 'std')

                fileNaming = [inputPrefix, stdNaming, '_sess_', icatb_returnFileIndex(sessionNumber(ii)), ...
                    COMPONENT_NAMING];
                % write data
                disp(['Writing std components over subjects for session ', num2str(sessionNumber(ii)), ' as ', ...
                    fileNaming]);
                fileNaming = [fileNaming, imageExtns];
                icatb_saveICAData(fileNaming, stdIm, stdA, mask_ind, numComp, HInfo, 'real', [], outputDir);

            elseif strcmpi(funcEval, 'tmap')

                fileNaming = [inputPrefix, tmapNaming, '_sess_', icatb_returnFileIndex(sessionNumber(ii)), ...
                    COMPONENT_NAMING];
                % write data
                disp(['Writing tmaps over subjects for session ', num2str(sessionNumber(ii)), ' as ', fileNaming]);
                fileNaming = [fileNaming, imageExtns];
                icatb_saveICAData(fileNaming, tmapIm, meanA, mask_ind, numComp, HInfo, 'real', [], outputDir);

            end
            % end for writing the data

        end
        % end for loop over sessions

end
% end for switch

try
    delete(helpHandle);
catch
end

% display the string
disp(['Components are written to directory: ', outputDir]);


function [meanIm, meanA] = calculate_mean(subjectNumber, sessionNumber, sesInfo)
% mean image

compSetFields = sesInfo.compSetFields;

timePoints = min(sesInfo.diffTimePoints);

count = 0;
for ii = 1:length(subjectNumber)
    for jj = 1:length(sessionNumber)

        count = count + 1;
        [tc, ic] = icatb_loadComp(sesInfo, (1:sesInfo.numComp), 'subjects', subjectNumber(ii), 'sessions', sessionNumber(jj), 'subject_ica_files', sesInfo.subjectICAFiles);
        ic = ic';

        if count == 1
            meanIm = zeros(size(ic));
            meanA = zeros(size(tc(1:timePoints, :)));
        end
        meanIm = ic + meanIm;
        % calculate meanA
        %if strcmpi(flagTimePoints, 'same_time_points')
        meanA = meanA + tc(1:timePoints, :);
        %end
        clear ic; clear tc;
    end
end
% calculate the mean
meanIm = meanIm ./ count;
meanA = meanA ./ count;


function [stdIm, stdA] = calculate_std(subjectNumber, sessionNumber, sesInfo, meanIm, meanA)
% standard deviation

compSetFields = sesInfo.compSetFields;

count = 0;
stdIm = zeros(size(meanIm));
stdA = zeros(size(meanA));

timePoints = min(sesInfo.diffTimePoints);

for ii = 1:length(subjectNumber)
    for jj = 1:length(sessionNumber)

        count = count + 1;
        [tc, ic] = icatb_loadComp(sesInfo, (1:sesInfo.numComp), 'subjects', subjectNumber(ii), 'sessions', sessionNumber(jj), 'subject_ica_files', sesInfo.subjectICAFiles);
        ic = ic';

        stdIm = ((ic - meanIm).^2) + stdIm;
        % calculate stdA
        %if strcmpi(flagTimePoints, 'same_time_points')
        stdA = ((tc(1:timePoints, :) - meanA).^2) + stdA;
        %end
        clear ic; clear tc;
    end
end
% standard deviation map
stdIm = sqrt(stdIm./(count - 1));
stdA = sqrt(stdA./(count - 1));


function [tmapIm] = calculate_tmap(meanIm, stdIm, count)
% calculate t-map

tmapIm = zeros(size(meanIm));
for comp = 1:size(tmapIm, 1)
    tmap_ind = find(stdIm(comp, :) > eps);
    divisor = squeeze(stdIm(comp, tmap_ind)) ./ sqrt(count - 1);
    tmapIm(comp, tmap_ind) = squeeze(meanIm(comp, tmap_ind)) ./ divisor;
end