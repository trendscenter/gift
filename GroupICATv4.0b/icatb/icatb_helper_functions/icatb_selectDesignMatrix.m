function icatb_selectDesignMatrix(handles)
% set the spm matrices

figureData = get(handles, 'userdata');

dispParameters = figureData.dispParameters;
% count for time points information
diffTimePoints = dispParameters.diffTimePoints;
% number of subjects and sessions
numOfSub = dispParameters.numOfSub;
numOfSess = dispParameters.numOfSess;
% parameter file
paramFile = dispParameters.paramFile;
% spm matrices
spmMatrices = dispParameters.spmMatrices;
% flag for the time points
flagTimePoints = dispParameters.flagTimePoints;

% load the parameter file
load(paramFile); % (contains sesInfo structure)
subjectFile = [sesInfo.userInput.prefix, 'Subject.mat'];
outputDir = fileparts(paramFile); % get the target directory

% select the state of the design matrix here
questionString = 'How do you want to select the regressors for subjects and sessions?';
optionsString = str2mat('One', 'All');

D(1).string = 'These options are for sorting components temporally. Each option is explained below:';
D(size(D, 2) + 1).string = '';
% D(size(D, 2) + 1).string = 'No - Design matrix is not selected for temporal sorting.';
% D(size(D, 2) + 1).string = '';
D(size(D, 2) + 1).string = 'Same Set of regressors for all subjects and sessions - One SPM2/SPM5/SPM8 model model is used for all the subjects and sessions. Same set of regressors will be applied to all data-sets.';
D(size(D, 2) + 1).string = '';
D(size(D, 2) + 1).string = 'Different set of regressors for all sessions - One SPM2/SPM5/SPM8 model is used for all subjects and sessions. Session specific regressors will be pulled during temporal sorting.';
D(size(D, 2) + 1).string = '';
D(size(D, 2) + 1).string = 'Different set of regressors for all subjects and sessions - One SPM2/SPM5/SPM8 model is used for each subject''s session. Regressors pulled during temporal sorting will be specific to subjects''s session.';

% help details
helpDetails.str = str2mat(D.string);
clear D;
helpDetails.title = 'Options for temporal sorting';
% Options include:
% 1. Same set of regressors for all subjects and sessions (One SPM.mat
% needs to be selected and the selected regressors will be correlated
% with each data-set).
% 2. Session specific regressors will be correlated with the
% corresponding time courses.
% 3. Regressors will be specific to each data-set.

% if the number of data sets is greater than one
if numOfSub*numOfSess > 1
    % if the number of sessions is greater than one
    if numOfSub == 1 & numOfSess > 1

        D(1).string = 'Same set of regressors for all subjects and sessions';
        D(length(D) + 1).string = 'Different set of regressors for sessions';
        % if the number of subjects is greater than one
    elseif numOfSess == 1 & numOfSub > 1
        D(1).string = 'Same set of regressors for all subjects and sessions';
        D(length(D) + 1).string = 'Different set of regressors for subjects and sessions';
    else
        % show all the three cases
        D(1).string = 'Same set of regressors for all subjects and sessions';
        D(length(D) + 1).string = 'Different set of regressors for sessions';
        D(length(D) + 1).string = 'Different set of regressors for subjects and sessions';
    end

    % get the value from the popup box
    InputHandle = icatb_getGraphics('Selecting Design Matrix', 'normal', 'SPM2/SPM5/SPM8 Design Matrix'); % figure handle
    % set menubar none
    set(InputHandle, 'menubar', 'none');
    giftHelpH = uimenu('parent', InputHandle, 'label', 'GIFT-Help');
    helpMenuH = uimenu('parent', giftHelpH, 'label', 'Design Matrix', 'callback', 'icatb_openHTMLHelpFile(''icatb_SPMFileSelection.htm'');');

    value = [];

    choiceString = str2mat(D.string);
    clear D;
    value = icatb_promptUI('popup', questionString, choiceString, 'numeric', InputHandle);
    close(InputHandle);

else
    % for only one subject and one session
    value = 1;
    choiceString = 'Same set of regressors for all subjects and sessions';
end

% file Info format
fileInfo.format = 'append';
fileInfo.fileName = fullfile(outputDir, subjectFile);

% display a dialog box showing the changed information
D(1).string = 'Design matrix is changed and saved in files: ';
D(size(D, 2) + 1).string = '';
D(size(D, 2) + 1).string = ['1. ', fullfile(outputDir, subjectFile)];
D(size(D, 2) + 1).string = '';
D(size(D, 2) + 1).string = ['2. ', paramFile];

% store the strings
newString = str2mat(D.string); % (string will be displayed when the design matrix information is changed)
if ~isempty(value)
    %answerString = deblank(optionsString(value, :));
    selectedString = deblank(choiceString(value, :));
    selectedString = lower(selectedString);

    if strcmp(selectedString, 'same set of regressors for all subjects and sessions')
        answerString = 'same_sub_same_sess';
    elseif strcmp(selectedString, 'different set of regressors for sessions')
        answerString = 'same_sub_diff_sess';
    elseif strcmp(selectedString, 'different set of regressors for subjects and sessions')
        answerString = 'diff_sub_diff_sess';
    else
        error('Unknown option');
    end

    switch lower(answerString)
        case 'same_sub_same_sess'
            if strcmp(lower(flagTimePoints), 'same_time_points')
                %subject_string.name = 'Design matrix for all data sets';
                [SPMFiles.name] = icatb_selectEntry('title', 'Select SPM2/SPM5/SPM8 design matrix', 'filter', '*.mat', ...
                    'typeEntity', 'file', 'typeSelection', 'single');

                % do spm checking
                icatb_loadSPM_new('spmName', SPMFiles.name, 'countTimePoints', ...
                    diffTimePoints(1), 'check_spm_design_matrix', 'yes', 'data_sessionNumber', 1);

                %[SPMFiles] = icatb_select_fMRIData('title', 'Select SPM2 design matrix', 'num_subjects', 1, ...
                %                 'num_sessions', 1, 'files_specification', 'equal', 'spm_check', 'yes', ...
                %                     'filter_string', '*.mat', 'type_file_selection', 'single', 'diffTimePoints', ...
                %                     diffTimePoints(1), 'startPath', outputDir, ...
                %                     'fileInfo', fileInfo, 'subject_string', subject_string, 'figure_menu', 'spm');
            else
                errorStr = ['The fMRI data have different number of images ', ...
                    'Thus identical regressors are not possible.'];
                % open the error dialog box
                error(errorStr);
            end

        case 'same_sub_diff_sess'
            % special case for handling different set of regressors over
            % sessions
            subject_string.name = 'Design matrix for all data sets';

            prevSession = zeros(1, numOfSess);
            currentSession = zeros(1, numOfSess);
            data_sets_count = 0;
            % loop over subjects
            for ii = 1:numOfSub
                % loop over sessions
                for jj = 1:numOfSess
                    data_sets_count = data_sets_count + 1;
                    currentSession(jj) = diffTimePoints(data_sets_count);
                end

                if ii == 1
                    prevSession = currentSession;
                end

                % end for sessions
                % compare the session time points
                if ii > 1
                    check = find(prevSession == currentSession);
                    if length(check) ~= numOfSess
                        error(['Time points or number of images are different over subjects.', ...
                            ' Select the option different regressors for subjects and sessions.']);
                    end
                end

            end
            % end for subjects

            [SPMFiles.name] = icatb_selectEntry('title', 'Select SPM2/SPM5/SPM8 design matrix', 'filter', '*.mat', ...
                'typeEntity', 'file', 'typeSelection', 'single');


            % loop over sessions
            for nSess = 1:numOfSess
                % check spm design matrix
                icatb_loadSPM_new('spmName', SPMFiles.name, 'countTimePoints', ...
                    diffTimePoints(1:numOfSess), 'check_spm_design_matrix', 'yes', 'data_sessionNumber', nSess);
            end
            % end for loop over sessions

            %             [SPMFiles] = icatb_select_fMRIData('title', 'Select SPM2 design matrix', 'num_subjects', 1, ...
            %                 'num_sessions', numOfSess, 'files_specification', 'equal', 'spm_check', 'yes', ...
            %                 'filter_string', '*.mat', 'type_file_selection', 'single', 'diffTimePoints', ...
            %                 prevSession, 'startPath', outputDir, 'fileInfo', fileInfo, 'subject_string', ...
            %                 subject_string, 'figure_menu', 'spm');

        case 'diff_sub_diff_sess'

            tempSPMFiles = icatb_selectEntry('title', 'Select SPM2/SPM5/SPM8 design matrix/matrices', 'filter', '*.mat', ...
                'typeEntity', 'file', 'typeSelection', 'multiple');

            drawnow;

            if (size(tempSPMFiles, 1) ~= numOfSub)
                error('Error:SPMDesign', 'Number of SPM.mat files (%d) doesn''t match the number of subjects (%d)\n', size(tempSPMFiles, 1), ...
                    numOfSub);
            end

            SPMFiles = repmat(struct('name', []), 1, numOfSub);
            for nSub = 1:numOfSub
                SPMFiles(nSub).name = deblank(tempSPMFiles(nSub, :));
            end

            %             for ii = 1:numOfSub
            %                 subject_string(ii).name = ['Subject ', num2str(ii)];
            %             end
            %             [SPMFiles] = icatb_select_data('title', 'Select SPM2/SPM5 design matrix', 'num_subjects', numOfSub, ...
            %                 'num_sessions', numOfSess, 'files_specification', 'equal', 'spm_check', 'yes', ...
            %                 'filter_string', '*.mat', 'type_file_selection', 'single', 'diffTimePoints', diffTimePoints, ...
            %                 'startPath', outputDir, 'fileInfo', fileInfo, 'subject_string', subject_string, 'figure_menu', 'spm');
    end

    try
        icatb_save(fullfile(outputDir, subjectFile), '-append', 'SPMFiles');

    catch
        % store the strings
        newString = ['Design matrix is changed and saved in file: ', paramFile]; % (string will be displayed when the design matrix information is changed)
    end

    % remove the design matrix field
    sesInfo.userInput = rmfield(sesInfo.userInput, 'designMatrix');
    % set the sesInfo user Input
    sesInfo.userInput.designMatrix = SPMFiles;
    sesInfo.userInput.spmMatFlag = answerString;
    icatb_save(paramFile, 'sesInfo');
    % set the figure data
    figureData.dispParameters.spmMatrices = SPMFiles;
    figureData.dispParameters.spmMatFlag = answerString;
    set(handles, 'userdata', figureData);
    % set the design matrix
    %designMatrixData.spmMatrices = SPMFiles;
    disp(newString);
else
    disp('Not changing the design matrix');
    return;
end