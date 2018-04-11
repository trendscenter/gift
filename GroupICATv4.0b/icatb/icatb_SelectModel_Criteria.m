function [modelTimecourse, nRegress, num_DataSets, refInfo] = icatb_SelectModel_Criteria(varargin)
% Purpose: Select model based on the regressor information, sorting criteria. 
% Inputs must be in pairs:
%
% Inputs:
% 1. statedesignmatrix - Options are 'same_sub_same_sess',
% 'same_sub_diff_sess' and 'diff_sub_diff_sess'.
% 2. keywd - Type of data-set like 'single_subject_session',
% 'single_session', 'single_subject', 'all_subjects_sessions'
% 3. numOfSess - Number of sessions
% 4. numOfSub - Number of subjects
% 5. sortingCriteria - 'Multiple regression', 'Maximum voxel',
% 'Correlation' and 'Kurtosis'
% 6. H - Figure handle
% 7. spmMatrices - SPM matrices
% 8. 'selected_design_matrix' - Selected design index for SPM matrices
% 9. 'difftimepoints' - 'same_time_points' or 'different_time_points'
% 10. 'counttimepoints' - Vector for time points of the selected subjects and
% sessions
% 11. 'all_time_points' - structure with field sub and sub with field sess
% 12. 'single_subject_session_number' - session number (important variable to check time
% points with the SPM design matrix).
% 13. 'sortingtextfile' - Text file used for sorting
%
% Output:
% 1. modelTimecourse - Model time course concatenated to equal the ICA time
% course
% 2. nRegress - Number of regressors
% 3. num_DataSets - Number of data sets
% 4. refInfo - refInfo structure containing SPM data


icatb_defaults;
global METHOD_ENTERING_REGRESSORS;
global TXTFILE_REGRESSORS;

% pass the variables in pairs
if mod(nargin, 2) ~= 0
    error('Arguments must be in pairs');
end

% Initialise session Index
sessionIndex = [];

% Loop over number of arguments
for ii = 1:2:nargin
    % get the state of the design matrix:
    % (same set of regressors over all data-sets, different regressors over
    % sessions, different regressors for all subjects and sessions)
    if strcmp(lower(varargin{ii}), 'statedesignmatrix')
        stateDesignMatrix = varargin{ii + 1};

        % get the key word (all datasets, single subject, single session,
        % single subject single session)
    elseif strcmp(lower(varargin{ii}), 'keywd')
        keywd = varargin{ii + 1};

        % original number of sessions
    elseif strcmp(lower(varargin{ii}), 'numofsess')
        numOfSess = varargin{ii + 1};
        % original number of subjects
    elseif strcmp(lower(varargin{ii}), 'numofsub')
        numOfSub = varargin{ii + 1};
        % sorting criteria
    elseif strcmp(lower(varargin{ii}), 'sortingcriteria')
        sortingCriteria = varargin{ii + 1};
        % pass the graphics handle
    elseif strcmp(lower(varargin{ii}), 'graphicshandle')
        H = varargin{ii + 1};
        % spm matrices file names
    elseif(strcmp(lower(varargin{ii}),'spmmatrices'))
        spmMatrices = varargin{ii + 1};
        % selected design matrix
    elseif strcmp(lower(varargin{ii}), 'selected_design_matrix')
        selectedDesignIndex = varargin{ii + 1};
        % flag for time points
    elseif strcmp(lower(varargin{ii}), 'difftimepoints')
        flagTimePoints = varargin{ii + 1};
        % time points information for all subjects and sessions
    elseif strcmp(lower(varargin{ii}), 'counttimepoints')
        countTimePoints = varargin{ii + 1};
        % time points information for all subjects and sessions
    elseif strcmp(lower(varargin{ii}), 'all_time_points')
        allTimePoints = varargin{ii + 1};
        % pass the session number
    elseif strcmp(lower(varargin{ii}), 'single_subject_session_number')
        sessionIndex = varargin{ii + 1};
        % pass the batch text file for regression or correlation
    elseif strcmp(lower(varargin{ii}), 'sortingtextfile')
        sortingTextFile = varargin{ii + 1};
    end

end
% end for checking the arguments

% state of the design matrix
if ~exist('stateDesignMatrix', 'var')
    stateDesignMatrix = 'no';
end

% By default the component set is one subject and one session
if ~exist('keywd', 'var')
    keywd = 'single_subject_session';
end

% number of sessions or subjects  doesn''t exist
if ~exist('numOfSess', 'var') | ~exist('numOfSub', 'var')
    error('Number of subjects or sessions doesn''t exist');
end

subStr = 'many';

% type of selection of list dialog
typeStr = 'single';

% Multiple regressors have the option for selecting more than one regressor
if strcmp(sortingCriteria, 'multiple regression')
    typeStr = 'multiple';
end

if ~exist('flagTimePoints', 'var')
    flagTimePoints = 'same_time_points';
end

if strcmp(lower(stateDesignMatrix), 'no')
    error('Please select the design matrix under menu Display GUI Options in displayGUI');
end

%% Concatenate model time courses here according to the one model
%% for all subjects or one model for each subject

refInfo.modelIndex = [];

refInfo.SPMFile.name = [];

refInfo.numSubjects = numOfSub;

refInfo.numSessions = numOfSess;

stateDesignMatrix = lower(stateDesignMatrix);

refInfo.spmMatFlag = stateDesignMatrix;

analType = lower(METHOD_ENTERING_REGRESSORS);

% check if the regressors are automatically selected
if strcmpi(analType, 'automatic')
    analType = 'batch';
end
    

% answer to the question is batch or gui
if strcmpi(analType, 'batch')
    if isempty(sortingTextFile)
        sortingTextFile = TXTFILE_REGRESSORS;
    end

    % text file for sorting
    refInfo.sortingTextFile = sortingTextFile;
end
% end for answer to the question

% if the text file is present for entering the regressors then read the
% information from the file and close the file


% Determine the type of concatentaion
switch stateDesignMatrix
    % One design matrix is selected
    case 'same_sub_same_sess'
        % design information

        % check the regressor information
        if strcmpi(analType, 'batch')

            % define a new function that takes input from the text file and returns the corresponding variable
            % spmData  or refAppData
            [spmData] = icatb_pullRegressors_file('txtfile', deblank(refInfo.sortingTextFile), 'numOfSub', numOfSub, ...
                'numOfSess', numOfSess, 'keywd', keywd, 'spmmatflag', stateDesignMatrix, ...
                'spmNames', spmMatrices(1).name, 'countTimePoints', countTimePoints(1), 'typeSelection', typeStr, ...
                'data_sessionNumber', 1, 'function_to_use', 'icatb_loadspm_new', 'subjectNumber', []);

        else

            % load model time courses and also check the SPM design matrix
            [spmData] = icatb_loadSPM_new('spmName', spmMatrices(1).name, 'countTimePoints', ...
                countTimePoints(1), 'typeStr', typeStr, 'data_sessionNumber', 1, ...
                'check_spm_design_matrix', 'yes', 'check_regressors', 'yes', 'flag_selecting_regressors', ...
                'yes', 'get_timecourses', 'yes');

        end
        % end for checking the regressor information

        % Reference Function Information
        refInfo.SPMFile.name = spmMatrices(1).name;
        refInfo.modelIndex = spmData.getIndex;
        % selected regressor information
        refInfo.selectedRegressors.name = str2mat(spmData.selectedRegressors);
        refInfo.spmAppData.data = spmData; %
        timecourse = spmData.timecourse; % model time course
        clear spmData;
        % number of regressors
        nRegress = size(timecourse, 2);
        % Initialise model time course
        modelTimecourse = zeros(sum(countTimePoints), nRegress);

        switch keywd

            case 'single_subject_session'
                % A particular subject session

                modelTimecourse = detrend(timecourse, 0); % (mean is removed)

                num_DataSets = 1; % (selected number of datasets)

            case 'single_session'
                % initialise begin
                begin = 1;

                % loop over subjects
                for numSub = 1:numOfSub
                    % loop over number of regressors
                    for numRegress = 1:nRegress
                        % concatentae model time course
                        modelTimecourse(begin:numSub*size(timecourse, 1), numRegress) = ...
                            detrend(timecourse(:, numRegress), 0);
                    end
                    % end loop over regressors
                    begin = begin + size(timecourse, 1); % update the number
                end

                num_DataSets = numOfSub; % return the number of data sets

            case 'single_subject'

                begin = 1;
                % loop over sessions
                for numSess = 1:numOfSess
                    % loop over regressors
                    for numRegress = 1:nRegress
                        % concatenate model time course
                        modelTimecourse(begin:numSess*size(timecourse, 1), numRegress) = ...
                            detrend(timecourse(:, numRegress), 0);
                    end
                    % end loop over regressors
                    begin = begin + size(timecourse, 1); % update the number
                end
                % end loop over sessions

                num_DataSets = numOfSess; % return the number of data sets

            case 'all_subjects_sessions'

                begin = 1;
                % loop over data sets
                for numDataSets = 1:numOfSub*numOfSess
                    % loop over regressors
                    for numRegress = 1:nRegress
                        modelTimecourse(begin:numDataSets*size(timecourse, 1), numRegress) = ...
                            detrend(timecourse(:, numRegress), 0);
                    end
                    % end loop over regressors
                    begin = begin + size(timecourse, 1); % update the number
                end
                % end loop over data sets

                num_DataSets = numOfSub*numOfSess; % return the number of data sets

        end
        % end for switch containing one design matrix selection


        % additional case deals with only one SPM.mat and different
        % sessions. There are four cases how the time courses can be
        % concatenated
    case 'same_sub_diff_sess'
        % combine both single subject and single session as well as

        % check the key word
        if strcmp(keywd, 'single_subject_session')

            % specific subject session time course will be correlated
            % with the model time course (pull session specific
            % regressors)

            [normalPos] = icatb_getfigPosition('normal');
            % width and height of the list box
            listsize = [0.9*normalPos(3) 0.95*normalPos(4)];

            % make the session index empty if it is zero
            if sessionIndex == 0
                sessionIndex = [];
            end

            % pull time points for a subject
            tempTimePoints = zeros(1, numOfSess);
            for nn = 1:numOfSess
                tempTimePoints(nn) = allTimePoints.sub(1).sess(nn);
            end

            if isempty(sessionIndex)
                tempTimePoints = tempTimePoints(1);
            end

            % check the regressor information
            if strcmpi(analType, 'batch')

                [spmData] = icatb_pullRegressors_file('txtfile', deblank(refInfo.sortingTextFile), 'numOfSub', ...
                    numOfSub, 'numOfSess', numOfSess, 'keywd', keywd, 'spmmatflag', stateDesignMatrix, ...
                    'spmNames', spmMatrices.name, 'countTimePoints', tempTimePoints, 'typeSelection', typeStr, ...
                    'data_sessionNumber', sessionIndex, 'function_to_use', 'icatb_loadspm_new', 'subjectNumber', 1);
            else

                % check the SPM design matrix
                [spmData] = icatb_loadSPM_new('spmName', spmMatrices.name, 'countTimePoints', ...
                    tempTimePoints, 'typeStr', typeStr, 'data_sessionNumber', sessionIndex, ...
                    'check_spm_design_matrix', 'yes', 'check_regressors', 'yes', 'flag_selecting_regressors', ...
                    'yes', 'get_timecourses', 'yes');

            end
            % end for checking batch or gui

            clear tempTimePoints;

            % Reference Function Information
            refInfo.SPMFile.name = spmMatrices.name;
            refInfo.modelIndex = spmData.getIndex;
            % selected regressors information
            refInfo.selectedRegressors.name = str2mat(spmData.selectedRegressors);
            refInfo.spmAppData.data = spmData;
            timecourse = spmData.timecourse; % get the time course

            % number of regressors
            nRegress = size(timecourse, 2);


            modelTimecourse = detrend(timecourse, 0);

            num_DataSets = 1; % return the number of data sets

        elseif strcmp(keywd, 'single_session')

            sessionNumber = selectedDesignIndex; % session number

            % pull time points for a subject
            tempTimePoints = zeros(1, numOfSess);
            for nn = 1:numOfSess
                tempTimePoints(nn) = allTimePoints.sub(1).sess(nn);
            end

            % check regressor information
            if strcmpi(analType, 'batch')

                [spmData] = icatb_pullRegressors_file('txtfile', deblank(refInfo.sortingTextFile), 'numOfSub', numOfSub, ...
                    'numOfSess', numOfSess, 'keywd', keywd, 'spmmatflag', stateDesignMatrix, ...
                    'spmNames', spmMatrices.name, 'countTimePoints', tempTimePoints, 'typeSelection', typeStr, ...
                    'data_sessionNumber', sessionNumber, 'function_to_use', 'icatb_loadspm_new', 'subjectNumber', 1);

            else

                % check the SPM design matrix
                [spmData] = icatb_loadSPM_new('spmName', spmMatrices.name, 'countTimePoints', ...
                    tempTimePoints, 'typeStr', typeStr, 'data_sessionNumber', sessionNumber, ...
                    'check_spm_design_matrix', 'yes', 'check_regressors', 'yes', ...
                    'flag_selecting_regressors', 'yes', 'get_timecourses', 'yes');

            end
            % end for checking regressor information

            % Reference Function Information
            refInfo.SPMFile.name = spmMatrices.name;
            refInfo.modelIndex = spmData.getIndex;
            % selected regressors information
            refInfo.selectedRegressors.name = str2mat(spmData.selectedRegressors);
            refInfo.spmAppData.data = spmData;
            timecourse = spmData.timecourse; % get the time course

            % number of regressors
            nRegress = size(timecourse, 2);

            % Initialise model time course
            modelTimecourse = zeros(sum(countTimePoints), nRegress);

            % concatenate over all subjects for a session

            % Concatenate the model time courses
            begin = 1;
            % loop over numnber of subjects
            for numSub = 1:numOfSub
                % loop over number of regressors
                for numRegress = 1:nRegress
                    modelTimecourse(begin : numSub*size(timecourse, 1), numRegress) = ...
                        detrend(timecourse(:, numRegress), 0);
                end
                % end loop over number of regressors
                begin = begin + size(timecourse, 1); % update the number
            end
            % end loop over number of subjects
            clear timecourse;
            num_DataSets = numOfSub; % return the number of data sets

        else

            % sessions for a particular subject will be pulled from
            % SPM.mat

            % subject string
            for ii = 1:numOfSess
                sub_String(ii).name = ['Session ', num2str(ii)];
            end

            % Reference Function Information
            refInfo.SPMFile.name = spmMatrices.name;
            % replicate spm matrices over sessions
            spmMatrices = repmat(spmMatrices, 1, numOfSess);

            % option for entering through a text file
            if strcmpi(analType, 'batch')

                refAppData = icatb_pullRegressors_file('txtfile', deblank(refInfo.sortingTextFile), 'numOfSub', ...
                    numOfSub, 'numOfSess', numOfSess, 'keywd', keywd, 'spmmatflag', ...
                    stateDesignMatrix, 'spmNames', str2mat(spmMatrices.name), 'countTimePoints', countTimePoints, ...
                    'typeSelection', typeStr, 'data_sessionNumber', [], ...
                    'function_to_use', 'icatb_selreffunc', 'subjectNumber', [], ...
                    'all_time_points', allTimePoints, 'numSubjects', 1, 'numSessions', numOfSess);

            else

                % load the reference functions
                refAppData = icatb_selRefFunc('spmNames', str2mat(spmMatrices.name), 'numSubjects', 1, ...
                    'numSessions', numOfSess, 'substring', sub_String, 'typeSelection', typeStr, ...
                    'store_model_timecourse', 'yes', 'count_time_points', countTimePoints, ...
                    'all_time_points', allTimePoints, 'numOfSess',  numOfSess, 'numOfSub',  numOfSub);

            end


            % Loop over number of subjects
            for jj = 1:numOfSess
                modelIndex(jj, :) = refAppData(jj).value;
            end

            % Reference Function Information
            refInfo.modelIndex = modelIndex;

            clear modelIndex;

            % number of regressors
            nRegress = size(refInfo.modelIndex, 2);

            % selected regressor information
            % loop over number of sessions
            for ii = 1:numOfSess
                % Reference Function Information
                refInfo.selectedRegressors(ii).name = str2mat(refAppData(ii).name);
                refInfo.spmAppData(ii).data = struct('selectedRegressors', refAppData(ii).spmData.selectedRegressors, 'SPM', refAppData(ii).spmData.SPM, ...
                    'flag_refFunc', refAppData(ii).spmData.flag_refFunc);
            end
            % end loop over number of sessions

            % Initialise model time course
            modelTimecourse = zeros(sum(countTimePoints), nRegress);

            % a particular subject sessions are selected
            if strcmp(keywd, 'single_subject')
                begin = 1;
                % loop over sessions
                for numSess = 1:numOfSess
                    % time course of a particular session
                    timecourse = refAppData(numSess).spmData.timecourse;
                    % loop over regressors
                    for numRegress = 1:nRegress
                        modelTimecourse(begin : begin + size(timecourse, 1) - 1, numRegress) = ...
                            detrend(timecourse(:, numRegress), 0);
                    end
                    % end loop over regressors
                    begin = begin + size(timecourse, 1); % update number
                    clear timecourse;
                end
                % end loop over number of sessions
                clear refAppData;
                num_DataSets = numOfSess; % return the number of data sets
            end
            % end for concatentating subject sessions case

            % when all data sets are concatenated
            if strcmp(keywd, 'all_subjects_sessions')
                % Concatenate the model time courses
                beginSub = 1; dataCount = 0;
                % loop over subjects
                for numSub = 1:numOfSub
                    begin = 1;
                    % loop over sessions
                    for ses = 1:numOfSess
                        dataCount = dataCount + 1;
                        % time course of a particular session
                        timecourse = refAppData(ses).spmData.timecourse;
                        % loop over regressors
                        for numRegress = 1:nRegress
                            sessionTimecourse(begin : begin + size(timecourse, 1) - 1, numRegress) =  ...
                                detrend(timecourse(:, numRegress), 0);
                        end
                        % end loop over regressors
                        begin = begin + size(timecourse, 1);
                        clear timecourse;
                    end
                    % end loop over sessions

                    modelTimecourse(beginSub:beginSub + size(sessionTimecourse, 1) - 1, 1:nRegress) = sessionTimecourse;

                    beginSub = beginSub + size(sessionTimecourse, 1);
                    clear sessionTimecourse;
                end
                % end for loop over subjects
                clear refAppData;
                num_DataSets = numOfSub*numOfSess; % return the data sets
            end
            % end for concatenating all data-sets

        end
        % end for checking the keywd


    case 'diff_sub_diff_sess'

        % condition if any subject is selected
        if strcmp(keywd, 'single_subject_session')
            selectStr = str2mat(spmMatrices.name);

            % selected design index
            if ~exist('selectedDesignIndex', 'var')
                selectedDesignIndex = 0;
            end

            if selectedDesignIndex == 0
                % subject design matrix
                for ii = 1:numOfSub
                    designmatrixStr(ii).str = selectStr(numOfSess*(ii - 1) + 1, :);
                end

                % matrix index
                [matrixIndex] = icatb_listdlg('PromptString', 'Select SPM2 Design Matrix', 'SelectionMode','single',...
                    'ListString', str2mat(designmatrixStr.str), 'movegui', 'east', 'windowStyle', 'modal');
                matrixIndex = matrixIndex*numOfSess;
            else
                matrixIndex = selectedDesignIndex;
            end

            % pull time points for a subject
            tempTimePoints = zeros(1, numOfSess);
            subIndex = ceil(matrixIndex/numOfSess);
            for nn = 1:numOfSess
                tempTimePoints(nn) = allTimePoints.sub(subIndex).sess(nn);
            end

            if sessionIndex == 0
                sessionIndex = [];
            end

            if isempty(matrixIndex)
                error('Design matrix is not selected');
            end

            if isempty(sessionIndex)
                tempTimePoints = tempTimePoints(1);
            end


            %             if ~isempty(sessionIndex)

            % check the regressor information
            if strcmpi(analType, 'batch')

                [spmData] = icatb_pullRegressors_file('txtfile', deblank(refInfo.sortingTextFile), 'numOfSub', numOfSub, ...
                    'numOfSess', numOfSess, 'keywd', keywd, 'spmmatflag', stateDesignMatrix, ...
                    'spmNames', spmMatrices(matrixIndex).name, 'countTimePoints', tempTimePoints, 'typeSelection', typeStr, ...
                    'data_sessionNumber', sessionIndex, 'function_to_use', 'icatb_loadspm_new', 'subjectNumber', subIndex);


            else

                % check the SPM design matrix
                [spmData] = icatb_loadSPM_new('spmName', spmMatrices(matrixIndex).name, 'countTimePoints', ...
                    tempTimePoints, 'typeStr', typeStr, 'data_sessionNumber', sessionIndex, 'check_spm_design_matrix', ...
                    'yes', 'check_regressors', 'yes', 'flag_selecting_regressors', 'yes', 'get_timecourses', 'yes'); % model Index doesn't contain actual indices

            end

            %             else
            %                 % load all regressors
            %                 % check the SPM design matrix
            %                 % check the SPM design matrix
            %                 [spmData] = icatb_loadSPM_new('spmName', spmMatrices(matrixIndex).name, 'countTimePoints', ...
            %                     tempTimePoints(1), 'typeStr', typeStr, 'data_sessionNumber', 1, 'check_spm_design_matrix', ...
            %                     'yes', 'check_regressors', 'yes', 'flag_selecting_regressors', 'yes', 'get_timecourses', 'yes'); % model Index doesn't contain actual indices
            %             end

            % Reference Function Information
            refInfo.SPMFile.name = spmMatrices(matrixIndex).name;
            refInfo.modelIndex = spmData.getIndex;
            % selected regressors information
            refInfo.selectedRegressors.name = str2mat(spmData.selectedRegressors);
            refInfo.spmAppData.data = spmData;
            timecourse = spmData.timecourse; % get the time course

            % number of regressors
            nRegress = size(timecourse, 2);

            modelTimecourse = detrend(timecourse, 0);

            num_DataSets = 1;

            % sessions for a particular subject
        elseif strcmp(keywd, 'single_subject')

            subjectNumber = selectedDesignIndex;

            designIndices = (subjectNumber - 1)*numOfSess + 1 : subjectNumber*numOfSess;

            for ii = 1:numOfSess
                spmNames(ii).name = spmMatrices(designIndices(ii)).name;
                refInfo.SPMFile(ii).name = spmNames(ii).name;
            end

            % selected design matrices
            selectStr = str2mat(spmNames.name);

            % subject string
            for ii = 1:numOfSess
                subString(ii).name = ['Session ', num2str(ii)];
            end

            % check the regressor information
            if strcmpi(analType, 'batch')


                % If spm2 model is selected for each session then use this
                % function
                %refAppData = icatb_selRefFunc(selectStr, typeStr, 1,
                %numOfSess, subString);
                % reference application data
                refAppData = icatb_pullRegressors_file('txtfile', deblank(refInfo.sortingTextFile), 'numOfSub', ...
                    numOfSub, 'numOfSess', numOfSess, 'keywd', keywd, 'spmmatflag', ...
                    stateDesignMatrix, 'spmNames', selectStr, 'countTimePoints', countTimePoints, ...
                    'typeSelection', typeStr, 'data_sessionNumber', [], ...
                    'function_to_use', 'icatb_selreffunc', 'subjectNumber', subjectNumber, ...
                    'all_time_points', allTimePoints, 'numSubjects', 1, 'numSessions', numOfSess);


            else

                refAppData = icatb_selRefFunc('spmNames', selectStr, 'numSubjects', 1, ...
                    'numSessions', numOfSess, 'substring', subString, 'typeSelection', typeStr, ...
                    'store_model_timecourse', 'yes', 'count_time_points', countTimePoints, 'all_time_points', ...
                    allTimePoints, 'numOfSess',  numOfSess, 'numOfSub',  numOfSub);

            end
            % end for checking the regressors information

            % Loop over number of subjects
            for jj = 1:numOfSess
                modelIndex(jj, :) = refAppData(jj).value;
            end

            % Reference Function Information
            refInfo.modelIndex = modelIndex;

            clear modelIndex;

            % number of regressors
            nRegress = size(refInfo.modelIndex, 2);

            % selected regressor information
            for ii = 1:numOfSess
                % Reference Function Information
                refInfo.selectedRegressors(ii).name = str2mat(refAppData(ii).name);
                refInfo.spmAppData(ii).data = struct('selectedRegressors', refAppData(ii).spmData.selectedRegressors, 'SPM', refAppData(ii).spmData.SPM, ...
                    'flag_refFunc', refAppData(ii).spmData.flag_refFunc);
                %                refInfo.spmAppData(ii).data.selectedRegressors =
                %                refAppData(ii).spmData.selectedRegressors;
            end

            begin = 1;
            % Initialise model time course
            modelTimecourse = zeros(sum(countTimePoints), nRegress);
            for numSess = 1:numOfSess
                % time course of a particular session
                timecourse = refAppData(numSess).spmData.timecourse;
                for numRegress = 1:nRegress
                    modelTimecourse(begin : begin + size(timecourse, 1) - 1, numRegress) = detrend(timecourse(:, numRegress), 0);
                end
                begin = begin + size(timecourse, 1);
                clear timecourse;
            end

            num_DataSets = numOfSess;

            clear refAppData;

            % condition if more than one subject is selected
        else

            if strcmp(keywd, 'single_session')

                sessionNumber = selectedDesignIndex;

                for jj = 1:numOfSub
                    subString(jj).name = ['Subject ', num2str(jj)];
                end
                % If spm2 model is selected for each function then use this
                % function

                for ii = 1:numOfSub
                    indexDesignMat = numOfSess*(ii - 1) + sessionNumber;
                    spmNames(ii).name = spmMatrices(indexDesignMat).name;
                end


                % read regressors from file
                if strcmpi(analType, 'batch')
                    % reference application data
                    refAppData = icatb_pullRegressors_file('txtfile', deblank(refInfo.sortingTextFile), 'numOfSub', ...
                        numOfSub, 'numOfSess', numOfSess, 'keywd', keywd, 'spmmatflag', ...
                        stateDesignMatrix, 'spmNames', str2mat(spmNames.name), 'countTimePoints', countTimePoints, ...
                        'typeSelection', typeStr, 'data_sessionNumber', sessionNumber, ...
                        'function_to_use', 'icatb_selreffunc', 'subjectNumber', [], ...
                        'all_time_points', allTimePoints, 'numSubjects', numOfSub, 'numSessions', 1);

                else

                    refAppData = icatb_selRefFunc('spmNames', str2mat(spmNames.name), 'numSubjects', numOfSub, ...
                        'numSessions', 1, 'substring', subString, 'typeSelection', typeStr, ...
                        'store_model_timecourse', 'yes', 'count_time_points', countTimePoints,  'data_sessionNumber', ...
                        sessionNumber, 'all_time_points', allTimePoints, 'numOfSess',  numOfSess, ...
                        'numOfSub',  numOfSub);

                end

                % Loop over number of subjects
                for jj = 1:numOfSub
                    modelIndex(jj, :) = refAppData(jj).value;
                end

                % Reference Function Information
                refInfo.modelIndex = modelIndex;

                % number of regressors
                nRegress = size(refInfo.modelIndex, 2);

                clear modelIndex;

                % selected regressor information
                for ii = 1:numOfSub
                    % Reference Function Information
                    refInfo.selectedRegressors(ii).name = str2mat(refAppData(ii).name);
                    refInfo.spmAppData(ii).data =  struct('selectedRegressors', refAppData(ii).spmData.selectedRegressors, 'SPM', refAppData(ii).spmData.SPM, ...
                        'flag_refFunc', refAppData(ii).spmData.flag_refFunc);
                end

                % Initialise model time course
                modelTimecourse = zeros(sum(countTimePoints), nRegress);

                % Concatenate the model time courses
                begin = 1;
                for numSub = 1:numOfSub
                    % Reference Function Information
                    refInfo.SPMFile(numSub).name = spmMatrices(numSub).name;
                    timecourse = refAppData(numSub).spmData.timecourse;
                    for numRegress = 1:nRegress
                        modelTimecourse(begin : begin + size(timecourse, 1) - 1, numRegress) = detrend(timecourse(:, numRegress), 0);
                    end
                    begin = begin + size(timecourse, 1);
                    clear timecourse;
                end

                num_DataSets = numOfSub;

                clear refAppData;

            elseif strcmp(keywd, 'all_subjects_sessions')

                % If spm2 model is selected for each function then use this
                % function
                dataCount = 0;
                for ii = 1:numOfSub
                    for jj = 1:numOfSess
                        dataCount = dataCount + 1;
                        subString(dataCount).name = ['Subject ', num2str(ii), ' Session ', num2str(jj)];
                    end
                end

                disp('Loading Reference Function(s) From SPM Matrices');

                % read regressors from file
                if strcmpi(analType, 'batch')
                    % reference application data
                    refAppData = icatb_pullRegressors_file('txtfile', deblank(refInfo.sortingTextFile), 'numOfSub', ...
                        numOfSub, 'numOfSess', numOfSess, 'keywd', keywd, 'spmmatflag', ...
                        stateDesignMatrix, 'spmNames', str2mat(spmMatrices.name), 'countTimePoints', countTimePoints, ...
                        'typeSelection', typeStr, 'data_sessionNumber', [], ...
                        'function_to_use', 'icatb_selreffunc', 'subjectNumber', [], ...
                        'all_time_points', allTimePoints, 'numSubjects', numOfSub, 'numSessions', numOfSess);

                else

                    % reference function for each dataset
                    refAppData = icatb_selRefFunc('spmNames', str2mat(spmMatrices.name), 'numSubjects', numOfSub, ...
                        'numSessions', numOfSess, 'substring', subString, 'typeSelection', typeStr, ...
                        'store_model_timecourse', 'yes', 'count_time_points', countTimePoints, 'all_time_points', allTimePoints,  'numOfSess',  ...
                        numOfSess, 'numOfSub',  numOfSub);

                end
                % end for if

                % Loop over number of subjects
                for jj = 1:numOfSub*numOfSess
                    modelIndex(jj, :) = refAppData(jj).value;
                end

                % Reference Function Information
                refInfo.modelIndex = modelIndex;

                nRegress = size(refInfo.modelIndex, 2);

                clear modelIndex;

                % selected regressor information
                for ii = 1:numOfSub*numOfSess
                    % Reference Function Information
                    refInfo.selectedRegressors(ii).name = str2mat(refAppData(ii).name);
                    %refInfo.spmAppData(ii).data = refAppData(ii).spmData;
                    refInfo.spmAppData(ii).data = struct('selectedRegressors', refAppData(ii).spmData.selectedRegressors, 'SPM', refAppData(ii).spmData.SPM, ...
                        'flag_refFunc', refAppData(ii).spmData.flag_refFunc);
                end


                % Initialise model time course
                modelTimecourse = zeros(sum(countTimePoints), nRegress);

                % Concatenate the model time courses
                beginSub = 1; dataCount = 0;
                for numSub = 1:numOfSub
                    begin = 1;
                    for ses = 1:numOfSess
                        dataCount = dataCount + 1;
                        % Reference Function Information
                        refInfo.SPMFile(dataCount).name = spmMatrices(dataCount).name;
                        timecourse = refAppData(dataCount).spmData.timecourse;
                        for numRegress = 1:nRegress
                            sessionTimecourse(begin : begin + size(timecourse, 1) - 1, numRegress) =  detrend(timecourse(:, numRegress), 0);
                        end
                        begin = begin + size(timecourse, 1);
                        clear timecourse;
                    end

                    modelTimecourse(beginSub:beginSub + size(sessionTimecourse, 1) - 1, 1:nRegress) = sessionTimecourse;

                    beginSub = beginSub + size(sessionTimecourse, 1);
                    clear sessionTimecourse;
                end
                clear refAppData; % clear the reference application data
                num_DataSets = numOfSub*numOfSess;
            end


        end
        % end for taking into account all the types for the sorting
        % criteria correlation
end

