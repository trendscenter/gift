function [regressSubSess] = icatb_readRegressors_file(varargin)
% read regressors from the text file depending upon the key word for the
% data sets and the spm mat flag
%
% Input:
% 1. txtFile - text file for getting the regressor information
% 2. numOfSub - Number of subjects.
% 3. numOfSess - Number of sessions.
% 4. keywd - Keyword for getting the regressors
% 5. spmMatFlag - SPM mat flag
%
% Output:
% regress_subjects_sessions structure containing the regressor names in a
% structure like regress_subjects_sessions.sub(1).sess(1).names
%
% loop over number of arguments
for ii = 1:2:nargin
    if strcmpi(varargin{ii}, 'txtfile')
        % text file
        txtFile = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'numofsub')
        % number of subjects
        numOfSub = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'numofsess')
        % number of sessions
        numOfSess = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'spmmatflag')
        % spm mat flag
        spmMatFlag = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'subjectnumber')
        % subject number
        subjectNumber = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'sessionnumber')
        % session number
        sessionNumber = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'keywd')
        % keywd for data sets
        keywd = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'typeSelection')
        % selection type
        typeSelection = varargin{ii + 1};
    end
end
% end loop over number of arguments

% if session index is empty
if isempty(sessionNumber)
    sessionNumber = 1;
end

% if subject index is empty
if isempty(subjectNumber)
    subjectNumber = 1;
end

% selection type
if ~exist('typeSelection', 'var')
    typeSelection = 'multiple';
end

% check the session number
if sessionNumber > numOfSess
    error('Session number cannot be greater than the number of sessions');
end

% check the session number
if subjectNumber > numOfSub
    error('Subject number cannot be greater than the number of subjects');
end

% initialise variable for storing the regressors information
regressSubSess.sub.sess.names = [];
% read the remaining statements

% steps:
% 1. Read regressor names and close the file.
% 2. open the icatb_selectModel_Criteria file
% 3. Return the indices of the regressors when comparing (if none print
% error).
% 4. Set the list box to the corresponding values in case for the
% icatb_selRefFunc file.


% compare the SPM mat flag
if strcmpi(spmMatFlag, 'same_sub_same_sess')
    % regressors for all subjects and sessions
    regressKeywd = 'regress_all';
    % store the keywd, dataType and vector in keywdsNeeded
    %%keywdsNeeded.tag = regressKeywd;  keywdsNeeded.flag = 'vector'; keywdsNeeded.dataType = 'string';
    % get the reference function names
    %regressSubSess.sub.sess.names = icatb_read_vector_file(fid, regressKeywd, 'string');
    inputData = icatb_read_variables(txtFile, regressKeywd, 'vector', 'string', 'no_display_error');
    regressSubSess.sub.sess.names = getfield(inputData, regressKeywd);
    clear inputData;

elseif strcmpi(spmMatFlag, 'same_sub_diff_sess')
    % regressors for sessions but same over subjects

    % consider cases for different data sets
    switch lower(keywd)

        case 'single_subject_session'

            regressKeywd = ['regress_sub1_sess', num2str(sessionNumber)];
            % get the reference function names
            inputData = icatb_read_variables(txtFile, regressKeywd, 'vector', 'string', 'no_display_error');
            regressSubSess.sub.sess.names = getfield(inputData, regressKeywd);
            clear inputData;

        case 'single_session'

            regressKeywd = ['regress_sub1_sess', num2str(sessionNumber)];
            % get the reference function names
            inputData = icatb_read_variables(txtFile, regressKeywd, 'vector', 'string', 'no_display_error');
            regressSubSess.sub.sess.names = getfield(inputData, regressKeywd);
            clear inputData;

        otherwise

            % loop over sessions
            for ii = 1:numOfSess
                % keywd for entering sessions
                regressKeywd =  ['regress_sub1_sess', num2str(ii)];
                % get the reference function names
                inputData = icatb_read_variables(txtFile, regressKeywd, 'vector', 'string', 'no_display_error');
                regressSubSess.sub(1).sess(ii).names = getfield(inputData, regressKeywd);
                clear inputData;

                if ii == 1
                    previousLength = length(regressSubSess.sub(1).sess(1).names);
                else
                    currentLength = length(regressSubSess.sub(1).sess(ii).names);
                    % check if equal number of regressors are not entered
                    if previousLength ~= currentLength
                        error(['Number of regressors should be equal over sessions. ', ...
                            ' Enter correct number of regressors for ', regressKeywd]);
                    end
                end

            end
            % end for loop over sessions

    end
    % end for cases for different data sets

else
    % regressors for subjects and sessions

    % consider cases for different data sets
    switch lower(keywd)

        case 'single_subject_session'

            % keywd for regressors
            regressKeywd = ['regress_sub', num2str(subjectNumber), '_sess', num2str(sessionNumber)];
            % get the reference function names
            inputData = icatb_read_variables(txtFile, regressKeywd, 'vector', 'string', 'no_display_error');
            regressSubSess.sub.sess.names = getfield(inputData, regressKeywd);
            clear inputData;

        case 'single_session'

            % loop over subjects
            for ii = 1:numOfSub
                % keyword for regressors
                regressKeywd = ['regress_sub', num2str(ii), '_sess', num2str(sessionNumber)];
                % get the reference function names
                inputData = icatb_read_variables(txtFile, regressKeywd, 'vector', 'string', 'no_display_error');
                regressSubSess.sub(ii).sess(1).names = getfield(inputData, regressKeywd);
                clear inputData;

                if ii == 1
                    previousLength = length(regressSubSess.sub(1).sess(1).names);
                else
                    currentLength = length(regressSubSess.sub(ii).sess(1).names);
                    % check if equal number of regressors are not entered
                    if previousLength ~= currentLength
                        error(['Number of regressors should be equal over subjects.',
                            ' Enter correct number of regressors for ', regressKeywd]);
                    end
                end
            end
            % end loop over subjects

        case 'single_subject'

            % loop over sessions
            for ii = 1:numOfSess
                % keyword for regressors
                regressKeywd = ['regress_sub', num2str(subjectNumber), '_sess', num2str(ii)];
                % get the reference function names
                inputData = icatb_read_variables(txtFile, regressKeywd, 'vector', 'string', 'no_display_error');
                regressSubSess.sub(1).sess(ii).names = getfield(inputData, regressKeywd);
                clear inputData;

                if ii == 1
                    previousLength = length(regressSubSess.sub(1).sess(1).names);
                else
                    currentLength = length(regressSubSess.sub(1).sess(ii).names);
                    % check if equal number of regressors are not entered
                    if previousLength ~= currentLength
                        error(['Number of regressors should be equal over sessions.', ...
                            ' Enter correct number of regressors for ', regressKeywd]);
                    end
                end
            end
            % end loop over subjects

        otherwise

            count = 0;
            % loop over subjects
            for ii = 1:numOfSub
                % loop over sessions
                for jj = 1:numOfSess
                    count = count + 1;
                    % keywd for entering sessions
                    regressKeywd =  ['regress_sub', num2str(ii), '_sess', num2str(jj)];
                    % get the reference function names
                    inputData = icatb_read_variables(txtFile, regressKeywd, 'vector', 'string', 'no_display_error');
                    regressSubSess.sub(ii).sess(jj).names = getfield(inputData, regressKeywd);
                    clear inputData;
                    if count == 1
                        previousLength = length(regressSubSess.sub(1).sess(1).names);
                    else
                        currentLength = length(regressSubSess.sub(ii).sess(jj).names);
                        % check if equal number of regressors are not entered
                        if previousLength ~= currentLength
                            error(['Number of regressors should be equal over data sets.', ...
                                ' Enter correct number of regressors for ', regressKeywd]);
                        end
                    end

                end
                % end for loop over sessions

            end
            % end for loop overs subjects

    end
    % end for cases for different data sets

end
% end for comparing the spm mat flag

% check the selection type
if ~strcmpi(typeSelection, 'multiple')
    % get only one regressor for sorting by correlation
    disp('Only one regressor selection is allowed. By default selecting the first regresor.');
    % loop over subjects
    for ii = 1:length(regressSubSess.sub)
        % loop over sessions
        for jj = 1:length(regressSubSess.sub(ii).sess)
            % get the current regressor names
            currentRegressors = regressSubSess.sub(ii).sess(jj).names;
            % get the first regressor name
            regressSubSess.sub(ii).sess(jj).names = {currentRegressors{1}};
            clear currentRegressors;
        end
        % end loop over sessions
    end
    % end loop over subjects
end
% end for checking the selection type