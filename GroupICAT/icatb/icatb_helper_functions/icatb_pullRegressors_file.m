function [refAppData] = icatb_pullRegressors_file(varargin)
% function returns reference application data variable refAppData that
% contains the necessary variable for returning the regressor time courses
%
% Inputs must be in pairs:
% 1. 'txtfile' - text file
% 2. 'numOfSub' - Number of subjects
% 3. 'numOfSess' - Number of sessions
% 4. 'keywd' - keyword for data sets
% 5. 'stateDesignMatrix' - State of design matrix
% 6. 'spmNames' - spmNames
% 7. 'countTimePoints' - Count for time points
% 8. 'typeSelection' - Single selection or multiple selection
% 9. 'data_sessionNumber' - Session number.
% 10. 'all_time_points' -
% 11. 'numSubjects' -
% 12. 'numSessions' -
% 13. 'subString' -
% 14. 'function_to_use' - icatb_loadSPM_new or icatb_selRefFunc
% 15. subjectNumber
% Output: refAppData


icatb_defaults;
global METHOD_ENTERING_REGRESSORS;

analType = METHOD_ENTERING_REGRESSORS;

% check if the regressors are automatically selected
if strcmpi(analType, 'automatic')
    flagUseRegressors = 'all';
else
    flagUseRegressors = 'user-specified';
end

subjectNumber = [];

% loop over arguments
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
    elseif strcmpi(varargin{ii}, 'keywd')
        % keyword
        keywd = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'spmmatflag')
        % state for design matrix (Same regressors over data sets,
        % Different regressors over sessions, Different regressors over
        % sessions and subjects)
        spmMatFlag = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'spmnames')
        % spm names
        spmNames = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'counttimepoints')
        % count for time points
        countTimePoints = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'typeselection')
        % selection type
        typeSelection = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'data_sessionnumber')
        % data session number
        data_sessionNumber = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'numsubjects')
        % number of subjects involved in sorting
        numSubjects = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'numsessions')
        % number of sessions involved in sorting
        numSessions = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'substring')
        % subject string
        subjectString = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'function_to_use')
        % function to use (icatb_loadSPM_new or icatb_selRefFunc)
        function_to_use = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'subjectnumber')
        % subject number
        subjectNumber = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'all_time_points')
        % all time points (structure containing time points information for
        % all data sets)
        allTimePoints = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'flaguseregressors')
        flagUseRegressors = varargin{ii + 1};
    end
end
% end for loop over arguments


% read regressor information from the text file
% [regress_subjects_sessions] = icatb_readRegressors_file('txtfile', 'Input_data_multiple_regression.txt', ...
% 'numOfSub', 4, 'numOfSess', 3, 'spmMatFlag', 'diff_sub_diff_sess', 'subjectNumber', 3, ...
% 'sessionNumber', 1, 'keywd', 'single_subject');

% function to use
if strcmpi(function_to_use, 'icatb_loadspm_new')

    % read regressors from file
    if ~strcmpi(flagUseRegressors, 'all')


        % read regressors from subjects and sessions
        [regress_subjects_sessions] = icatb_readRegressors_file('txtfile', deblank(txtFile), 'numOfSub', numOfSub, ...
            'numOfSess', numOfSess, 'spmmatflag', spmMatFlag, 'keywd', keywd, 'sessionNumber', ...
            data_sessionNumber, 'subjectNumber', subjectNumber, 'typeSelection', typeSelection);

    end

    % load model time courses and also check the SPM design matrix
    [refAppData] = icatb_loadSPM_new('spmName', spmNames, 'countTimePoints', ...
        countTimePoints, 'typeStr', typeSelection, 'data_sessionNumber', data_sessionNumber, ...
        'check_spm_design_matrix', 'yes', 'check_regressors', 'yes');

    allNames = refAppData.sel_refnames;
    allNames = deblank(allNames);

    % use regr
    if ~strcmpi(flagUseRegressors, 'all')
        % reference function names from the text file
        refFunNames = regress_subjects_sessions.sub.sess.names;

        % check the index
        for ii = 1:length(refFunNames)
            tempIndex = [];
            for jj = 1:size(allNames, 1)
                % remove the trailing parts
                currentRefName = deblank(allNames(jj, :));
                %tempIndex = strmatch(lower(refFunNames), lower(allNames), 'exact');
                if strcmpi(currentRefName, refFunNames{ii})
                    %error(['Specified reference function with name ', lower(refFunNames{ii}), ' doesn''t exist']);
                    tempIndex = jj;
                end
            end
            % end for loop for allnames
            if isempty(tempIndex)
                error(['Specified reference function with name ', lower(refFunNames{ii}), ' doesn''t exist']);
            end

            % store the indices
            getIndex(ii) = tempIndex;
        end

        % reference function names
        refFunNames = str2mat(refFunNames);
    else
        % automatically select regressors from SPM.mat
        % single file selection
        if strcmpi(typeSelection, 'multiple')
            % get all the regressor indices
            getIndex = [1:1:size(allNames, 1)];
            % reference function names
            refFunNames = allNames;
        else
            % get all the regressors
            getIndex = 1;
            % reference function names
            refFunNames = deblank(allNames(getIndex, :));
            fprintf('\n');
            disp(['Regressor selected for temporal sorting is ', refFunNames]);
        end
        % end for file selection

    end
    % end for checking

    refAppData.selectedRegressors = refFunNames; % store the selected regressors
    refAppData.getIndex = getIndex; % store the selected information in spm data

    % get the selected time course information
    [refAppData] = icatb_loadSPM_new('spmData', refAppData, 'get_timecourses', 'yes', ...
        'check_spm_design_matrix', 'no', 'check_regressors', 'no', 'flag_selecting_regressors', 'no');

else


    % selecting regressors function
    [refAppData] = icatb_selRefFunc('txtfile', deblank(txtFile), 'numOfSub', numOfSub, ...
        'numOfSess', numOfSess, 'spmmatflag', spmMatFlag, 'keywd', keywd, 'data_sessionnumber', ...
        data_sessionNumber, 'subjectNumber', subjectNumber, 'handle_visibility', 'off', 'numSubjects', ...
        numSubjects, 'numSessions', numSessions, 'spmNames', spmNames, 'typeSelection', typeSelection, ...
        'store_model_timecourse', 'yes' , 'count_time_points', countTimePoints, 'all_time_points', ...
        allTimePoints, 'flagUseRegressors', flagUseRegressors);


end
% end for checking