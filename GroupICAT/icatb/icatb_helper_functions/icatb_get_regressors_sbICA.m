function sesInfo = icatb_get_regressors_sbICA(inputFile, sesInfo)
% Get the reference functions or regressor names for Semi-blind ICA
% algorithm
%
% Inputs:
% sesInfo - structure containing parameter information
%
% Output:
% sesInfo - Updated structure containing parameter information

if ~isempty(inputFile)

    % pulling the time courses from the semi-blind ICA
    if sesInfo.userInput.algorithm == 9

        % check total datasets
        if sesInfo.userInput.numOfSub*sesInfo.userInput.numOfSess == 1
            totalFiles = sesInfo.userInput.files(1).name;
            % get the count for files
            [diffTimePoints] = sesInfo.userInput.diffTimePoints;

            % check if the design matrix is specified for semi-blind ICA
            % check whether the index is passed or not,
            if isfield(sesInfo.userInput, 'designMatrix')
                % open the design matrix file
                checkFile = sesInfo.userInput.designMatrix.name;
                if isempty(checkFile)
                    error('Design matrix should be specified to run semi-blind ICA');
                end

                % do error checking for SPM design matrix
                [spmData] = icatb_loadSPM_new('spmName', checkFile, 'countTimePoints', ...
                    diffTimePoints(1), 'data_sessionNumber', 1, ...
                    'check_spm_design_matrix', 'yes', 'check_regressors', 'yes');
                allNames = spmData.sel_refnames;
                allNames = deblank(allNames);
                keywd = 'refFunNames';
                [inputData, status, message] = icatb_read_variables(inputFile, keywd, 'vector', 'cell');
                % read Reference function names here
                refFunNames = getfield(inputData, keywd);
                if length(refFunNames) > 2
                    disp('-- By default using only two reference function names');
                    for jj = 1:2
                        temp{jj} = refFunNames{jj};
                    end
                    clear refFunNames;
                    refFunNames = temp;
                    clear temp;
                end

                % check the index
                for ii = 1:length(refFunNames)
                    tempIndex = strmatch(lower(refFunNames{ii}), lower(allNames), 'exact');
                    if isempty(tempIndex)
                        error(['Specified reference function with name ', lower(refFunNames{ii}), ' doesn''t exist']);
                    end
                    % store the indices
                    getIndex(ii) = tempIndex;
                end

                refFunNames = str2mat(refFunNames);

                spmData.selectedRegressors = refFunNames; % store the selected regressors
                spmData.getIndex = getIndex; % store the selected information in spm data

                % get the selected time course information
                [spmData] = icatb_loadSPM_new('spmData', spmData, 'get_timecourses', 'yes', ...
                    'check_spm_design_matrix', 'no', 'check_regressors', 'no', 'flag_selecting_regressors', 'no');

                timecourse = spmData.timecourse; % get the selected time courses information

                clear spmData;

                numICAOptions = length(sesInfo.userInput.ICA_Options);
                % ICA Options
                ICAOptions = sesInfo.userInput.ICA_Options;
                % Append the ICA Options
                ICAOptions{numICAOptions + 1} = 'TC';
                ICAOptions{numICAOptions + 2} = timecourse;
                % Pass the updated ICA Options and append the structure
                % with all the time courses that are present
                sesInfo.userInput.ICA_Options = ICAOptions;
            else
                error('The field designMatrix doesn''t exist for sesInfo');
            end
            % end for checking design matrix

        else

            error('Presently Semi-blind ICA works with only one data set');

        end
        % end for checking total number of data sets

    end
    % check semi-blind ICA

end