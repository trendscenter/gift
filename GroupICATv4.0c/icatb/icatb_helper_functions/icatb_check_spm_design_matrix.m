function icatb_check_spm_design_matrix(designMatrix, numOfSub, numOfSess, diffTimePoints, spmMatFlag, inputFile)
% Check if the design matrix matches the data
%
% Inputs:
% 1. designMatrix - designMatrix structure
% 2. numOfSub - Number of subjects
% 3. numOfSess - Number of sessions
% 4. diffTimePoints - time points vector
% 5. inputFile - input file
%


if ~isempty(designMatrix(1).name)
    % consider the three cases here

    % check the design matrix flag variable
    if strcmpi(spmMatFlag, 'same_sub_same_sess')
        % check the number of time points
        checkTimePoints = find(diffTimePoints == diffTimePoints(1));
        if length(checkTimePoints) ~= length(diffTimePoints)
            error(['Please select same number of scans or time points if you want ' ...
                'to select same set of regressors for all subjects. Check other options for keyword_designMatrix' ...
                ' variable in file: ', inputFile]);
        end
        % do spm checking
        [spmData] = icatb_loadSPM_new('spmName', deblank(designMatrix(1).name), 'countTimePoints', ...
            diffTimePoints(1), 'check_spm_design_matrix', 'yes', 'data_sessionNumber', 1);

    elseif strcmpi(spmMatFlag, 'same_sub_diff_sess')
        % check the number of time points
        startTp = 1; endTp = 0; tpOverSessions = diffTimePoints(1:numOfSess);
        % loop over subjects
        for nSub = 1:numOfSub
            % end point
            endTp = nSub*numOfSess;
            checkTimePoints = diffTimePoints(startTp:endTp);
            checkTimePoints = (tpOverSessions == checkTimePoints);
            checkTp = find(checkTimePoints == 1);
            if length(checkTp) ~= numOfSess
                error('Error:spmCheck', ['Subject %s images are [%s].\nThe number of images cannot be different over subjects.\nUse ', ...
                    ' diff_sub_diff_sess value for keyword_designMatrix variable in file:\n %s if you want to run the analysis with different design matrices.'], ...
                    num2str(nSub), num2str(diffTimePoints(startTp:endTp)), inputFile);
            end
            startTp = endTp + 1; % increment the starting point
        end
        % end for loop over subjects

        % loop over sessions
        for nSess = 1:numOfSess
            % check spm design matrix
            [spmData] = icatb_loadSPM_new('spmName', deblank(designMatrix(1).name), 'countTimePoints', ...
                diffTimePoints(1:numOfSess), 'check_spm_design_matrix', 'yes', 'data_sessionNumber', nSess);
        end
        % end for loop over sessions

    else
        % loop over number of subjects
        for ii = 1:numOfSub
            startTp = (ii - 1)*numOfSess + 1; % start time point
            endTp = ii*numOfSess; % end time point
            % loop over sessions
            for jj = 1:numOfSess
                [spmData] = icatb_loadSPM_new('spmName', deblank(designMatrix(ii).name), ...
                    'countTimePoints', diffTimePoints(startTp:endTp), 'data_sessionNumber', numOfSess, ...
                    'check_spm_design_matrix', 'yes');
            end
            % end for loop over sessions
        end
        % end for loop over subjects

    end
    % end for checking the design matrix flag variable
end
% end for checking the empty state of the design matrix
