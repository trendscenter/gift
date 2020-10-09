function [spmData] = icatb_loadSPM_new(varargin)

% Procedure: There are four steps:
% First Step: Check SPM design matrix and do error checking.
% ************** Store the data in spmData ******************
% Store the following fields: nscan, flag_refFunc ('all' or
% 'session_specfic') and data_sessionNumber.
% Update SPM.Sess and store in spmData as field SPM.
%
% Second Step: Get session specific regressors if possible.
% ************** Store the data in spmData ******************
% store sel_refnames and sel_indices as fields in spmData
%
% Third Step: Select regressors
% get model time course information as modelT
% store modelT in spmData with the following fields:
% 1. designMatrixName
% 2. selectedRegressors
% 3. name (all regressors available)
% 4. nsess - modified length of scans
% 5. actualIndices W.r.t SPM.xX.name as sel_indices
% ************** Store the data in spmData ******************
% store modelT information
%
% Fourth Step: Getting time course information
% Use SPM and modelT fields in spmData to get the time courses and modified
% SPM structure
%
% Inputs:
% 1. spmName - full file name of the spm design matrix
% 2. countTimePoints - length of time points for each session and is a
% vector tor like [220 220 ...] depending on the number of sessions used.
% 3. optional argument: sessionNumber to pull out the reference function names
% time courses of that session.
%
% Output:
% spmData containing the information about selected time courses and
% regressors

% Initialise vars
spmName = [];
countTimePoints = [];
data_sessionNumber = []; % session number is used to pull the reference function names
flag_refFunc = 'all'; % by default pull all the reference functions names and timecourses
analType = 'gui'; % default is GUI
prompt = 'Select Timecourse'; % title of the listbox
typeStr = 'single'; % selection is single by default
maxSelection = []; % Maximum selection
check_spm_design_matrix = 'yes';
check_regressors = 'no';
flag_selecting_regressors = 'no';
get_timecourses = 'no';
title_fig = 'Reference Time courses';

% loop over number of input arguments
for ii = 1:2:nargin
    if strcmp(lower(varargin{ii}), 'spmname')
        spmName = varargin{ii + 1};
    elseif strcmp(lower(varargin{ii}), 'counttimepoints')
        countTimePoints = varargin{ii + 1};
    elseif strcmp(lower(varargin{ii}), 'data_sessionnumber')
        data_sessionNumber = varargin{ii + 1};
    elseif strcmp(lower(varargin{ii}), 'analtype')
        analType = varargin{ii + 1};
    elseif strcmp(lower(varargin{ii}), 'typestr')
        typeStr = varargin{ii + 1};
    elseif strcmp(lower(varargin{ii}), 'maxselection')
        maxSelection = varargin{ii + 1};
    elseif strcmp(lower(varargin{ii}), 'flag_selecting_regressors')
        flag_selecting_regressors = varargin{ii + 1};
    elseif strcmp(lower(varargin{ii}), 'check_spm_design_matrix')
        check_spm_design_matrix = lower(varargin{ii + 1});
    elseif strcmp(lower(varargin{ii}), 'check_regressors')
        check_regressors = varargin{ii + 1};
    elseif strcmp(lower(varargin{ii}), 'get_timecourses')
        get_timecourses = varargin{ii + 1};
    elseif strcmp(lower(varargin{ii}), 'spmdata')
        spmData = varargin{ii + 1};
    elseif strcmp(lower(varargin{ii}), 'prompt')
        prompt = varargin{ii + 1};
    elseif strcmp(lower(varargin{ii}), 'title_fig')
        title_fig = varargin{ii + 1};
    end
end

try

    if isempty(data_sessionNumber)
        data_sessionNumber = 1;
    end

    % check the design matrix
    if strcmp(check_spm_design_matrix, 'yes')

        % load the spm file
        % Note: storing in spmVar as loading SPM design matrix produces SPM variable
        % which is same as function name
        spmVar = load(spmName);

        % Check the existence of SPM variable
        if ~exist('spmVar', 'var')
            error('SPM  variable doesn''t exist. Please check the .MAT file.');
        end

        fieldNamesOfSPM = fieldnames(spmVar); % get the field names

        spmMat = getfield(spmVar, fieldNamesOfSPM{1}); % get the SPM design matrix here

        clear spmVar;

        % get the number of scans
        if isfield(spmMat , 'nscan')
            nscan = spmMat.nscan;
        else
            error('Field nscan doesn''t exist in structure SPM .');
        end

        % number of reference functions
        if isfield(spmMat , 'xX')
            numRefFunc = length(spmMat.xX.name);
        else
            error('Field xX doesn''t exist in structure SPM .');
        end

        errorString = ['SPM design matrix doesn''t match with the data'];

        % Modifying nscan
        % check the nscan with the time points
        % takes care of different or same session timings
        % automatically converts to vector form.
        if length(nscan) == 1
            if nscan == sum(countTimePoints)
                nscan = countTimePoints; % update nscan with time points
                SessInfo = spmMat.Sess;   % session information
                flag_refFunc = 'session_specific'; % pass flag
            else
                if nscan == countTimePoints(data_sessionNumber)
                    data_sessionNumber = 1; % update data_sessionNumber
                else
                    error(errorString);
                end
            end
        else
            lengthScan = length(nscan);
            lengthTimePoints = length(countTimePoints);
            if sum(nscan) == sum(countTimePoints) & lengthScan == lengthTimePoints
                if nscan == countTimePoints
                    flag_refFunc = 'session_specific';
                else
                    error(errorString);
                end
            else
                matchedIndex = find(nscan == countTimePoints(data_sessionNumber)); % find the matching entries
                if isempty(matchedIndex)
                    error(errorString);
                elseif length(matchedIndex) == 1
                    data_sessionNumber = matchedIndex; % update the data_sessionNumber
                    flag_refFunc = 'session_specific'; % pass flag
                end
            end
        end

        if exist('SessInfo', 'var')
            % get correct onsets from the vectorized session information
            if length(nscan) > length(SessInfo) &  length(SessInfo) == 1
                temp = repmat(SessInfo, 1, length(nscan)); % replicate with length of the nscan
                % loop over number of onsets
                for jj = 1:length(SessInfo.U)
                    % get the concatenated onset
                    currentOnset = SessInfo.U(jj).ons;
                    startTp = 1; numberToSubtract = 0;
                    % loop over number of nscan
                    for ii = 1:length(nscan)
                        endTp = sum(nscan(1:ii));
                        subtractedOnset = currentOnset - numberToSubtract; % get the corrected onset for that session
                        temp(ii).U(jj).ons = subtractedOnset(find(subtractedOnset < nscan(ii) & subtractedOnset >= 0)); % get the onsets
                        startTp = endTp + 1;
                        numberToSubtract = endTp;
                    end
                    % end loop for number of nscan
                end
                % end loop for number of onsets
                spmMat.Sess = temp;
                clear temp;
            end
        end

        % clear variables in SPM structure to save memory
        xBF.UNITS = spmMat.xBF.UNITS;
        repetition_time = spmMat.xY.RT;
        spmMat = rmfield(spmMat, 'xBF');
        spmMat.xBF = xBF;
        clear xBF;
        spmMat = rmfield(spmMat, 'xY');
        spmMat.xY.RT = repetition_time;

        % remove all the other fields
        fieldsReq1 = {'nscan', 'xX', 'xBF', 'xY', 'Sess'};
        % include only the required fields
        spmMat = icatb_includeFields_struct(spmMat, fieldsReq1);
        clear fieldsReq1;
        Sess_new = spmMat.Sess;
        spmMat = rmfield(spmMat, 'Sess');

        if isfield(Sess_new, 'U')

            if ~isempty(Sess_new(1).U)
                Sess_new = icatb_includeFields_struct(Sess_new, {'U'});
                fieldsReq1 = {'name', 'ons'};
                % include only the required fields
                for ii = 1:length(Sess_new)
                    Sess_new(ii).U =  icatb_includeFields_struct(Sess_new(ii).U, fieldsReq1);
                end
                clear fieldsReq1;
            end

        end

        spmMat.Sess = Sess_new;
        clear Sess_new;

        % remove unnecessary fields
        XField = spmMat.xX.X;
        nameField = spmMat.xX.name;
        spmMat = rmfield(spmMat, 'xX');
        spmMat.xX.X = XField;
        spmMat.xX.name = nameField;
        clear XField; clear nameField;


        spmMat.nscan = nscan;
        % store the data in spmData
        spmData.SPM = spmMat;
        spmData.nscan = nscan;
        spmData.flag_refFunc = flag_refFunc;
        spmData.data_sessionNumber = data_sessionNumber;
        clear spmMat;
    end
    % end for checking SPM design matrix

    % get session specific regressors if possible
    if strcmp(lower(check_regressors), 'yes')
        if ~exist('spmData', 'var')
            error('spmData structure should be passed');
        end
        spmMat = spmData.SPM; % number of scans
        nscan = spmData.nscan;
        flag_refFunc = spmData.flag_refFunc;
        data_sessionNumber = spmData.data_sessionNumber;
        %%%%%%%%%%%%%%%%%%%
        % selection process

        allNames = spmMat .xX.name; % all reference function names
        % Initialise selected reference function names and indices
        sel_refnames = {};
        sel_indices = [];
        if strcmp(lower(flag_refFunc), 'session_specific') & ~isempty(data_sessionNumber)
            count = 0;
            % pull the session specific or all regressors names
            for ii = 1:numRefFunc
                firstBracket = find(allNames{ii} == '(');
                endBracket = find(allNames{ii} == ')');
                currentNumber = str2num(allNames{ii}(firstBracket(1) + 1:endBracket(1) - 1)); % get the session number
                if currentNumber ==  data_sessionNumber
                    count = count + 1;
                    sel_refnames{count} = allNames{ii};
                    sel_indices(count) = ii;
                end
            end
        end
        % check if the selected reference function names is empty
        if isempty(sel_refnames)
            sel_refnames = allNames;
            sel_indices = (1:length(sel_refnames));
        end

        for ii = 1:length(sel_refnames)
            options(ii).string = sel_refnames{ii};
        end
        sel_refnames = str2mat(options.string); % convert cell to matrix
        clear options;
        % store the data
        spmData.sel_refnames = sel_refnames; % selected reference function names
        spmData.sel_indices = sel_indices; % selected session specific indices
    end
    % end for checking regressors

    % get the regressors of interest
    if strcmp(lower(flag_selecting_regressors), 'yes')
        if ~exist('spmData', 'var')
            error('spmData structure should be passed');
        end
        sel_refnames = spmData.sel_refnames; % selected reference function names
        sel_indices = spmData.sel_indices; % selected session specific indices
        spmMat = spmData.SPM;
        nscan = spmData.nscan;

        % get the indices of the selected regressor names
        if strcmp(lower(analType), 'batch')
            if isfield(design_Info, 'refFunNames')
                selectedRefFunNames = design_Info.refFunNames;
            else
                error('Specify the reference function names for the batch analysis');
            end
            % check the index
            for ii = 1:length(selectedRefFunNames)
                tempIndex = strmatch(lower(selectedRefFunNames{ii}), lower(sel_refnames), 'exact');
                if isempty(tempIndex)
                    error(['Specified reference function with name ', selectedRefFunNames{ii}, ' doesn''t exist']);
                end
                % store the indices
                getIndex(ii) = tempIndex;
            end
        else
            if strcmp(typeStr, 'single')
                [getIndex, name_button] = icatb_listdlg('PromptString', prompt, 'SelectionMode','single',  'ListString', sel_refnames, ...
                    'movegui', 'center', 'windowStyle', 'modal', 'title_fig', title_fig);
            else
                if isempty(maxSelection)
                    [getIndex, name_button] = icatb_listdlg('PromptString', prompt, 'SelectionMode','multiple',...
                        'ListString', sel_refnames, 'movegui', 'center', 'windowStyle', 'modal', ...
                        'title_fig', title_fig);
                else
                    [getIndex, name_button] = icatb_listdlg('PromptString', prompt, 'SelectionMode','multiple',...
                        'ListString', sel_refnames, 'movegui', 'center', 'windowStyle', 'modal', ...
                        'maxSelection', maxSelection, 'title_fig', title_fig);
                end
            end
            if isempty(getIndex)
                error('Reference function name should be selected in order to select the model time course');
            end
        end
        clear options;
        %get information from SPM design matrix
        spmData.selectedRegressors = sel_refnames(getIndex, :); % store the selected regressors
        spmData.getIndex = getIndex;
        spmData.nsess = length(nscan); % get the number of sessions in design matrix

    end
    % end for getting the regressors of interest

    % get the time courses of interest
    if strcmp(lower(get_timecourses), 'yes')
        if ~exist('spmData', 'var')
            error('spmData structure should be passed');
        end
        spmMat = spmData.SPM; % get the SPM structure
        numRefFunc = length(spmMat.xX.name); % number of reference functions
        nscan = spmMat.nscan; % number of scans
        data_sessionNumber = spmData.data_sessionNumber; % session number
        flag_refFunc = spmData.flag_refFunc;
        sel_refnames = spmData.sel_refnames; % selected reference function names
        sel_indices = spmData.sel_indices; % selected session specific indices
        getIndex = spmData.getIndex;
        selectedRegressors = spmData.selectedRegressors;

        % check the size of the time course
        if isfield(spmMat.xX, 'X')
            % pull the time course
            checkX = spmMat.xX.X;
            % check the size of the design matrix
            sizeX = size(checkX);

            if sizeX(1) == 1 & sizeX(2) == numRefFunc*sum(nscan)
                checkX = checkX';
                % check the size of the design matrix
                sizeX = size(checkX);
            end

            if sizeX(1) == numRefFunc & sizeX(2) == sum(nscan)
                checkX = checkX';
                checkX = reshape(checkX, prod(size(checkX)), 1);
                % check the size of the design matrix
                sizeX = size(checkX);
            end

            % modify the model time course as nrows == sum(nscan), ncols = num
            % reference functions
            %disp('Checking SPM.Mat file ...');
            if sizeX(2) == numRefFunc
                %disp('SPM.mat is compatible with the GIFT.');
            elseif sizeX(2) == 1 & sizeX(1) == sum(nscan)*numRefFunc
                numberScans = sum(nscan);
                temp = zeros(numberScans, numRefFunc);
                for ii = 1:numRefFunc
                    temp(:, ii) = checkX(numberScans*(ii - 1) + 1 : numberScans*ii, 1);
                end
                clear checkX;
                %SPMFile = fullfile(inputDir, 'newSPM.mat');
                spmMat.xX.X = temp;
                clear temp;
                disp('Modified the model time course to make compatible with the GIFT');
                % save the modified file
                %save(SPMFile, 'SPM');
                %disp(['Modified the model time course. Please use the new SPM design matrix ', SPMFile, ' to do temporal sorting.']);
            else
                disp('The specified design matrix is not compatible with the GIFT software.');
            end

        else
            error('Field X doesn''t exist in structure SPM .xX');
        end

        % converted SPM time course to (sum(nscans) by numRefFunc) format
        % check if the reference functions are session_specific or all
        if strcmp(lower(flag_refFunc), 'session_specific') & ~isempty(data_sessionNumber)
            selTimecourses = spmMat.xX.X(:, sel_indices(getIndex));
            startTp = 1;
            for ii = 1:data_sessionNumber
                endTp = sum(nscan(1:ii));
                temp = selTimecourses(startTp:endTp, :);
                startTp = endTp + 1;
            end
            timecourse = temp; % modified time course
            spmData.nsess = 1; % number of sessions
            clear temp;
        else
            for ii = 1:length(getIndex)
                % get the data session number
                [sesNum] = icatb_get_onsets(deblank(selectedRegressors(ii, :)), spmMat);
                if sesNum == 1
                    startTp = 1;
                    endTp = nscan(1);
                else
                    startTp = sum(nscan(1 : sesNum - 1)) + 1; % sum until length of nscans - 1
                    endTp = sum(nscan(1:sesNum)); % sum until length of nscans
                end
                timecourse(:, ii) = spmMat.xX.X(startTp:endTp, sel_indices(getIndex(ii)));
            end
        end
        spmData.timecourse = timecourse;
        clear timecourse;
        % Fields required:
        fieldsRequired = {'xBF', 'xY', 'Sess', 'nscan'};
        spmTemp = struct;
        for ii = 1:length(fieldsRequired)
            spmTemp = setfield(spmTemp, fieldsRequired{ii}, getfield(spmMat, fieldsRequired{ii}));
        end
        clear spmMat;
        %     % clear variables in SPM structure to save memory
        spmData.SPM = spmTemp;

    end
    % end for getting time courses of interest

catch

    icatb_displayErrorMsg;

end