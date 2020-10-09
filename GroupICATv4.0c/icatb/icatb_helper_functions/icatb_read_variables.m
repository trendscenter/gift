function [inputData, status, message] = icatb_read_variables(file_name, keywd, varFlag, varCriteria, dispError)
% reads variable or variables from file with the specified keywd

% get the file name and extensions
[pathstr, fName, extns] = fileparts(file_name);

% initialise input data structure
inputData = struct;

% old directory
oldDir = pwd;

currentVar = {};

if ~exist('dispError', 'var')
    dispError = 'display_error';
end

% check the input file extension
if ~strcmpi(extns, '.m')

    fid = fopen(file_name, 'r');

    if fid == -1
        error(['File ', file_name, ' doesn''t exist']);
    end

    try
        % read current variable
        if strcmpi(varFlag, 'scalar')
            currentVar = icatb_read_scalar_file(fid, keywd, 'string');
            % get the status and message
            [status, message] = icatb_errorCheck(currentVar, varCriteria, keywd);
            if status == 0
                error(message);
            end
        else
            currentVar = icatb_read_vector_file(fid, keywd, 'string');
        end

        if strcmpi(varCriteria, 'numeric') | strcmpi(varCriteria, 'integer') | strcmpi(varCriteria, 'float')
            currentVar = str2num(currentVar);
        end

        % set the field to the input data
        inputData = setfield(inputData, keywd, currentVar);
        clear currentVar;

    catch
        % close the file
        fclose(fid);
        icatb_displayErrorMsg;
    end

else

    if ~isempty(pathstr)
        % change to new directory
        cd(pathstr);
    end

    % evaluate the input text
    eval(fName);
    % change to old directory
    cd(oldDir);

    try
        currentVar = eval(keywd);
    catch
        error(['Error:', keywd], ...
            'Variable %s doesn''t exist or \nis not of the required data type \nin file: %s', ...
            keywd, file_name);
    end

    if isnumeric(currentVar)
        if strcmpi(varFlag, 'scalar')
            if ~isempty(currentVar)
                currentVar = currentVar(1);
            end
        end
    end

    checkVar = currentVar;

    % check if is a valid integer or not
    if strcmpi(varCriteria, 'integer') | strcmpi(varCriteria, 'float') | strcmpi(varCriteria, 'numeric')
        checkVar = num2str(checkVar);
    end

    % get the status and message
    [status, message] = icatb_errorCheck(checkVar, varCriteria, keywd);
    if (status == 0) & strcmpi(dispError, 'display_error')
        error(message);
    end
    % end for getting the correct data type
    inputData = setfield(inputData, keywd, currentVar);
    clear currentVar; clear checkVar;

end