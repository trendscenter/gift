function [inputText] = icatb_checkInputFile(inputFile, inputText)
% Read from the input file the variables that match the tag in the input
% text and set the answer string for edit or value for the popup controls.
%
% Input:
% 1. inputFile - M file containing the necessary variables to run the batch
% script
% 2. inputText - inputText is the structure containing the information
% regarding user interface controls and the values
%
% Output:
% inputText - Updated inputText in the field answerString for edit or value
% for popup

[pathstr, inputF, extn] = fileparts(inputFile);

oldDir = pwd;

controlsTobeRead = find([inputText.read_from_file] == 1);

% Loop over the required controls
for ii = 1:length(controlsTobeRead)
    % current tag
    currentTag = inputText(controlsTobeRead(ii)).tag;
    % data type
    dataType  = inputText(controlsTobeRead(ii)).dataType;
    % UI type
    uiType = inputText(controlsTobeRead(ii)).uiType;
    % answer string
    answerString = str2mat(inputText(controlsTobeRead(ii)).answerString);
    % Current flag (Vector or scalar)
    currentFlag = inputText(controlsTobeRead(ii)).flag;
    
    if strcmp(currentTag, 'prefix')
       dataType = 'output_prefix';
    end
    
    % Read input data
    [inputData, status, message] = icatb_read_variables(inputFile, currentTag, currentFlag, dataType, 'no_display_error');
    currentData = getfield(inputData, currentTag);
    
    if strcmp(currentTag, 'prefix')
        if status == 0
            error(message);
        end
    end
    
    if isempty(currentData) & ~strcmp(currentTag, 'maskFile')
        disp(message);
        error(['Error:', currentTag], ['Check input file: %s \nto make sure that variable %s exists or is of the required data type'], ...
            inputFile, currentTag);
    end
    
    % For edit controls
    if strcmpi(uiType, 'edit')
        % convert to string
        if isnumeric(currentData)
            currentData = num2str(currentData);
        end
        inputText(controlsTobeRead(ii)).answerString = currentData;
        
    elseif strcmpi(uiType, 'popup') | strcmpi(uiType, 'listbox')
        
        % Set the value for popup control
        if strcmp(currentTag, 'maskFile')
            
            if isnumeric(currentData) & ~isempty(currentData)
                error(['Error:', currentTag], ['Check input file: %s \nto make sure that variable %s is either empty or a valid file'], ...
                    inputFile, currentTag);
            end
            
            if ~isempty(currentData)
                if ~strcmp(currentData, '[]')
                    inputText(controlsTobeRead(ii)).value = 2;
                else
                    inputText(controlsTobeRead(ii)).value = 1;
                end
            else
                inputText(controlsTobeRead(ii)).value = 1;
            end
            
        elseif strcmp(currentTag, 'numReductionSteps')
            % Special case for data reduction steps                
            matchedIndex = strmatch(num2str(currentData), answerString, 'exact');
            
            inputText(controlsTobeRead(ii)).value = matchedIndex;
            
        elseif strcmp(currentTag, 'scaleType')
            % Special case for scale type                
            
            if isnumeric(currentData)
                
                inputText(controlsTobeRead(ii)).value = currentData + 1;
                
                if size(answerString, 1) < currentData
                    disp(message);
                   error(['Error:', currentTag], ['Value for the variable %s \nin the input file: %s exceeds the maximum'], ...
                        currentTag, inputFile);
                end
                
            else
                disp(message);
                error(['Error:', currentTag], ['Value for the variable %s \nin the input file: %s must be an integer'], ...
                        currentTag, inputFile);
                
            end
            
        else
            
            % Check the data type for the pop up control
            if isnumeric(currentData)
                
                if size(answerString, 1) < currentData
                    disp(message);
                    error(['Error:', currentTag], ['Value for the variable %s \nin the input file: %s exceeds the maximum'], ...
                        currentTag, inputFile);
                elseif currentData == 0
                    inputText(controlsTobeRead(ii)).value = 1;
                else
                    inputText(controlsTobeRead(ii)).value = round(currentData);
                end
            else
                % for popup set the value
                matchIndex = strmatch(lower(currentData), lower(answerString), 'exact');
                if ~isempty(matchIndex)
                    inputText(controlsTobeRead(ii)).value = matchIndex;
                end
            end
            % end for checking the data type for the pop up control
        end
        % end for setting the value for popup control
    end
    % end for checking controls
    
end
% End loop over the required controls