function [Value_vector] = icatb_read_vector_file(fid, keywd, returnType)
% reads vectors from a file

% loop until end of file is encountered
while ~feof(fid)

    % Read the current line
    tline = fgets(fid);

    % read the string
    tlineNew = strread(tline, '%s');

    if ~isempty(tlineNew)
        tlineNew = tlineNew{1};

        % check again the character array
        if ~isempty(tlineNew)

            % Check if the lines are commented or not
            if ~strcmp(tlineNew(1), '%')

                % check if the existing keyword exists
                FindChar = findstr(lower(tline), lower(keywd));

                if ~ isempty(FindChar)

                    equalPos = findstr(tline, '=');
                    % read the entire text from '=' to the end
                    % Including white space characters
                    [currentStr] = sscanf(tline(equalPos + 1:length(tline)), '%c');
                    currentStr = deblank(currentStr);

                    % collect the values in a vector depending on the return type
                    if strcmp(lower(returnType), 'integer')
                        [StrVar] = strread(currentStr, '%s', 'delimiter', ',');

                        % Check if the values entered are integers
                        for ii = 1:length(StrVar)
                            currentVar = check_variable_integer(deblank(StrVar{ii}));
                            if isempty(currentVar)
                                error('The variable entered is not a integer');
                            end
                            Value_vector(ii, 1) = currentVar;
                        end

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    elseif strcmp(lower(returnType), 'float')
                        [Value_vector] = strread(currentStr, '%f', 'delimiter', ',');
                    else
                        [Value_vector] = strread(currentStr, '%s', 'delimiter', ',');
                    end
                    % end for checking the class of the variable

                end

            end
            % end for checking the commented lines

        end
        % end for checking character tlineNew

    end
    % end for checking empty state of tlineNew

end

% Place the position of the file at the beginning of the record
frewind(fid);

if ~exist('Value_vector', 'var')
    %fclose(fid);
    error(['Cannot find the keyword ', keywd, ' in the given file.']);
end