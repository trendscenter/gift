function [Value_scalar] = icatb_read_scalar_file(fid, keywd, returnType)
% reads parameters file and gets the scalar variables

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
                    charPos = length(tline) - equalPos - 1;

                    if strcmp(lower(returnType), 'integer')
                        [currentStr, strVar] = strread(tline, ['%', num2str(equalPos), 'c %s']);

                        %%% Check if the variable is a integer
                        Value_scalar = icatb_check_variable_integer(deblank(strVar{1}));
                        if isempty(Value_scalar)
                            error('The variable entered is not a integer');
                        end

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    elseif strcmp(lower(returnType), 'float')
                        [currentStr, Value_scalar] = strread(tline, ['%', num2str(equalPos), 'c %f']);
                    else
                        [currentStr, Value_scalar] = strread(tline, ['%', num2str(equalPos), 'c %', num2str(charPos), 'c']);
                        Value_scalar = deblank(Value_scalar);
                        if strcmpi(returnType, 'numeric')
                            Value_scalar = str2num(Value_scalar);
                        end
                    end

                end

            end
            % end for checking the commented lines or not

        end
        % end for checking character tlineNew

    end
    % end for checking empty state of tlineNew variable

end

% Place the position of the file at the beginning of the record
frewind(fid);

if ~exist('Value_scalar', 'var')
    %fclose(fid);
    error(['Cannot find the keyword ', keywd, ' in the given file.']);
end