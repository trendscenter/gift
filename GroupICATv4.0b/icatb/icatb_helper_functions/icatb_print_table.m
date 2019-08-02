function icatb_print_table(tbl, file_name, file_mode, closeFile)
%% Print table to file
%
% Inputs:
% 1. tbl - Table of type cell array
% 2. file_name - File name
% 3. file_mode - File mode
% 4. closeFile - Options are 0 or 1

%% Initial check
if ~iscell(tbl)
    error('Table must be a cell array');
end

if ~exist('file_mode', 'var')
    file_mode = 'w+';
end

if ~exist('closeFile', 'var')
    closeFile = 1;
end

%% Open file
if (ischar(file_name))
    fid = fopen(file_name, file_mode);
    if (fid == -1)
        error('Error:FileOpen', 'File %s cannot be opened in mode %s\n', file_name, file_mode);
    end
else
    fid = file_name;
end

%% Add one empty line for append mode
if (strcmpi(file_mode, 'a+'))
    fprintf(fid, '\n');
end

% Format for string
formatStr = repmat('%s\t\t', 1, size(tbl, 2));

try

    %% Initialise new table
    newTbl = cell(size(tbl, 1), size(tbl, 2));

    %% Loop over columns
    for nCol = 1:size(tbl, 2)
        currentCol = str2mat(tbl(:, nCol));
        % Loop over rows
        for nRow = 1:size(currentCol, 1)
            newTbl{nRow, nCol} = currentCol(nRow, :);
        end
        % End loop over rows
    end
    %% End loop over columns

    clear tbl;

    %% Loop over rows
    for nRow = 1:size(newTbl, 1)
        currentRow = newTbl(nRow, :);
        fprintf(fid, formatStr, currentRow{:});
        fprintf(fid, '\n');
    end
    %% End loop over rows

    if (closeFile)
        fclose(fid);
    end

catch

    %% Rethrow the error
    msg = lasterror;

    if (closeFile)
        try
            fclose(fid);
        catch
        end
    end

    rethrow(msg);

end
