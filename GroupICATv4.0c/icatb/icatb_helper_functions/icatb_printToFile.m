function icatb_printToFile(fileName, varStruct, titleStr, printType, optional, additionalInfo)
%
% fields for varStruct are tag and value

try
    if ~exist('printType', 'var')
        printType = 'column_wise';
    end

    if ~exist('optional', 'var')
        optional = 'overwrite';
    end

    if ~exist('varStruct', 'var')
        error('varStruct variable is not passed');
    end

    if ~isstruct(varStruct)
        error('varStruct shoule be of structure data type');
    end
    
    if ~exist('additionalInfo', 'var')
        additionalInfo = [];
    end

    if strcmpi(optional, 'overwrite')
        fid = fopen(fileName, 'w');        
        fprintf(fid, '%s\n\n', titleStr);                
    else
        % append
        fid = fopen(fileName, 'a');
        fprintf(fid, '\n');
        fprintf(fid, '%s\n\n', titleStr);
    end
    
    % additional info
    if ~isempty(additionalInfo)
        for ii = 1:size(additionalInfo, 1)
            fprintf(fid, '%s\n', deblank(additionalInfo(ii, :)));        
        end
    end

    % get the tags
    for ii = 1:length(varStruct)
        tags{ii} = varStruct(ii).tag;
    end

    % if print type is column wise
    if strcmpi(printType, 'column_wise')
        for ii = 1:length(varStruct)
            tempVar = varStruct(ii).value;
            if ~ischar(tempVar)            
                tempVar = num2str(tempVar, 8);
            end
            stringToPrint(ii).string =  str2mat(tags{ii}, tempVar);
            clear tempVar;
        end
    else
        firstStr = str2mat(tags);
        stringToPrint(1).string = firstStr;
        firstVal = varStruct(1).value;

        % collect entries from each variable and form a string
        for ii = 1:size(firstVal, 1)
            for jj = 1:length(tags)
                % temporary variable
                tempVar = varStruct(jj).value(ii, :);
                if ~ischar(tempVar)                         
                    tempVar = num2str(tempVar, 8);
                end
                concatVar(jj).name = tempVar;
                clear tempVar;
            end
            stringToPrint(ii + 1).string = str2mat(concatVar.name);
            clear concatVar;
        end
    end
    % end for checking the printType

    %%%%%%%%%%%%% Print strings row by row %%%%%%%%%%%%%%%%
    % Loop over number of rows involved
    % print the number of columns with tab after each string is printed
    for ii = 1:size(stringToPrint(1).string, 1)
        % loop over the columns
        for jj = 1:length(stringToPrint)
            % get the first row of each string
            fprintf(fid, '%s\t', stringToPrint(jj).string(ii, :));
        end
        if strcmpi(printType, 'column_wise') & ii == 1
            fprintf(fid, '\n');
        end
        fprintf(fid, '\n');
    end
    fclose(fid);
catch
    fclose(fid);
end