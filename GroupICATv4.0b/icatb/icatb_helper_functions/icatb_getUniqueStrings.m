function [uniqueStr, countVec] = icatb_getUniqueStrings(inputStr)
% returns unique strings from a set of repeated strings

if ~isempty(inputStr)

    if ~iscell(inputStr)
        % make a cell array
        if ischar(inputStr)
            inputStr = cellstr(inputStr);
        else
            error('Not a string or cell array');
        end
    end
    
    % sort by rows
    [sortedRows, index] = sortrows(str2mat(inputStr));
    kk = 1; differentIndex = 0; 
    matchedIndices = {};
    jj = 1;
    while jj <= length(inputStr)
        % return matching indices
        matchedIndices{1} = strmatch(sortedRows(jj, :), sortedRows, 'exact');
        % Length of the matching indices
        differentIndex = length(matchedIndices{1});
        indx = matchedIndices{1}(1);        
        % store the first of the matched indices
        countVec(kk) = indx;        
        kk = kk + 1;
        % Update jj
        jj = jj + differentIndex;
        clear matchedIndices; clear differentIndex;
    end    
    % vector and unique strings
    countVec = sort(index(countVec));    
    uniqueStr = cell(1, length(countVec));
    % get the unique strings
    for ii = 1:length(countVec)
		uniqueStr{ii} = inputStr{countVec(ii)};
    end
    % end for getting the unique strings    
else
    uniqueStr = {};
    countVec = 0;
end