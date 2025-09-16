function [dataStr, uniqueInd] = icatb_uniqueStr(dataStr, type)
% Returns unique strings

if ~isempty(dataStr)
    useLowerCase = 0;
    if exist('type', 'var')
        if strcmpi(type, 'dir')
            if ispc
                useLowerCase = 1;
            end
        end
    end

    if ischar(dataStr)
        dataStr = cellstr(str2mat(dataStr));
    end

    for nD = 1:length(dataStr)
        if ~useLowerCase
            ind = strmatch(dataStr{nD}, dataStr, 'exact');
        else
            ind = strmatch(lower(dataStr{nD}), lower(dataStr), 'exact');
        end
        if length(ind) > 1
            if exist('type', 'var') && strcmpi(type, 'dir')
                dataStr(ind(1:end-1)) = {''};
            else
                dataStr(ind(2:end)) = {''};
            end
        end
    end

    indices = icatb_good_cells(dataStr);
    uniqueInd = find(indices);
    dataStr = dataStr(uniqueInd);
end