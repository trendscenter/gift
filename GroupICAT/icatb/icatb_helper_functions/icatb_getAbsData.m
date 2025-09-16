function [data] = icatb_getAbsData(data)
% return the absolute value of the data if the data is complex

% check the data
if isstruct(data)
    if isfield(data, 'mag')
    data = data.mag;
    else
        error('Magnitude field is missing in data');
    end
else
    if ~isreal(data)
        data = abs(data);
    end
end
% end for checking data type