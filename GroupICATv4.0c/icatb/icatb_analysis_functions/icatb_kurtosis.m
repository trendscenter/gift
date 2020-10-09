function [kurtosisValue, detrended_data] = icatb_kurtosis(data, num_DataSets, diffTimePoints, detrendNumber)
% compute kurtosis or fourth order statistic
% check if the different time points option is passed

[nRows, nCols] = size(data);

% convert to column vector
if nRows == 1
    data = data';
end

[nRows, nCols] = size(data);

nPoints = size(data, 1); % total number of points (including the concatenated set along rows)
 
% give the number of data sets involved in the
% observation vector
if ~exist('num_DataSets', 'var')
    num_DataSets = 1;
end

% get diffTimePoints vector
% otherwise replicate the diffTimepoints with the 
% number of datasets
if ~exist('diffTimePoints', 'var')  
    numTimePoints = ceil(nPoints / num_DataSets);
    numTimePoints = repmat(numTimePoints, 1, num_DataSets);
else
    if length(diffTimePoints) ~= num_DataSets
        error('length of different time points must be the same as the number of data sets');
    end
    numTimePoints = diffTimePoints; % contains the length of the time points of each vector along rows
    clear diffTimePoints;
end

if exist('detrendNumber', 'var')
    % detrending data based on the DETRENDNUMBER specified in icatb_defaults.m
    [detrended_data] = icatb_detrend(data, num_DataSets, numTimePoints, detrendNumber);
else
    [detrended_data] = icatb_detrend(data, num_DataSets, numTimePoints);
end

% loop over the size of the data
for ii = 1:size(data, 2)   
    kurtosisValue(ii) =  sum((detrended_data(:, ii)./std(data(:, ii))).^4) / (nRows) - 3;
end 

