function [correlationValue, timecourse, b, ModelIndices] = icatb_correlateFunctions(modelTimecourse, timecourse, num_DataSets, diffTimePoints, ...
    detrendNumber)
% computes the correlation between the model time course and the observed
% time course. Current version uses multiple regression to calculate the correlation
% values. Sign is determined by multiplying model time course with ICA
% 
% Input: model time course, ica time course, number of data sets used to
% concatenate data, different time points
% 
% Output: correlation value 

% convert to column vector if possible
if size(modelTimecourse, 1) == 1
    modelTimecourse = modelTimecourse';
end

% convert to column vector if possible
if size(timecourse, 1) == 1
    timecourse = timecourse';
end

nPoints = size(timecourse, 1);

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

% sum x * sum y
numerator = modelTimecourse'*timecourse; 
% Multiple regression is used to determine the r-square statistic as it has
% more general way of incorporating regressors and nuisance parameters
if exist('detrendNumber', 'var')
    [rSquare_stat, b, ModelIndices, otherIndices, linearRegress, removeTrend, timecourse] = icatb_multipleRegression(modelTimecourse, timecourse, 1, num_DataSets, numTimePoints, detrendNumber);
else
    [rSquare_stat, b, ModelIndices, otherIndices, linearRegress, removeTrend, timecourse] = icatb_multipleRegression(modelTimecourse, timecourse, 1, num_DataSets, numTimePoints);
end
% correlation value is determined by multiplying sign with the square root
% of r-square statistic.
correlationValue = sign(numerator)*sqrt(rSquare_stat);
