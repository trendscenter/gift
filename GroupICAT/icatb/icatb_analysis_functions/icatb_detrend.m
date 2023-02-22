function [y] = icatb_detrend(y, num_DataSets, diffTimePoints, detrendNumber)
% performs detrending based on the detrend number selected
% in icatb_defaults file selected
%
% Input: observation vector or matrix
%
% Output: detrended observations based on the detrend number specified in
% the defaults

if ~exist('detrendNumber', 'var')
    % use the defaults here
    icatb_defaults;
    global DETRENDNUMBER;
    detrendNumber = DETRENDNUMBER;
end

if exist('diffTimePoints', 'var')
    if isempty(diffTimePoints)
        clear diffTimePoints;
    end
end

% observations must be a column vector
if size(y, 1) == 1
    y = y';
end

nPoints = size(y, 1); % total number of points (including the concatenated set along rows)

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

% different time points must be a vector
if prod(size(numTimePoints)) ~= length(numTimePoints)
    error('different time points must be a vector');
end

if ~exist('detrendNumber', 'var')
    detrendNumber = 0;
    disp('-- DETRENDNUMBER doesn''t exist. By default selecting DETRENDNUMBER = 0.');
end

if detrendNumber > 3 | round(detrendNumber) - detrendNumber ~=0
    detrendNumber = 0;
    disp('-- DETRENDNUMBER selected is not in options or is not a valid integer. By default selecting DETRENDNUMBER = 0.');
end

% Initialise remove Trend
removeTrend = zeros(size(y));

% Start grouping the nuisance parameters as
% models
switch detrendNumber
    % Detrend with 0 or remove mean
    case 0
        nTerms = 1;
        rowStart = 1;
        % Initialise beta coefficients
        b = zeros(nTerms*num_DataSets, size(y, 2));
        % Loop over number of data sets
        for nDataSets = 1:num_DataSets
            rowEnd = sum(numTimePoints(1:nDataSets));
            X = [ones(numTimePoints(nDataSets), 1)];
            % loop over dimension 2
            for nColumns = 1:size(y, 2)
                % solve beta coeff
                b(nTerms*(nDataSets - 1) + (1:nTerms), nColumns) = icatb_regress(y(rowStart:rowEnd, nColumns), X);
                % Calculating how much is to be subtracted from time course
                removeTrend(rowStart:rowEnd, nColumns) = X*b(nTerms*(nDataSets - 1) + (1:nTerms), nColumns);
            end
            %end for loop over columns
            % update the row starting
            rowStart = 1 + rowEnd;
            clear X;
        end
        % Removes the mean and linear trend
    case 1
        nTerms = 2;
        rowStart = 1;
        % Initialise beta coefficients
        b = zeros(nTerms*num_DataSets, size(y, 2));
        % Loop over number of data sets
        for nDataSets = 1:num_DataSets
            rowEnd = sum(numTimePoints(1:nDataSets));
            rampFun = icatb_unitRamp((1:numTimePoints(nDataSets))');
            rampFun = detrend(rampFun, 0); % remove the mean before using ramp Function
            X = [rampFun, ones(numTimePoints(nDataSets), 1)];
            % loop over dimension 2
            for nColumns = 1:size(y, 2)
                % solve beta coeff
                b(nTerms*(nDataSets - 1) + (1:nTerms), nColumns) = icatb_regress(y(rowStart:rowEnd, nColumns), X);
                % Calculating how much is to be subtracted from time course
                removeTrend(rowStart:rowEnd, nColumns) = X*b(nTerms*(nDataSets - 1) + (1:nTerms), nColumns);
            end
            % update the row starting
            rowStart = 1 + rowEnd;
            clear X;
            clear rampFun;
        end
        % Uses sine and cosine one cycle, removes mean and
        % linear trend
    case 2
        nTerms = 4;
        rowStart = 1;
        % Initialise beta coefficients
        b = zeros(nTerms*num_DataSets, size(y, 2));
        % Loop over number of data sets
        for nDataSets = 1:num_DataSets
            rowEnd = sum(numTimePoints(1:nDataSets));
            % steps in sin and cosine
            unitPoint_sin2phi = (2*pi/(numTimePoints(nDataSets) - 1));
            vec2phi = (0:unitPoint_sin2phi:2*pi)';
            rampFun = icatb_unitRamp((1:numTimePoints(nDataSets))');
            rampFun = detrend(rampFun, 0); % remove the mean before using ramp Function
            X = [detrend(sin(vec2phi), 0), detrend(cos(vec2phi), 0), rampFun, ...
                ones(numTimePoints(nDataSets), 1)];
            % loop over dimension 2
            for nColumns = 1:size(y, 2)
                % solve beta coeff
                b(nTerms*(nDataSets - 1) + (1:nTerms), nColumns) = icatb_regress(y(rowStart:rowEnd, nColumns), X);
                % Calculating how much is to be subtracted from time course
                removeTrend(rowStart:rowEnd, nColumns) = X*b(nTerms*(nDataSets - 1) + (1:nTerms), nColumns);
            end
            % update the row starting
            rowStart = 1 + rowEnd;
            clear X;
            clear rampFun;
        end
        % Uses sine and cosine two cycles plus sine and cosine one cycle, removes mean and
        % linear trend
    case 3
        nTerms = 6;
        rowStart = 1;
        % Initialise beta coefficients
        b = zeros(nTerms*num_DataSets, size(y, 2));
        % Loop over number of data sets
        for nDataSets = 1:num_DataSets
            rowEnd = sum(numTimePoints(1:nDataSets));
            % steps in sin and cosine
            unitPoint_sin2phi = (2*pi/(numTimePoints(nDataSets) - 1));
            unitPoint_sin4phi = (4*pi/(numTimePoints(nDataSets) - 1));
            vec2phi = (0:unitPoint_sin2phi:2*pi)';
            vec4phi = (0:unitPoint_sin4phi:4*pi)';
            rampFun = icatb_unitRamp((1:numTimePoints(nDataSets))');
            rampFun = detrend(rampFun, 0); % remove the mean before using ramp Function
            X = [detrend(sin(vec4phi), 0), detrend(cos(vec4phi), 0), detrend(sin(vec2phi), 0), ...
                detrend(cos(vec2phi), 0), rampFun, ones(numTimePoints(nDataSets), 1)];
            % loop over dimension 2
            for nColumns = 1:size(y, 2)
                % solve beta coeff
                b(nTerms*(nDataSets - 1) + (1:nTerms), nColumns) = icatb_regress(y(rowStart:rowEnd, nColumns), X);
                % Calculating how much is to be subtracted from time course
                removeTrend(rowStart:rowEnd, nColumns) = X*b(nTerms*(nDataSets - 1) + (1:nTerms), nColumns);
            end
            % end for loop
            % update the row starting
            rowStart = 1 + rowEnd;
            clear X;
            clear rampFun;
        end

end
% end switch for detrendNumber

% Trend is removed in timecourse
y = y - removeTrend;