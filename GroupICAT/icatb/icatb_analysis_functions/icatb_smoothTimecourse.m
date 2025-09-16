function icaTimecourse = icatb_smoothTimecourse(icaTimecourse, diffTimePoints, num_DataSets, numComp, ...
    smoothPara, smoothValue)
% apply smoothing defaults to the timecourse

% observations must be a column vector
if size(icaTimecourse, 1) == 1
    icaTimecourse = icaTimecourse';
end

nPoints = size(icaTimecourse, 1); % total number of points (including the concatenated set along rows)

if exist('diffTimePoints', 'var')
    if isempty(diffTimePoints)
        clear diffTimePoints;
    end
end

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

% smooth the time courses
if strcmp(lower(smoothPara), 'yes')
    tempVar = zeros(size(icaTimecourse)); % initialise vars
    for nComp = 1:numComp
        startTp = 1;
        for tp = 1:length(numTimePoints)
            endTp = sum(numTimePoints(1:tp));
            % apply the gaussian smoothing to the time courses
            tempVar(startTp:endTp, nComp) = icatb_gauss_smooth1D(icaTimecourse(startTp:endTp, nComp), smoothValue);
            % update the starting time point
            startTp = endTp + 1;
        end
    end
    icaTimecourse = tempVar;
    clear tempVar;
end