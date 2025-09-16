function tc_normalize = icatb_normalizeTimecourse(tc_normalize, ref_vector, diffTimePoints, num_DataSets)
% Normalize time course

% observations must be a column vector
if size(tc_normalize, 1) == 1
    tc_normalize = tc_normalize';
end

nPoints = size(tc_normalize, 1); % total number of points (including the concatenated set along rows)

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

% check the dimensions
if size(ref_vector, 1) == 1
    ref_vector = ref_vector';
end

tempVar = zeros(size(tc_normalize)); % initialise vars
startTp = 1;
for tp = 1:length(numTimePoints)
    endTp = sum(numTimePoints(1:tp));
    % normalize the time courses
    I = ref_vector(startTp:endTp, :); % reference vector
    meanI = mean(I); % mean of reference vector
    normI = norm(I - meanI); % norm of the detrended ref. vector
    % loop over number of regressors
    for ii = 1:size(tc_normalize, 2)
        if all(tc_normalize(startTp:endTp, ii) == 1)
            %zeroMM = modelTimecourse(:, ii);
            zeroMM = tc_normalize(startTp:endTp, ii) - mean(tc_normalize(startTp:endTp, ii));
            tempVar(startTp:endTp, ii) = zeroMM;
        else
            zeroMM = tc_normalize(startTp:endTp, ii) - mean(tc_normalize(startTp:endTp, ii));
            normM = norm(zeroMM);
            if normM == 0
                disp(['norm of the column ',  num2str(ii), ' is zero. Not calculating normalization for this column.']);
                return;
            else
                tempVar(startTp:endTp, ii) = meanI + (zeroMM/normM * normI);
            end
            % end for checking
        end
    end
    % update the starting time point
    startTp = endTp + 1;
    clear I;
    clear meanI; 
    clear zeroMM; 
    clear normM; clear normI;
end
tc_normalize = tempVar;
clear tempVar;