function [rSquare_stat, b, ModelIndices, otherIndices, linearRegress, removeTrend, ica, subject_partial_corr, ...
    partialCorrSlopes] = ...
    icatb_multipleRegression(model, ica, num_Regress, num_DataSets, diffTimePoints, detrendNumber)
% Performs multiple regression
% Outputs: R-Square Statistic, Coefficients, Model Indices, Other than the
% model indices, line fit to model observed data

% inputs are modelTimecourse and icaTimecourse, number of regressors,
% number of data sets

% By default DETRENDNUMBER is 0

% DETRENDNUMBER - case 0 - Removes the mean
% DETRENDNUMBER - case 1 - Removes the mean and linear trend
% DETRENDNUMBER - case 2 - Uses sine and cosine one cycle, removes mean and
% linear trend
% DETRENDNUMBER - case 3 - Uses sine and cosine two cycles plus sine and cosine one cycle, removes mean and
% linear trend


if ~exist('detrendNumber', 'var')
    % use the defaults here
    icatb_defaults;
    global DETRENDNUMBER;
    detrendNumber = DETRENDNUMBER;
end

if nargin < 2
    error('Need atleast two arguments for the multiple linear regression');
end


% Convert to column vector
if size(ica, 1) == 1
    ica = ica';
end

% number of points taken into consideration
nPoints = size(ica, 1);

% check the dimensions
if size(model, 1) ~= nPoints & size(model, 2) ~= nPoints
    error('Model dimensions doesn''t match with that of observed data');
elseif size(model, 1) ~= nPoints
    model = model';
end


% If the number of regressors doesn't exist
if ~exist('num_Regress', 'var')
    num_Regress = size(model, 2);
    %     disp('-- Number of Regressors is not selected. By default checking the size of the design matrix');
end

% By default number of data sets is 1
if ~exist('num_DataSets', 'var')
    num_DataSets = 1;
    %     disp('-- Number of data sets is not selected. By default selecting data sets equal to 1.');
end

if ~exist('detrendNumber', 'var')
    detrendNumber = 0;
    disp('-- detrendNumber doesn''t exist. By default selecting detrendNumber = 0.');
end

if detrendNumber > 3 | round(detrendNumber) - detrendNumber ~=0
    detrendNumber = 0;
    disp('-- detrendNumber selected is not in options or is not a valid integer. By default selecting detrendNumber = 0.');
end

if ~exist('diffTimePoints', 'var')
    % number of Time Points
    numTimePoints = nPoints / num_DataSets;
    numTimePoints = repmat(numTimePoints, 1, num_DataSets);
else
    numTimePoints = diffTimePoints;
end

%% Concatenate the design matrix using the number of
%% regressors and the data sets and depending on the detrend number

removeTrend = zeros(size(ica));

linearRegress = zeros(size(ica));

switch detrendNumber
    case 0
        nTerms = 1;
    case 1
        nTerms = 2;
    case 2
        nTerms = 4;
    otherwise
        nTerms = 6;
end

totalTerms = nTerms + num_Regress;
residual = zeros(size(ica));
% Initialise beta coefficients
b = zeros(totalTerms*num_DataSets, 1);
% collect the indices other than the model indices
otherIndices = zeros(1, nTerms*num_DataSets);
ModelIndices = zeros(1, num_Regress*num_DataSets);
rowStart = 1;
% Initialise individual subject partial correlations
subject_partial_corr = zeros(num_DataSets*num_Regress, 1);
partialCorrSlopes = subject_partial_corr;

countDataSet = 0;
% Loop over number of data sets
for nDataSets = 1:num_DataSets

    rowEnd = sum(numTimePoints(1:nDataSets));

    % Regressor matrix
    [X] = icatb_modelX(model(rowStart:rowEnd, :), numTimePoints(nDataSets), detrendNumber);

    % Calculate regression
    [ica, b, removeTrend, linearRegress, rowStart, rowEnd, ModelIndices, otherIndices] = ...
        returnParameters(nDataSets, X, num_Regress, totalTerms, ica, removeTrend, b, linearRegress, ...
        rowStart, rowEnd, nTerms, otherIndices, ModelIndices);

    % get the corresponding model indices
    model_ind = ModelIndices((nDataSets - 1)*num_Regress + 1 : num_Regress*nDataSets);

    residual(rowStart:rowEnd) = ica(rowStart:rowEnd) - X(:, 1:num_Regress)*b(model_ind, 1);

    %%%%%%%%%%%%% Individual Partial Correlations %%%%%%%%%%%%

    % Loop over regressors
    for nRegress = 1:num_Regress
        countDataSet = countDataSet + 1;
        % get the other indices
        other_indices = find(model_ind(nRegress) ~=  model_ind);

        if ~isempty(other_indices)
            % Remove other regressor contribution
            currentTc = ica(rowStart:rowEnd) - X(:, other_indices)*b(model_ind(other_indices), 1);
        else
            currentTc = ica(rowStart:rowEnd);
        end
        % Partial correlation slopes and correlation
        [partialCorrSlopes(countDataSet), partial_rsquare] = icatb_regress(currentTc, X(:, nRegress));
        subject_partial_corr(countDataSet) = sqrt(abs(partial_rsquare));
        clear currentTc;
        clear other_indices;
    end
    % end loop over regressors

    clear model_ind;
    %%%%%%%%%%%%% End for Individual Partial Correlations %%%%%%%%%%%%

    % update the row starting
    rowStart = 1 + rowEnd;


    clear X;
    clear rampFun;

end
% End loop over data sets


% R-square statistic
rSquare_stat = 1 - (norm(residual, 2) / norm(ica - mean(ica), 2))^2;


function rSquare_stat = rSquareStat(x, model, b)
% Rsquare stat


residual = x - model*b;

rSquare_stat = 1 - (norm(residual, 2) / norm(x - mean(x), 2))^2;


function [ica, b, removeTrend, linearRegress, rowStart, rowEnd, ModelIndices, otherIndices] = ...
    returnParameters(nDataSets, X, num_Regress, totalTerms, ica, removeTrend, b, linearRegress, ...
    rowStart, rowEnd, nTerms, otherIndices, ModelIndices)

% modelX which contains the concatenated model
%modelX(rowStart:rowEnd, (nDataSets - 1)*size(X, 2) + 1 :  nDataSets*size(X, 2)) = X;
% calculate other indices
otherIndex = totalTerms*(nDataSets - 1) + (num_Regress + 1 : nTerms + num_Regress);
% calculate model indices
modelIndex = totalTerms*(nDataSets - 1) + (1:num_Regress);
% solve beta coeff
b(totalTerms*(nDataSets - 1) + (1:totalTerms), 1) = icatb_regress(ica(rowStart:rowEnd, 1), X);

% Calculating how much is to be subtracted from time course
removeTrend(rowStart:rowEnd, 1) = X(:, num_Regress + 1:totalTerms)*b(otherIndex);

% Trend is removed in timecourse
ica(rowStart:rowEnd, 1) = ica(rowStart:rowEnd, 1) - removeTrend(rowStart:rowEnd, 1);

% use the detrended ica time course and solve the coefficients
% for model
b(modelIndex, 1) = icatb_regress(ica(rowStart:rowEnd, 1), X(:, 1:num_Regress));

% calculate linefit
linearRegress(rowStart:rowEnd, 1) = X(:, 1:num_Regress)*b(modelIndex, 1);

otherIndices(1, (nDataSets-1)*nTerms + 1 : nTerms*nDataSets) = totalTerms*(nDataSets - 1) + ...
    (num_Regress + 1 : nTerms + num_Regress);
ModelIndices(1, (nDataSets-1)*num_Regress + 1 : num_Regress*nDataSets) = totalTerms*(nDataSets - 1) + (1:num_Regress);