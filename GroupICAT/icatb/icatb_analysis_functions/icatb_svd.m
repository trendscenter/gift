function [V, Lambda] = icatb_svd(data, numpc, varargin)
%% Do singular value decomposition on data using economy size decomposition
%
% Inputs:
% 1. data - Data is 2D double array
% 2. numpc - No. of principal components
%
% Outputs:
% 1. V - Eigen vectors
% 2. Lambda - Eigen values in a diagonal matrix
%

solverType = 'selective';
minDim = min(size(data));

verbose = 1;

%% Loop over arguments
for i = 1:2:length(varargin)
    if (strcmpi(varargin{i}, 'solver'))
        solverType = lower(varargin{i + 1});
    elseif (strcmpi(varargin{i}, 'verbose'))
        verbose = varargin{i + 1};
    end
end

if (~exist('numpc', 'var'))
    numpc = minDim;
end

if (verbose)
    disp('Using Singular Value Decomposition to do PCA');
end

if (numpc > minDim)
    if (verbose)
        fprintf('The selected no. of components (%d) is larger than the minimum no. of data dimensions (%d)\n', numpc, minDim);
    end
    numpc = minDim;
end

if (numpc == minDim)
    solverType = 'all';
end

eigRetainedDisp = 0;
if (strcmpi(solverType, 'selective') && isa(data, 'double'))
    opts = struct('disp', 0, 'maxit', 1000);
    [U, Lambda, V] = svds(data, numpc, 'L', opts);
else
    eigRetainedDisp = 1;
    if (icatb_get_matlab_version > 13)
        [U, Lambda, V] = svd(data, 'econ');
    else
        [U, Lambda, V] = svd(data, 0);
    end
end

%% Sort eigen vectors in Ascending order
Lambda = diag(Lambda);
Lambda = Lambda ./ sqrt(size(data, 1) - 1);
[Lambda, inds] = sort(Lambda.^2);
V = V(:, inds);

sumAll = sum(Lambda);

%% Return only the extracted components
V = V(:, end-numpc+1:end);
Lambda = Lambda(end-numpc+1:end);
sumUsed = sum(Lambda);
Lambda = diag(Lambda);

retained = (sumUsed/sumAll)*100;

if (eigRetainedDisp && verbose)
    fprintf('%g%% of (non-zero) eigenvalues retained.\n', retained);
end

if (verbose)
    fprintf('Done\n\n');
end