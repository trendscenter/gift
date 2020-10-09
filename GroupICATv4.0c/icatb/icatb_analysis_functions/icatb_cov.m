function A = icatb_cov(X, removeMean)
%% Compute covariance matrix. If out of memory error occurs, covariance matrix is
% computed with least memory requirements.
%
% Inputs:
% 1. X - Data of size m by n
% 2. removeMean - Remove mean. Options are 0 and 1.
%
% Outputs:
%
% 1. A - Covariance matrix
%

if (~exist('X', 'var'))
    error('Data variable is not passed');
end

%% Use defaults
if (~exist('removeMean', 'var'))
    removeMean = 1;
end

[m, n] = size(X);

%% Remove mean
if (removeMean)
    X = icatb_remove_mean(X, 1);
end

try
    A = X'*X;
catch
    disp('Using slower way to compute covariance matrix ...');
    blocks = 1000;
    loops = ceil(m/blocks);
    if (isa(X, 'single'))
        A = zeros(n, n, 'single');
    else
        A = zeros(n, n);
    end
    endBlock = 0;
    for nblock = 1:loops
        startBlock = endBlock + 1;
        endBlock = min([endBlock + blocks, m]);
        A = A + (X(startBlock:endBlock, :)'*X(startBlock:endBlock, :));
    end
end

%% Return the final result
if (m > 1)
    A = A ./ (m - 1);
else
    A = A ./ m;
end
