function [X, orig_mean] = icatb_remove_mean(X, verbose)
%% Remove mean of each column and return original mean values in a vector.
%
% Inputs:
% 1. X - X is 2D double array.
% 2. verbose - Options are 0 and 1.
%
% Outputs:
% 1. X - Data with zero mean
% 2. orig_mean - Original mean values
%

if (~exist('verbose', 'var'))
    verbose = 1;
end

orig_mean = mean(X);

try
    if (~isempty(which('bsxfun.m')))
        %% Use bsxfun if possible
        X = bsxfun(@minus, X, orig_mean);
    else
        tempVar = repmat(orig_mean, size(X, 1), 1);
        X = X - tempVar;
        clear tempVar;
    end
catch
    clear tempVar;
    if (verbose)
        disp('Using slow but less memory for removing mean in the data');
    end
    %% Remove mean of sample one at a time
    for nM = 1:size(X, 2)
        X(:, nM) = X(:, nM) - orig_mean(nM);
    end
end

