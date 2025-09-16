function d = icatb_pdist(X, flagP, dist_metric)
%% Pairwise distance
%
% Inputs:
% 1. X - Matrix of dimensions m x n.
% 2. flagP - Options are 1 and 2.
%   1 - All combinations of rows
%   2 - Successive rows
% 3. dist_metric - Options are 1 and 2.
%   1 - L1 norm
%   2 - L2 norm


%% Succesive rows or all combinations
if (~exist('flagP', 'var') || isempty(flagP))
    flagP = 1;
end

if (ischar(flagP))
    if (strcmpi(flagP, 'successive'))
        flagP = 2;
    else
        % all
        flagP = 1;
    end
end

%% Check Euclidean (L2) or Manhattan (L1) distance
if (~exist('dist_metric', 'var') || isempty(dist_metric))
    dist_metric = 1;
end


if (ischar(dist_metric))
    if (strcmpi(dist_metric, 'manhattan'))
        dist_metric = 1;
    else
        % Euclidean distance
        dist_metric = 2;
    end
end

if (flagP == 1)
    pairs = nchoosek(1:size(X, 1), 2);
else
    pairs = [(1:size(X, 1)-1)', (2:size(X, 1))'];
end

d = zeros(1, size(pairs, 1));
for np = 1:length(d)
    d(np) = norm(X(pairs(np, 1), :) - X(pairs(np, 2), :), dist_metric);
end

