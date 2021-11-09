function [t, p, stats] = icatb_mT_reduce(stats_all, terms, term, options)
%% T-test using collection of stats files run using mancova
%
%
% Inputs:
% 1. stats_all: Cell array of stats structures generated from mT
% 2. terms: Terms in the model
% 3. term: Term to test
%
% Outputs:
% 1. t - T-stat
% 2. p - P-values
% 3. stats - stats structure
%

if ~exist('options', 'var')
    options = cell(0);
end


if ~isempty(strmatch('verbose', options, 'exact'))
    fprintf('\n')
    fprintf('Computing t-statistics ...\n')
    fprintf('\n')
    tMT = tic;
end

%% Aggregate cov(X, X) and cov(X, Y)
n = 0;
for nS = 1:length(stats_all)
    stats_u = stats_all{nS};
    n = n + size(stats_u.X, 1);
    if (nS == 1)
        covXX = zeros(size(stats_u.X, 2), size(stats_u.X, 2));
        covXY = zeros(size(stats_u.X, 2), size(stats_u.B, 2));
    end
    tmpXX = stats_u.X'*stats_u.X;
    covXX = covXX + tmpXX;
    covXY = covXY + (tmpXX*stats_u.B);
end


%% Estimate beta weights and SSE
invXTX = pinv(covXX);
B = covXX \ covXY;
SSE = 0;
X = cell(length(stats_all));
for nS = 1:length(stats_all)
    stats_u = stats_all{nS};
    res = stats_u.X*B - stats_u.X*stats_u.B;
    SSE = SSE + stats_u.SSE + sum(res.^2);
    X{nS} = stats_u.X;
end

X = cat(1, X{:});


r = min([n, size(stats_u.X, 2)]);
I = mFindTerms(term, terms);

if isempty(I)
    error('The specified term was not found in the model.')
end

% Create a matrix of indices for contrasts.

contrasts = [ I(:) , zeros(length(I), 1) ];

% Create a matrix of indices for pairs of contrasts.

combinations = mNC2(length(I));

% Concatenate all contrasts into a single matrix of indices.

contrasts = cat(1, contrasts, I(combinations));

DFE = n - r;
MSE = SSE / DFE;

t = NaN(size(contrasts,1), size(B, 2));
p = NaN(size(contrasts,1), size(B, 2));

for j = 1 : size(B, 2)
    
    s = MSE(j) * invXTX;
    
    for i = 1 : size(contrasts, 1)
        
        if contrasts(i, 2) == 0
            
            % Test if a regression coefficient is different from zero.
            
            if ~isempty(strmatch('verbose', options, 'exact')) && j == 1
                fprintf('\tT-Test %d for Y(:, %d): B(%d, %d) == 0 ... ', ...
                    i, j, contrasts(i, 1), j)
            end
            
            t(i, j) = B(contrasts(i, 1), j) / sqrt(s(contrasts(i, 1), contrasts(i,1)));
            p(i, j) = 2 * (1 - mTCDF(abs(t(i, j)), n - r));
            
            if ~isempty(strmatch('verbose', options, 'exact')) && j == 1
                fprintf('p = %g\n', p(i, j));
            end
            
        else
            
            if ~isempty(strmatch('verbose', options, 'exact')) && j == 1
                fprintf('\tT-Test %d for Y(:, %d): B(%d, %d) == B(%d, %d) ... ', ...
                    i, j, contrasts(i, 1), j, contrasts(i, 2), j)
            end
            
            % Test if two regression coefficinets are different from each other.
            
            t(i, j) = (B(contrasts(i, 1), j) - B(contrasts(i, 2), j)) / sqrt(s(contrasts(i, 1), contrasts(i, 1)) + ...
                s(contrasts(i, 2), contrasts(i, 2)) - 2 * s(contrasts(i, 1), contrasts(i, 2)));
            
            p(i, j) = 2 * (1 - mTCDF(abs(t(i, j)), n - r));
            
            if ~isempty(strmatch('verbose', options, 'exact')) && j == 1
                fprintf('p = %g\n', p(i, j));
            end
            
        end
    end
    
    if ~isempty(strmatch('verbose', options, 'exact')) && j == 1
        fprintf('\n')
    end
    
end

levels = contrasts;
levels(levels > 0) = levels(levels > 0) - levels(1,1) + 1;

stats.Terms  = terms;
stats.X      = X;
stats.B      = B;
stats.SSE    = SSE;
stats.DFE    = DFE;
stats.MSE    = MSE;
stats.Term   = term;
stats.Levels = levels;

if ~isempty(strmatch('verbose', options, 'exact'))
    fprintf('... finished in %f seconds.\n', toc(tMT))
    fprintf('\n')
end