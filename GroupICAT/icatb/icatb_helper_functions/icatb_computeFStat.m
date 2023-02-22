function stats = icatb_computeFStat(U, X, terms, names)
%% Compute Mancova F-statistic
%
% Inputs:
% U - data of dimensions subjects x voxels
% X - Design matrix
% terms - Terms
% names - Regressor names
%
% Outputs:
% stats - Statistics Table
%

if (numel(U) == length(U))
    U = U(:);
end

YTY = zeros(1, size(U, 2));
SSEF = zeros(1, size(U, 2));
SSER = zeros(length(names), size(U, 2));
B = zeros(size(X, 2), size(U, 2));
FANCOVAN = repmat(NaN, length(names), size(U, 2));
pANCOVAN = repmat(NaN, length(names), size(U, 2));
B0_a = cell(size(U, 2), length(names));

%% Fit at each voxel
for j = 1:size(U, 2)
    u = U(:, j);
    Xr = X(~isnan(u), :);
    u = u(~isnan(u));
    stats_r = compute_fstat(u, Xr, terms, names);
    FANCOVAN(:, j) = stats_r.FANCOVAN;
    pANCOVAN(:, j) = stats_r.pANCOVAN;
    uniqueTerms = stats_r.terms;
    SSEF(j) = stats_r.SSEF;
    SSER(:, j) = stats_r.SSER;
    YTY(j) = stats_r.YTY;
    dFR = stats_r.dFR;
    dFF = stats_r.dFF;
    B(:, j) = stats_r.B;
    B0_a(j, :) = stats_r.B0;
end


B0 = cell(1, length(names));
for i = 1:length(B0)
    tmp = B0_a(:, i);
    B0{i} = cat(2, tmp{:});
end



%% Stats output
stats.terms = uniqueTerms;
stats.regressors = names;
stats.FANCOVAN = FANCOVAN;
stats.pANCOVAN = pANCOVAN;
stats.SSEF = SSEF;
stats.dFF = dFF;
stats.SSER = SSER;
stats.dFR = dFR;
stats.YTY = YTY;
stats.B = B;
stats.B0 = B0;
stats.X = X;

function stats = compute_fstat(U, X, terms, names)
 
 
[~, inds]= unique(names);
names = names(sort(inds));

uniqueTerms = mUnique(terms(2:end));
M        = X * pinv(X' * X) * X';

FANCOVAN = repmat(NaN, length(uniqueTerms), size(U, 2));
pANCOVAN = repmat(NaN, length(uniqueTerms), size(U, 2));
SSER = pANCOVAN;
dFR = zeros(1, length(uniqueTerms));

B0 = cell(1, length(uniqueTerms));
dFF = size(X, 1) - rank(X);
SSEF = sum(U .* ((eye(size(M)) - M) * U));

B = pinv(X'*X)*X'*U;

for i = 1 : length(uniqueTerms)
    I0 = mTerms(uniqueTerms{i}, terms);
    X0 = X(:, I0);
    M0 = X0 * pinv(X0' * X0) * X0';
    SSER(i, :) = sum(U .* ((eye(size(M0)) - M0) * U));
    dFR(i) = size(X0, 1) - rank(X0);
    FANCOVAN(i, :) = ((SSER(i, :) - SSEF) ./ (dFR(i) - dFF)) ./ (SSEF / dFF);
    try
        pANCOVAN(i, :) = 1 - mFCDF(FANCOVAN(i, :), dFR(i) - dFF, dFF);
    catch
    end
    B0{i} = pinv(X0'*X0)*X0'*U;
    %YTXO{i} = U'*X0;
    %     for j = 1 : size(U, 2)
    %         [ FANCOVAN(i, j), pANCOVAN(i, j), SSEF(i, j), dFF(i, j), SSER(i, j), dFR(i, j)] = mF(U(:, j), X, X0, M, M0);
    %     end
end

%% Stats output
stats.terms = uniqueTerms;
stats.regressors = names;
stats.FANCOVAN = FANCOVAN;
stats.pANCOVAN = pANCOVAN;
stats.SSEF = SSEF;
stats.dFF = dFF;
stats.SSER = SSER;
stats.dFR = dFR;
stats.YTY = sum(U.^2);
stats.B = B;
stats.B0 = B0;
%stats.YTX = YTX;
%stats.YTXO = YTXO;
stats.X = X;