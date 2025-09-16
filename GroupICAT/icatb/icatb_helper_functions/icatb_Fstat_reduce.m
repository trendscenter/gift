function stats = icatb_Fstat_reduce(stats_all)
%% Compute Mancova F-statistic
%
% Inputs:
% 1. stats_all: Cell array of stats structures generated from mT
% 2. terms: Terms in the model
% 3. names - Regressor names
%
% Outputs:
% stats - Statistics Table
%


uniqueTerms = stats_all{1}.terms;
names = stats_all{1}.regressors;

terms = [{0}, uniqueTerms];

X = cell(length(stats_all));
for nS = 1:length(stats_all)
    X{nS} = stats_all{nS}.X;
end

X = cat(1, X{:});

numVoxels = size(stats_all{1}.YTY, 2);

FANCOVAN = repmat(NaN, length(uniqueTerms), numVoxels);
pANCOVAN = repmat(NaN, length(uniqueTerms), numVoxels);
SSER = pANCOVAN;
dFR = zeros(1, length(uniqueTerms));

dFF = size(X, 1) - rank(X);

[YTMY, YTY, X] = getSSE(stats_all);

SSEF = YTMY;

for i = 1 : length(uniqueTerms)
    I0 = mTerms(uniqueTerms{i}, terms);
    X0 = X(:, I0);
    dFR(i) = size(X0, 1) - rank(X0);
    tmp = getSSE(stats_all, I0, i);
    SSER(i, :) = tmp;
    dFR(i) = size(X0, 1) - rank(X0);
    FANCOVAN(i, :) = ((SSER(i, :) - SSEF) ./ (dFR(i) - dFF)) ./ (SSEF / dFF);
    try
        pANCOVAN(i, :) = 1 - mFCDF(FANCOVAN(i, :), dFR(i) - dFF, dFF);
    catch
    end
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
%stats.B = B;
%stats.B0 = B0;
stats.X = X;


function varargout = getSSE(stats_all, I0, termNo)
%% Get sum of squared residual
%

betaField = 'B';
SSEField = 'SSEF';
if (exist('I0', 'var'))
    betaField = 'B0';
    SSEField = 'SSER';
end

%% Aggregate cov(X, X) and cov(X, Y)
n = 0;
for nS = 1:length(stats_all)
    stats_u = stats_all{nS};
    XM = stats_u.X;
    if (exist('I0', 'var'))
        XM = XM(:, I0);
    end
    Bhat = stats_u.(betaField);
    if (exist('termNo', 'var'))
        Bhat = Bhat{termNo};
    end
    n = n + size(XM, 1);
    
    if (nS == 1)
        covXX = zeros(size(XM, 2), size(XM, 2));
        covXY = zeros(size(XM, 2), size(Bhat, 2));
        if (nargout > 1)
            YTY = zeros(1, size(Bhat, 2));
        end
    end
    tmpXX = XM'*XM;
    covXX = covXX + tmpXX;
    covXY = covXY + (tmpXX*Bhat);
    if (nargout > 1)
        YTY = YTY + stats_u.YTY;
    end
end


%% Estimate beta weights and SSE
B = covXX \ covXY;
SSE = 0;
if (nargout > 1)
    X = cell(length(stats_all));
end
for nS = 1:length(stats_all)
    stats_u = stats_all{nS};
    XM = stats_u.X;
    if (exist('I0', 'var'))
        XM = XM(:, I0);
    end
    Bhat = stats_u.(betaField);
    sse = stats_u.(SSEField);
    if (exist('termNo', 'var'))
        Bhat = Bhat{termNo};
        sse = sse(termNo, :);
    end
    res = XM*B - XM*Bhat;
    SSE = SSE + sse + sum(res.^2);
    if (nargout > 1)
        X{nS} = XM;
    end
end

if (nargout > 1)
    X = cat(1, X{:});
end

varargout{1} = SSE;
if (nargout > 1)
    varargout{2} = YTY;
end

if (nargout > 2)
    varargout{3} = X;
end
