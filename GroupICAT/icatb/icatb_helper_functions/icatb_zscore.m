function tc = icatb_zscore(tc)
%% Convert data to z-scores
%

if isvector (tc)
    tc = detrend(tc, 0) / std(tc);
    return;
end

for n = 1:size(tc, 2)
    tc(:, n) = detrend(tc(:, n), 0) ./ std(tc(:, n));
end