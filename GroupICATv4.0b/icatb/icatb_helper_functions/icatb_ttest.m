function [p_value, t_value, df] = icatb_ttest(x, tail)
%% Compute one sample t-test of sample x
%
% Inputs:
% 1. x - Sample
% 2. tail - Options are 0 (two_tailed), 1 (right tailed) or -1 (left
% tailed)
%
% Outputs:
% 1. p_value - P value
% 2. t_value - T value
% 3. df - Degrees of freedom

% Convert x to column vector
if (numel(x) == length(x))
    x = x(:);
end

%% Remove NaN's
%chkNaN = prod(x, 2);
%x(isfinite(chkNaN)==0, :)=[];
%x(isnan(x)) = [];

if isempty(x)
    error('Data is empty or it has all NaN''s');
end

if ~exist('tail', 'var')
    tail = 0;
end

% Degrees of freedom
%df = size(x, 1) - 1;

N = sum(isfinite(x));
df = N - 1;

%% T-value
%t_value = mean(x) ./ (std(x) ./ sqrt(size(x, 1)));
t_value = icatb_nanmean(x) ./ (icatb_nanstd(x) ./ sqrt(N));

%% P value
p_value = icatb_get_pvalue(t_value, df, tail);

