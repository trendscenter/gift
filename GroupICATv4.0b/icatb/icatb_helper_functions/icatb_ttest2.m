function [p_value, t_value, df] = icatb_ttest2(x, y, tail)
%% Compute two sample t-test between samples x and y
%
% Inputs:
% 1. x - Sample A
% 2. y - Sample B
% 3. tail - Options are 0 (two_tailed), 1 (right tailed) or -1 (left
% tailed)
%
% Outputs:
% 1. p_value - P value
% 2. t_value - T Value
% 3. df - Degrees of freedom

%% Convert both samples to column vectors
x = x(:);
y = y(:);

if ~exist('tail', 'var')
    tail = 0;
end

%% Remove NaN's in x
x(isnan(x)) = [];

if isempty(x)
    error('Data (x) is empty or it has all NaN''s');
end

%% Remove NaN's in y
y(isnan(y)) = [];

if isempty(y)
    error('Data (y) is empty or it has all NaN''s');
end

% Sample size of x and y
n = length(x);
m = length(y);

a = (n - 1)*(std(x)^2);
b = (m - 1)*(std(y)^ 2);

% Numerator
num = (sqrt(n*m))*(mean(x) - mean(y));

% Denominator
denom = sqrt((a + b)*(n + m) / (n + m - 2));

%% T-value
t_value = num/denom;

% Degrees of freedom
df = (n + m - 2);

%% P value
p_value = icatb_get_pvalue(t_value, df, tail);