function pval = icatb_get_pvalue(tval, df, tail)
%% Compute p-value significance
%
% Inputs:
% 1. tval - Tvalue
% 2. df - Degrees of freedom
% 3. tail -  Options are 0 (two_tailed), 1 (right tailed) or -1 (left
% tailed)
%
% pval - P value significance
%

%% P value
if (tail == 0)
    % Two tailed t-test
    pval = 2*icatb_spm_Tcdf(-abs(tval), df);
elseif (tail == 1)
    % Right tailed ttest
    pval = icatb_spm_Tcdf(-tval, df);
elseif (tail == -1)
    % Left tailed ttest
    pval = icatb_spm_Tcdf(tval, df);
else
    error('Tail must be either 0, 1 or -1');
end