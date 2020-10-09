function anova_tbl = icatb_anova_table(anovaResults)
%% Function to generate table from anova results data structure
%
% Inputs:
% 1. anovaResults - Anova Results data structure
%
% Outputs:
% anova_tbl - Anova table
%

%% No. of rows contain Title + (Main Effects + Interaction) + Error Term + Total
num_rows = length(anovaResults.factor) + 3;

%% Column headings
col_headings = {'Source', 'Sum Sq. (Type III)', 'df', 'Mean Sq.', 'F', 'p-value'};

%% No. of cols
num_cols = length(col_headings);

%% Initialise cell array
anova_tbl = repmat({''}, num_rows, num_cols);

%% Headings of table
anova_tbl(1, :) = col_headings;

%% Total sum of squares and degrees of freedom
anova_tbl(end, 1:3) = {'Total', num2str(anovaResults.sum_sq_corrected_total, '%0.3f'), num2str(anovaResults.df)};

%% Error term, degrees of freedom, Mean Sq.
anova_tbl(end - 1, 1:4) = {'Error', num2str(anovaResults.sum_sq_err, '%0.3f'), num2str(anovaResults.erdf), num2str(anovaResults.mean_sq_err, '%0.3f')};

%% Loop over terms from 2 to end - 2
for nF = 1:length(anovaResults.factor)
    anova_tbl(nF + 1, :) = {anovaResults.factor(nF).name, num2str(anovaResults.factor(nF).sum_sq, '%0.3f'), num2str(anovaResults.factor(nF).df), ...
            num2str(anovaResults.factor(nF).mean_sq, '%0.3f'), num2str(anovaResults.factor(nF).Fstat, '%0.3f'), ...
            num2str(anovaResults.factor(nF).pval, '%0.6f')};
end
%% End loop over terms