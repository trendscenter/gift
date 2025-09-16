function [contrastResults, con_table] = icatb_anova_contrasts(anovaResults, terms, all_contrasts, contrast_sources, tail)
%% Function to do anova contrasts
%
% Inputs:
% 1. anovaResults - Anova Results data structure
% 2. terms - Factors to include
% 3. all_contrasts - Contrast matrix of dimensions unique levels by
% the no. of contrasts
% 4. contrast_sources - Contrast names in a cell array or character array.
% 5. tail - Options are 0 (two_tailed), 1 (right tailed) or -1 (left
% tailed)
%
% Outputs:
% 1. contrastResults - Data structure containing the results of contrast tests
% 2. Contrast table - Contrast table.

%% Get variables from anovaResults data structure
data = anovaResults.data;
anovaResults = rmfield(anovaResults, 'data');
num_main_factors = anovaResults.num_main_factors;

if ~exist('tail', 'var')
    tail = 0;
end


if (max(terms) > num_main_factors)
    error('Error:MaxTerms', 'Max value (%d) in terms variable is greater than the no. of factors (%d)\n', max(terms), num_main_factors);
end

terms = sort(terms(:)');
numSamples = [anovaResults.factor(terms).numSamples];
numSamples = numSamples(:);
sum_sq_err = anovaResults.sum_sq_err;
mean_sq_error = anovaResults.mean_sq_err;
erdf = anovaResults.erdf;

% Total no. of levels
totalLevels = sum([anovaResults.factor(terms).cols_contrast_matrix]);

if (size(all_contrasts, 1) ~= totalLevels)
    error('Error:TotalLevels', 'No. of rows (%d) of contrast Matrix must equal the no. of total levels (%d)\n', size(all_contrasts, 1), ...
        totalLevels);
end

if ~exist('contrast_sources', 'var')
    contrast_sources = gen_con_names(anovaResults, terms, all_contrasts);
else
    % Convert to cell array
    contrast_sources = cellstr(contrast_sources);
    if (length(contrast_sources) ~= size(all_contrasts, 2))
        error('Error:ContrastSources', 'Number of contrast source names (%d) must match the number of contrasts (%d)\n', length(contrast_sources), ...
            size(all_contrasts, 2));
    end
end

%% Beta weights
beta_weights = [anovaResults.factor(terms).beta_weights];
beta_weights = beta_weights(:);

%% Means for each contrast
contrast_means = all_contrasts'*beta_weights;
contrast_means = contrast_means(:)';

%% Standard error for each contrast
std_error = zeros(1, size(all_contrasts, 2));
for nCon = 1:size(all_contrasts, 2)
    currentContrast = all_contrasts(:, nCon);
    std_error(nCon) = sqrt(mean_sq_error*sum((currentContrast.^2)./numSamples));
end

%% T value
tval = contrast_means./std_error;

%% P value
pval = icatb_get_pvalue(tval, erdf, tail);

%% Contrast results
contrastResults.contrast_sources = contrast_sources;
contrastResults.contrast_means = contrast_means;
contrastResults.std_error = std_error;
contrastResults.tval = tval;
contrastResults.pval = pval;

%% Compute F-statistic for multi contrasts if only one term is present
if (length(terms) == 1)
    
    degrees_freedom = rank(all_contrasts);
    
    % Initialise contrast matrix
    contrastMatrix = zeros(length(data), size(all_contrasts, 2));
    
    currentfactor = [anovaResults.factor(terms).currentfactor];
    [unique_factor, unique_inds] = unique(currentfactor);
    unique_factor = currentfactor(sort(unique_inds));
    
    % Loop over contrasts
    for nCon = 1:size(all_contrasts, 2)
        % Loop over levels
        for nL = 1:length(unique_factor)
            if isnumeric(unique_factor)
                temp_ind = (currentfactor == unique_factor(nL));
            else
                temp_ind = strcmpi(currentfactor, unique_factor{nL});
            end
            contrastMatrix(temp_ind, nCon) = all_contrasts(nL, nCon);
        end
        % End loop over levels
    end
    % End loop over contrasts
    
    
    %% Do Regression
    b = pinv(contrastMatrix)*data;
    
    %% line fit
    line_fit = contrastMatrix*b;
    
    %% Sum of squares, degrees of freedom, F stat and p-value
    ss = sum(line_fit.^2);
    ms = ss/degrees_freedom;
    residual = data - line_fit;
    Fstat = ms/mean_sq_error;
    pval = 1 - icatb_spm_Fcdf(Fstat, degrees_freedom, erdf);
    
    %% Multicontrast results
    contrastResults.multi.beta_weights = b(:)'; % betas
    contrastResults.multi.line_fit = line_fit; % Line fit
    contrastResults.multi.sum_sq = ss; % Sum of squares
    contrastResults.multi.df = degrees_freedom; % Degrees of freedom
    contrastResults.multi.mean_sq = ms; % Mean squares
    contrastResults.multi.residual = residual; % residual
    contrastResults.multi.erdf = erdf; % Error degrees of freedom
    contrastResults.multi.sum_sq_err = sum_sq_err; % Sum of squared error
    contrastResults.multi.mean_sq_err = mean_sq_error; % Mean square of error
    contrastResults.multi.Fstat = Fstat; % F statistic
    contrastResults.multi.pval = pval; % P value
    contrastResults.multi.name = 'Contrast';
    
end

%% Table for contrasts
tbl = cell(size(all_contrasts, 2) + 1, 5);

tbl(1, :) = {'Contrast Source', 'Contrast Estimate', 'Standard Error', 'p-value', 'T-value'};
tbl(2:end, 1) = contrast_sources;
tbl(2:end, 2) = cellstr(num2str(contrastResults.contrast_means(:), '%0.3f'));
tbl(2:end, 3) = cellstr(num2str(contrastResults.std_error(:), '%0.3f'));
tbl(2:end, 4) = cellstr(num2str(contrastResults.pval(:), '%0.6f'));
tbl(2:end, end) = cellstr(num2str(contrastResults.tval(:), '%0.3f'));

con_table = {tbl};

clear tbl;

%% Store multi-contrasts if it exists
if isfield(contrastResults, 'multi')
    tbl = repmat({''}, 3, 6);
    tbl(:, 1) = {'Source'; 'Multi Contrast'; 'Error'};
    tbl(:, 2) = {'Sum of Squares'; num2str(contrastResults.multi.sum_sq, '%0.3f'); num2str(contrastResults.multi.sum_sq_err, '%0.3f')};
    tbl(:, 3) = {'df'; num2str(contrastResults.multi.df); num2str(contrastResults.multi.erdf)};
    tbl(:, 4) = {'Mean Square'; num2str(contrastResults.multi.mean_sq, '%0.3f'); num2str(contrastResults.multi.mean_sq_err, '%0.3f')};
    tbl(1:2, 5) = {'F'; num2str(contrastResults.multi.Fstat, '%0.3f')};
    tbl(1:2, end) = {'p-value'; num2str(contrastResults.multi.pval, '%0.6f')};
    con_table{end + 1} = tbl;
end


function contrast_sources = gen_con_names(anovaResults, terms, all_contrasts)
%% Generate contrast source names
%
% Inputs:
% 1. anovaResults - anovaResults data structure
% 2. all_contrasts - All contrasts
%
% Outputs:
% contrast_sources - Contrast source names
%

%% Initialise contrast names
contrast_sources = repmat({''}, size(all_contrasts, 2), 1);

%% Loop over contrasts
for nCon = 1:size(all_contrasts, 2)
    % Initialise some variables
    startTp = 1;
    endTp = 0;
    countT = 0;
    % Loop over terms
    for nTerm = 1:length(terms)
        endTp = endTp + length(anovaResults.factor(terms(nTerm)).numSamples);
        indices = (startTp:endTp);
        conInd = find(all_contrasts(indices, nCon) ~= 0);
        % Loop over current contrast indices
        for nInd = 1:length(conInd)
            countT = countT + 1;
            currentConValue = all_contrasts(indices(conInd(nInd)), nCon);
            
            % Operator
            op = '+';
            if (sign(currentConValue) == -1)
                op = '-';
            end
            
            temp_name = ['(', anovaResults.factor(terms(nTerm)).name, '(', num2str(conInd(nInd)), ') = ', ...
                    num2str(abs(currentConValue)), ')'];
            
            % Add negative sign for the first character
            if (countT == 1)
                if (sign(currentConValue) == -1)
                    temp_name = [op, temp_name];
                end
                temp_con_name = temp_name;
            else
                temp_con_name = [temp_con_name, ' ', op, ' ', temp_name];
            end
            
        end
        % End loop over current contrast indices
        startTp = endTp + 1;
    end
    % End loop over terms
    contrast_sources{nCon} = temp_con_name;
end
%% End loop over contrasts
