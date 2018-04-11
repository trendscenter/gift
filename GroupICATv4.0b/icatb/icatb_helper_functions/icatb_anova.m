function anovaResults = icatb_anova(data, factors, varargin)
%% Compute N-way anova
%
% Inputs:
% 1. data - Data sample must be a vector. Otherwise data will be treated as
% vector
% 2. factors - factors must be in a cell array.
% 3. varargin - Arguments after factors must come in pairs like
% 'model', 'interaction', 'var_names', {'X1'; 'X2'; 'X3'}
%
%
% Outputs:
% anovaResults - Anova results data structure.
%
% Example:
% y = [52.7 57.5 45.9 44.5 53.0 57.0 45.9 44.0]';
% g1 = [1 2 1 2 1 2 1 2];
% g2 = {'hi';'hi';'lo';'lo';'hi';'hi';'lo';'lo'};
% g3 = {'may'; 'may'; 'may'; 'may'; 'june'; 'june'; 'june'; 'june'};
% anovaResults = icatb_anova(y, {g1, g2, g3}, 'model', 'interaction', 'var_names', {'G1'; 'G2'; 'G3'});
% tbl = icatb_anova_table(anovaResults);
%

if (nargin < 2)
    error('Anova requires atleast two parameters. Arguments after factors must come in pairs.');
end

num_parameters = nargin - 2;

if (mod(num_parameters, 2) ~= 0)
    error('Arguments must come in pairs after factors variable');
end

if ~iscell(factors)
    factors = {factors};
end

% Convert to column vector
data = data(:);

check_len_factors = cellfun('length', factors);

if (any(check_len_factors ~= length(data)))
    error('Please check factors variable as length of vector in each factor must match the length of the data');
end

%% Initialise variables
interaction = 0;
varNames = strcat('X', num2str((1:length(factors))'));
varNames = cellstr(varNames);

%% Loop over parameters
for nVar = 1:2:num_parameters
    if strcmpi(varargin{nVar}, 'model')
        if strcmpi(varargin{nVar + 1}, 'interaction')
            interaction = 1;
        end
    elseif strcmpi(varargin{nVar}, 'var_names')
        varNames = varargin{nVar + 1};
    end
end
%% End loop over parameters

varNames = cellstr(varNames);

if (length(varNames) ~= length(factors))
    error('Error:VarNames', 'Number of factor names (%d) must match the number of factors (%d)\n', length(varNames), length(factors));
end

%% Remove NaN's in the data
checkNaNs = isnan(data);
data(checkNaNs) = [];
if isempty(data)
    error('Data is empty or it has all NaN''s');
end

%% Initialise Anova Results data structure
anovaResults = struct;
anovaResults.data = data;

%% Remove mean and compute sum of squares
sum_sq_total = sum(data.^2);
data = detrend(data, 0);
sum_sq_corrected_total = sum(data.^2);
intercept_term = sum_sq_total - sum_sq_corrected_total;

num_main_factors = length(factors);
anovaResults.num_main_factors = num_main_factors;
anovaResults.sum_sq_total = sum_sq_total;
anovaResults.sum_sq_corrected_total = sum_sq_corrected_total;
anovaResults.df = length(data) - 1;
anovaResults.ss_intercept = intercept_term;
anovaResults.intercept_df = 1;

%% Form design matrix
num_combinations = length(factors);
modelMatrix = eye(num_main_factors);
if interaction
    if (num_main_factors > 1)
        allFactors = (1:num_main_factors);
        % All Possible combinations
        all_comb = nchoosek(allFactors, 2);
        num_combinations = num_combinations + size(all_comb, 1);
        tempM = modelMatrix;
        clear modelMatrix;
        % Initialise Model matrix
        modelMatrix = zeros(num_combinations, num_main_factors);
        % Set the Main effects to identity matrix
        modelMatrix(1:num_main_factors, :) = tempM;
        clear tempM;
        % Loop over interactions
        for nM = 1:size(all_comb, 1)
            % Set the other combinations to 1 as well
            modelMatrix(num_main_factors + nM, all_comb(nM, :)) = 1;
        end
        % End loop over interactions
    end
end

%% Initialise factor data structure
anovaResults.factor = repmat(struct('contrastMatrix', [], 'cols_contrast_matrix', [], 'numSamples', [], 'currentfactor', [], 'beta_weights', ...
    [], 'sum_sq', [], 'df', [], 'mean_sq', [], 'name', [], 'Fstat', [], 'pval', []), ...
    1, num_combinations);

%% Main effects
for nF = 1:length(factors)
    
    currentfactor = factors{nF};
    
    if isnumeric(currentfactor)
        currentfactor = currentfactor(:);
    else
        % Treat them as characters
        if ~iscell(currentfactor)
            currentfactor = cellstr(currentfactor);
        end
    end
    
    % Remove the levels that has NaN's in the data
    currentfactor(checkNaNs) = [];
    
    if (length(currentfactor) ~= length(data))
        error('Error:DataLength', 'Factor %d length (%d) must be the same as data length (%d)\n', nF, length(currentfactor), length(data));
    end
    
    % Get unique factors
    [unique_factor, unique_inds] = unique(currentfactor);
    unique_factor = currentfactor(sort(unique_inds));
    
    if (length(unique_factor) <= 1)
        error('Error:Levels', 'You need to specify atleast two levels for factors %d to work with anova', nF);
    end
    
    % Loop over unique levels
    contrastMatrix = logical(zeros(length(data), length(unique_factor)));
    numSamples = zeros(1, length(unique_factor));
    cell_means = zeros(1, length(unique_factor));
    for nU = 1:length(unique_factor)
        if isnumeric(unique_factor)
            temp_ind = (currentfactor == unique_factor(nU));
        else
            temp_ind = strcmpi(currentfactor, unique_factor{nU});
        end
        numSamples(nU) = length(find(temp_ind == 1));
        temp_ind = temp_ind(:);
        cell_means(nU) = mean(anovaResults.data(temp_ind));
        contrastMatrix(:, nU) = temp_ind;
    end
    % End loop over unique levels
    
    % Store some important fields
    anovaResults.factor(nF).currentfactor = currentfactor;
    anovaResults.factor(nF).cell_means = cell_means;
    anovaResults.factor(nF).contrastMatrix = contrastMatrix;
    anovaResults.factor(nF).numSamples = numSamples;
    anovaResults.factor(nF).cols_contrast_matrix = size(contrastMatrix, 2);
    
    anovaResults.factor(nF).df = length(unique_factor) - 1;
    anovaResults.factor(nF).terms = repmat(nF, 1, length(numSamples));
    anovaResults.factor(nF).name = varNames{nF};
    
end

%% Model interactions
if interaction
    if (length(factors) > 1)
        
        allFactors = (1:length(factors));
        
        nFactor = 2;
        %while (nFactor <= length(factors))
        % All Possible combinations
        %all_comb = nchoosek(allFactors, nFactor);
        
        countFactor = num_main_factors;
        %% Loop over combinations
        for nComb = 1:size(all_comb, 1)
            countFactor = countFactor + 1;
            clear contrastMatrix;
            % Current combination
            currentComb = all_comb(nComb, :);
            temp_varname = [];
            for nM = 1:length(currentComb) - 1
                if (nM == 1)
                    contrastMatrix = icatb_elem_by_elem_mult(anovaResults.factor(currentComb(nM)).contrastMatrix, ...
                        anovaResults.factor(currentComb(nM+1)).contrastMatrix);
                    temp_varname = [varNames{currentComb(nM)}, ' * ', varNames{currentComb(nM + 1)}];
                else
                    contrastMatrix = icatb_elem_by_elem_mult(contrastMatrix, anovaResults.factor(currentComb(nM+1)).contrastMatrix);
                    temp_varname = [temp_varname, ' * ', varNames{currentComb(nM + 1)}];
                end
            end
            
            % Degrees of freedom
            df = prod([anovaResults.factor(currentComb).df]);
            anovaResults.factor(countFactor).cols_contrast_matrix = size(contrastMatrix, 2);
            anovaResults.factor(countFactor).contrastMatrix = contrastMatrix;
            anovaResults.factor(countFactor).df = df;
            anovaResults.factor(countFactor).terms = repmat(countFactor, 1, size(contrastMatrix, 2));
            anovaResults.factor(countFactor).name = temp_varname;
            
        end
        %% End loop over combinations
        
        %   nFactor = nFactor + 1;
        %end
        
    end
    
end

%% Make anova constraint matrix to handle overdetermined system of
% equations
constraint_matrix = icatb_anova_constraint(modelMatrix, [anovaResults.factor(1:num_main_factors).cols_contrast_matrix]);

% Get full design
contrastMatrix = [anovaResults.factor.contrastMatrix];

% Remove contrast matrix
anovaResults.factor = rmfield(anovaResults.factor, 'contrastMatrix');

% Add Intercept term
contrastMatrix = [contrastMatrix, ones(length(data), 1)];

%% Fit the model with the data using QR decomposition
[beta_weights, line_fit] = icatb_constrained_lse(contrastMatrix, data, constraint_matrix);

beta_weights = beta_weights(:)';

% Sum of squares of full model
sum_sq_full_model = sum(line_fit.^2);

anovaResults.beta_weights = beta_weights;
anovaResults.line_fit = line_fit;
anovaResults.sum_sq_full_model = sum_sq_full_model;

% Add intercept term
terms = [anovaResults.factor.terms, 0];

%% Compute Type III sum of squares of each term
for nResults = 1:length(anovaResults.factor)
    all_terms = terms;
    current_terms = (all_terms == nResults);
    anovaResults.factor(nResults).beta_weights = beta_weights(current_terms);
    all_terms = find(current_terms == 0);
    df = anovaResults.factor(nResults).df;
    if isempty(all_terms)
        sum_sq = sum_sq_full_model;
    else
        %% Fit the model with the data using QR decomposition
        [temp_beta_weights, temp_line_fit] = icatb_constrained_lse(contrastMatrix(:, all_terms), data, constraint_matrix(:, all_terms));
        % Sum of squares of reduced model
        temp_sum_sq_full_model = sum(temp_line_fit.^2);
        sum_sq = sum_sq_full_model - temp_sum_sq_full_model;
        clear temp_fit_para;
    end
    anovaResults.factor(nResults).sum_sq = sum_sq;
    anovaResults.factor(nResults).mean_sq = sum_sq/df;
end

%% Error term
anovaResults.sum_sq_err = sum_sq_corrected_total - sum_sq_full_model;
anovaResults.erdf = length(data) - 1 - sum([anovaResults.factor.df]);
anovaResults.mean_sq_err  = anovaResults.sum_sq_err/anovaResults.erdf;

%% Compute F statistic and P - value for each term
for nResults = 1:length(anovaResults.factor)
    % F statistic
    Fstat = anovaResults.factor(nResults).mean_sq/anovaResults.mean_sq_err;
    % P value significance
    pval = 1 - icatb_spm_Fcdf(Fstat, anovaResults.factor(nResults).df, anovaResults.erdf);
    anovaResults.factor(nResults).Fstat = Fstat;
    anovaResults.factor(nResults).pval = pval;
end


%% Store contrast and constraint matrix
anovaResults.all_terms = terms;
anovaResults.contrastMatrix = contrastMatrix;
anovaResults.constraint_matrix = constraint_matrix;