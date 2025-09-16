function constraint_matrix = icatb_anova_constraint(model, levels_vector)
%% Make Anova constraint matrix using the model and no. of levels in each
% factor.
%
% Inputs:
% 1. model - Model matrix
% 2. levels_vector - Number of levels in each factor
%
% Outputs:
% constraint_matrix - Constraint matrix
%

% Initialise contrast matrix
constraint_matrix = cell(size(model, 1), 1);

% Loop over model terms
for nM = 1:size(model, 1)
    currentTm = model(nM, :);
    % Get factors involved
    currentFactors = find(currentTm == 1);
    count = 0;
    for nF = 1:length(currentFactors)
        count = count + 1;
        cF = currentFactors(nF);
        num_levels = levels_vector(cF);
        % Construct the term name and constraints matrix
        if (count == 1)
            con = ones(1, num_levels);
        else
            % Interaction terms
            con = [kron(con, eye(num_levels));
                kron(eye(size(con, 2)), ones(1, num_levels))];
            %temp_cell = repmat({ones(1, num_levels)}, 1, num_levels);
            %con = [repmat(eye(num_levels, num_levels), 1, num_levels);
            %    blkdiag(temp_cell{:})];
        end
    end
    % End loop over factors involved
    constraint_matrix{nM} = con;
end
% End loop over model terms

% Generate block diagonal matrix (Add intercept term at the end)
constraint_matrix = blkdiag(constraint_matrix{:}, 0);