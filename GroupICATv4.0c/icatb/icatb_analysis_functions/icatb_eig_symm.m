function [V, D] = icatb_eig_symm(L, N, varargin)
%% Compute eigen values and vectors of a real symmetric matrix
%
% Inputs:
% 1. L
%      a. Symmetric matrix N by N for full storage
%                   OR
%      b. Lower triangular portion of symmetric matrix entered as a vector for packed storage.
%
% 2. N - Number of rows of symmetric matrix.
%
% 3. varargin - Arguments passed in pairs like 'num_eigs', 10,
% 'eig_solver', 'selective', 'create_copy', 1
%
%
% Outputs:
%
% 1. V - Eigen vectors in ascending order of eigen values
%
% 2. D - Eigen values diagonal matrix in ascending order of eigen values
%

%% Defaults
num_eigs = N;
eig_solver = 'selective';
create_copy = 1;
verbose = 1;

% Loop over input arguments
for i = 1:2:length(varargin)
    if strcmpi(varargin{i}, 'num_eigs')
        num_eigs = varargin{i + 1};
    elseif strcmpi(varargin{i}, 'eig_solver')
        eig_solver = varargin{i + 1};
    elseif strcmpi(varargin{i}, 'create_copy')
        create_copy = varargin{i + 1};
    elseif strcmpi(varargin{i}, 'verbose')
        verbose = varargin{i + 1};
    end
end
%% End of loop over input arguments

if (num_eigs > N)
    error('Error:EigenSolver', 'Number of eigen values (%d) desired is greater than the number of rows of matrix (%d) \n', num_eigs, N);
end

if (verbose)
    disp(['Covariance matrix size is ', num2str(N), ' ^2']);
    disp('Calculating eigendecomposition');
end


if (size(L, 1) == size(L, 2))
    %% Full storage
    
    if (size(L, 1) ~= N)
        error('Error:EigenSolver', 'Please check the dimensions of the real symmetric matrix as it doesn''t match the number of rows (%d)\n', N);
    end
    
    if (isa(L, 'double') && strcmpi(eig_solver, 'selective'))
        opts = struct('issym', 1, 'disp', 0, 'maxit', 1000);
        [V, D] = eigs(L, num_eigs, 'lm', opts);
    else
        % All eigen values
        [V, D] = eig(L, 'nobalance');
        sumAll = sum(diag(D));
    end
    
    if (verbose)
        disp('Sorting eigenvalues');
    end
    
    %% Sort in ascending order
    [D, inds] = sort(diag(D));
    
else
    %% Packed storage
    
    L = L(:);
    if (numel(L) ~= (N*(N + 1)/2))
        error('Number of elements of matrix must be equal to N*(N+1)/2');
    end
    
    if (strcmpi(eig_solver, 'selective'))
        % Selective eigen solver
        [V, D] = icatb_eig_symm_sel(L, N, num_eigs, create_copy);
    else
        % All eigen values
        [V, D] = icatb_eig_symm_all(L, N, create_copy);
        sumAll = sum(D);
    end
    
    if (verbose)
        disp('Sorting eigenvalues');
    end
    
    %% Sort in ascending order
    [D, inds] = sort(D);
    
end

if (verbose)
    disp('Selecting Desired Eigenvalues');
end

%% Return top few eigen values
inds = inds(end-num_eigs+1:end);
V = V(:, inds);
D = diag(D(end-num_eigs+1:end));

if (exist('sumAll', 'var') && verbose)
    sumUsed = sum(diag(D));
    retained = (sumUsed/sumAll)*100;
    fprintf('%g%% of (non-zero) eigenvalues retained.\n', retained);
end

% make sure eigen vectors have unit norm
V =  V*diag(1./sqrt(sum(V.^2)));