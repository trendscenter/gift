function [V, D, CO] = icatb_calculate_em_pca(data, C, varargin)
%% Use Expectation maximization to compute PCA
% (http://www.cs.nyu.edu/~roweis/papers/empca.pdf)
%
% Inputs:
% 1. data - Voxels by Time
% 2. numpc - No. of comps
% 3. varargin - Arguments must be passed in pairs like maximum no. of
% iterations and tolerance.
%
% Outputs:
% 1. V - Eigen vectors
% 2. D - Eigen values diagonal matrix
% 3. CO - Final Transformation matrix
%

tol = 1.0e-4;
verbose = 1;
max_iter = 1000;

%% Loop over options
for i = 1:2:length(varargin)
    if strcmpi(varargin{i}, 'tolerance')
        tol = varargin{i + 1};
    elseif strcmpi(varargin{i}, 'max_iter')
        max_iter = varargin{i + 1};
    elseif strcmpi(varargin{i}, 'verbose')
        verbose = varargin{i + 1};
    end
end

%% No. of components and timepoints
numpc = size(C, 1);
tp = size(data, 2);

count = 0;

residual_err = 1;

%% Convert transformation matrix to unit norm
C = norm2(C);

if (verbose)
    disp('Using Expectation Maximization to do PCA');
end

usePinv = 0;
while ((residual_err > tol) && (count <= max_iter))
    
    count = count + 1;
    
    X = data*C'*pinv(C*C');
    
    X = double(X);
    
    if (count == 1)
        rank_initial = rank(X);
        if (rank_initial < size(C, 1))
            usePinv = 1;
            if (verbose)
                disp('Initial projection of transformation matrix on to the data is found to be rank deficient.');
                disp('Using pseudo-inverse to solve linear equations for rank deficient systems ...');
            end
        end
    end
    
    C_old = C;
    
    %% Use double precision when solving linear equations
    if (~usePinv)
        denom = X'*data;
        denom = double(denom);
        C = (X'*X) \ denom;
        clear denom;
    else
        C = pinv(X)*data;
    end
    
    clear X;
    
    C = norm2(C);
    
    residual_err = norm_resid(C, C_old);
    
    if (mod(count, 5) == 0)
        if (verbose)
            disp(['Step No: ', num2str(count), ' Norm of residual error: ', num2str(residual_err, '%0.6f')]);
        end
    end
    
end

clear C_old;

fprintf('\n');
if (verbose)
    if (residual_err <= tol)
        disp(['No of iterations required to converge is ', num2str(count)]);
    else
        disp(['Reached max iterations (', num2str(count), ')']);
        disp(['Residual error is ', num2str(residual_err)]);
    end
end

CO = C;

C = C';

%% Compute eigen vectors and values
[C, dd] = svd(C, 0);
clear dd;

cov_m = data*C;

cov_m = (cov_m'*cov_m)/(size(cov_m, 1) - 1 );

[V, D] = eig(cov_m, 'nobalance');

V = C*V;

if (verbose)
    disp('Done');
    fprintf('\n');
end


function C = norm2(C)
%% Normalize transformation matrix using slower way
%

for nC = 1:size(C, 1)
    C(nC, :) =  C(nC, :) ./ norm(C(nC, :), 2);
end

function residual_err = norm_resid(C, C_old)
%% Use norm2 of residual error
%

residual_err = 0;
for nC = 1:size(C, 1)
    res = C(nC, :) - C_old(nC, :);
    residual_err = residual_err + sum(res.^2);
end

residual_err = sqrt(residual_err);