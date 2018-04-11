function [a, R2, residual] = icatb_regress(y, X)
% Inputs are the observed data and the model matrix (one or more models)
% Assuming the rows are greater in number than columns.

% determine the rank of the X matrix (Assuming the observations are in
% columns)
rank_X = rank(X);

% for rank deficient systems pseudo-inverse will be used to determine the
% coefficients
if rank_X < size(X, 2)
    
    a = pinv(X)*y;
    
else
    
    % solve normal equations
    a = (X'*X) \ (X'*y); 
    
end

% calculating R-square statistic
if nargout > 1
    residual = y - X*a;
    R2 = 1 - (norm(residual, 2) / norm(y - mean(y), 2))^2;
end