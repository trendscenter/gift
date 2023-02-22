function c = icatb_corrcov(X)
%% Compute correlation from covariance matrix
%
sigma = sqrt(diag(X));
sigma = sigma*sigma';
c = X./(sigma + eps);