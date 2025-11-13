function [X_new, W_whiten, W_dewhiten] = pca_whitening(X, N)
%%% PCA whitening a set of NxT matrices
% the resulting matrix has covariance equal T*I_n

[P, T, K] = size(X);

if nargin < 2
    N = P;
end

X_new = zeros(N, T, K);
W_whiten = zeros(N, P, K);
W_dewhiten = zeros(P, N, K);

for k = 1:K
    Xk = X(:, :, k);
    Xk = Xk - mean(Xk, 2); % remove mean is necessary
    [U, d] = eig(Xk*Xk'/(T - 1), 'vector'); % eig on X*X' is faster than svd on X
    [d, idx] = sort(d, 'descend');
    U = U(:, idx); % PxP

    U1 = U(:, 1:N); % PxN
    s1 = sqrt(d(1:N)); % Nx1
    W_whiten(:, :, k) = U1' ./ s1; % NxP
    W_dewhiten(:, :, k) = U1 .* s1'; % PxN
    X_new(:, :, k) = W_whiten(:, :, k) * Xk; % NxT
end

X_new = squeeze(X_new);
W_whiten = squeeze(W_whiten);
W_dewhiten = squeeze(W_dewhiten);

end