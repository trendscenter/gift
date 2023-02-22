function di = nc_dist(X)
%% Euclidean distance between columns
%

pairs = nchoosek(1:size(X, 2), 2);

d = zeros(1, size(pairs, 1));
for np = 1:length(d)
    d(np) = norm(X(:, pairs(np, 1)) - X(:, pairs(np, 2)), 2);
end


di = zeros(size(X, 2), size(X, 2));
chk = tril(ones(size(X, 2), size(X, 2))) - eye(size(X, 2), size(X, 2));
di(chk == 1) = d;
di = di + di';