function M = icatb_low2Mat(vec, N, low_inds)
%% Construct matrix from lower triangular portion of matrix given the lower triangular indices
%

if (numel(vec) == length(vec))
    vec = vec(:)';
end

M = zeros(size(vec, 1), N, N);

for nF = 1:size(vec, 1)
    %tmp = vec(nF, :);
    tmp = zeros(N, N);
    tmp(low_inds) = vec(nF, :);
    tmp = tmp + tmp';
    % average diagonals
    tmp(eye(size(tmp))==1) = tmp(eye(size(tmp))==1)/2;
    M(nF, :, :) = tmp;
end