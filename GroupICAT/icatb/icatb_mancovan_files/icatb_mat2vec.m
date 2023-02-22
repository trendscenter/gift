function [vec, IND] = icatb_mat2vec(mat)
% [vec, IND] = mat2vec(mat)
% Returns the lower triangle of mat
% mat should be square [m x m], or if 3-dims should be [n x m x m]

if ndims(mat) == 2
    
    [n,m] = size(mat);
    if n ~=m
        error('mat must be square!')
    end
elseif ndims(mat) == 3
    
    [n,m,m2] = size(mat);
    if m ~= m2
        error('2nd and 3rd dimensions must be equal!')
    end
end

temp = ones(m);
%% find the indices of the lower triangle of the matrix
IND = find((temp-triu(temp))>0);
if ndims(mat) == 2
    vec = mat(IND);
else
    mat = reshape(mat, n, m*m2);
    vec = mat(:,IND);
end

