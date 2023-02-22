function mat = icatb_vec2mat(vec, full)
%returns the correlation matrix from a vector
% if full == 1, returns the symmetric matrix
if nargin <2
    full = 0;
end

N = length(vec);
n = 1/2 + sqrt(1+8*N)/2;
mat = nan*ones(n);

temp = ones(n);
%% find the indices of the lower triangle of the matrix
IND = find((temp-triu(temp))>0);
mat(IND) = vec;

if full
    tempmat = flipud(rot90(mat));
    tempmat(IND) = vec;  
    mat = flipud(rot90(tempmat));
end