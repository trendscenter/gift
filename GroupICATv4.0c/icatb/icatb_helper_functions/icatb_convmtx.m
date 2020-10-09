function convoluted_matrix = icatb_convmtx(z, n)
% Convoluted matrix
% 
% Inputs:
%
% 1. z - Row vector or column vector
% 2. n - Length of another vector
%
% Output:
% convoluted_matrix - Convoluted matriz of size (z + n - 1) by n for column
% vector or n by (z + n - 1) for row vector

[numRows, numCols] = size(z);

% Return column vector
z = z(:);

size_new_z = length(z) + n - 1;

% convoluted matrix
convoluted_matrix = zeros(size_new_z, n);

% Initialise vars
startVec = 1;
endVec = length(z);

% Calculate convoluted matrix
for ii = 1:n
    convoluted_matrix(startVec:endVec, ii) = z;
    startVec = startVec + 1;
    endVec = endVec + 1;
end
% End for calculating convoluted matrix


if numRows < numCols
    convoluted_matrix = convoluted_matrix.';
end