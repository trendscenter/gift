function C = icatb_elem_by_elem_mult(A, B)
%% Function to do element by element multiplication of two matrices
% which maybe unequal in size.
%
% Inputs:
% 1. A - Matrix A of dimensions M by N
% 2. B - Matrix B of dimensions M by Q
%
% Outputs:
% C - Matrix of dimensions M by N*Q
%

if (size(A, 1) ~= size(B, 1))
    error('No. of rows of matrices A and B must be equal');
end

convertToLogical = 0;
if (islogical(A) && islogical(B))
    convertToLogical = 1;
end

% Initialise output matrix
C = zeros(size(A, 1), size(A, 2)*size(B, 2));

count = 0;
% Loop over columns of A
for nColA = 1:size(A, 2)
    % Loop over columns of B
    for nColB = 1:size(B, 2)
        count = count + 1;
        C(:, count) = A(:, nColA).*B(:, nColB);
    end
    % End of loop over columns of B
end
% End of loop over columns of A

if (convertToLogical)
    C = logical(C);
end