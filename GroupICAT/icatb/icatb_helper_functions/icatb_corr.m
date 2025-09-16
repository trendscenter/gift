function C = icatb_corr(A, B)
%% Correlations

A = convert_to_col_vector(A);

if (~exist('B', 'var'))
    B = A;
else
    B = convert_to_col_vector(B);
end

A = A.*repmat(1./sqrt(sum(A.^2)), size(A, 1), 1);
B = B.*repmat(1./sqrt(sum(B.^2)), size(B, 1), 1);

C = A'*B;

function A = convert_to_col_vector(A)

if (length(A) == numel(A))
    A = A(:);
end

A = icatb_remove_mean(A);