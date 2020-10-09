function matrixFormed = icatb_formMatrix(matrixToForm, DIM, maskVec, optional)
% Function to reshape the data using the maskVec as input 
% In old versions of GIFT pca all the voxels are stored. In the current
% version only the voxels that are in the mask.
%
% Inputs:
% 1. matrixToForm - matrix to form.
% 2. DIM - dimensions that the matrix needs to be transformed.
% 3. maskVec - mask indices
% 4. optional - 
%
% Output:
% matrixFormed - Reshaped Matrix

size_data = size(matrixToForm);

if ~exist('optional', 'var')
    optional = '2d';
end

if size(size_data, 1) ~= prod(DIM(1:3))
    % Output matrix
    matrixFormed = zeros(prod(DIM(1:3)), DIM(4));
    matrixFormed(maskVec, :) = matrixToForm;
    clear matrixToForm;
    if strcmpi(optional, '4d')
        matrixFormed = reshape(matrixFormed, DIM);
    end
else
    if strcmpi(optional, '4d')
        matrixFormed = reshape(matrixToForm, DIM);
        clear matrixToForm;
    else
        matrixFormed = matrixToForm;
        clear matrixForm;
    end
end





