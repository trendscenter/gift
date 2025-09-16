function newMatrix = icatb_v_selcol (oldMatrix, maskVector);
% newMatrix = selcol(oldMatrix, maskVector);
%
% Selects the columns of the matrix that marked by one in the given vector.
% The maskVector is a column vector.
%
% This function is needed by PCAMAT

% 15.3.1998

if size (maskVector, 1) ~= size (oldMatrix, 2),
  error ('The mask vector and matrix are of uncompatible size.');
end

numTaken = 0;

% for i = 1 : size (maskVector, 1),
%   if maskVector (i, 1) == 1,
%     takingMask (1, numTaken + 1) = i;
%     numTaken = numTaken + 1;
%   end
% end

% take only the indices in the mask
takingMask = find(maskVector == 1);

numTaken = length(takingMask);

newMatrix = oldMatrix(:, takingMask);
