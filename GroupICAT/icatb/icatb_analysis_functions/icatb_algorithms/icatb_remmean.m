function [newVectors, meanValue] = icatb_remmean(vectors);
%REMMEAN - remove the mean from vectors
%
% [newVectors, meanValue] = remmean(vectors);
%
% Removes the mean of row vectors.
% Returns the new vectors and the mean.
%
% This function is needed by FASTICA and FASTICAG

% 24.8.1998
% Hugo G�vert

newVectors = zeros (size (vectors));
meanValue = mean (vectors')';
newVectors = vectors - meanValue * ones (1,size (vectors, 2));
