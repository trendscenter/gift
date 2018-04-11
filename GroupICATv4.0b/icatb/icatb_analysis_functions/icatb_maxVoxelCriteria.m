function [maxValue, maxVoxelPos] = icatb_maxVoxelCriteria(tempMap, compMap)

% Computes the maximum value other than the zeros
% Sort IC based on maximum criteria. For every component 
% check the indices that have more than zero value
% get the maximum of the voxel value.


% Make the compMap to be column vector
if size(compMap, 1) == 1
    compMap = compMap';
end

% Make template column vector
if size(tempMap, 1) == 1
    tempMap = tempMap';
end

% Replace Nan with zeros
compMap(isnan(compMap)) = 0;
tempMap(isnan(tempMap)) = 0;

index = find(tempMap ~= 0);

if ~isempty(index)    
    % Get the maximum of all the voxels
    [maxValue, maxVoxelPos] = max(abs(compMap(index)));   
    maxVoxelPos = index(maxVoxelPos);
else
    error('Select a suitable mask.');
end
    