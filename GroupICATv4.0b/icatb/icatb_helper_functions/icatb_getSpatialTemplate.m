function [newSpatialTemplate, sizeSpatialTemplate] = icatb_getSpatialTemplate(spatialImage)

% Threshold image, convert to z scores

% Number of images
nImages = size(spatialImage, 1);

% Initialise vars
newSpatialTemplate = cell(nImages, 1);
spatialTemplate = cell(nImages, 1);
sizeSpatialTemplate = cell(nImages, 1);

% Storing all the templates in a cell array 
for ii = 1:nImages
    % Load spatial image
    spatialTemplate{ii} = icatb_loadData(spatialImage(ii, :));
end

% get the size of the array
for jj = 1:nImages
    sizeSpatialTemplate{jj} = size(spatialTemplate{jj});
end

% Store all the templates in a cell array
for ii = 1:nImages
    
    % reshape the array into single vector
    icasig = reshape(spatialTemplate{ii}, 1, prod(sizeSpatialTemplate{ii}));

    newSpatialTemplate{ii} = icasig;
    
    clear icasig;
end
