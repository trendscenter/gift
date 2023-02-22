function [overlayedImages, maxICAIM, minICAIM, minInterval, maxInterval] = icatb_check_overlayComp(icasig, ...
    structuralImage, imageValues, numComp)
% case 1: overlays components on the structural image
% case 2: overlays real comp on real part of structural image and similarly
% imaginary part on imaginary part of strcutural image
% case 3: overlays magnitude part of comp on magnitude part of structural
% and similarly phase part on phase part of structural

% check the data type of the component
if isreal(icasig)
    % real data type
    % real image
    titleStr = 'Computing Spatial Maps';
    % function to overlay images on the specified structural image
    [overlayedImages, maxICAIM, minICAIM, maxInterval, minInterval] = ...
        overlayComp(icasig, structuralImage, imageValues, numComp, titleStr);
elseif isa(icasig, 'complex_data')
    % handle complex_data class
    titleStr = 'Computing Spatial Maps';
    % function to overlay images on the specified structural image
    %[overlayedImages, maxICAIM, minICAIM, minInterval, maxInterval] = icatb_overlayComp(icasig, structuralImage, imageValues, numComp, titleStr);
    % get the field names
    fieldNames = fieldnames(icasig);
    
    % Initialise objects
    overlayedImages = complex_data;
    maxICAIM = overlayedImages; minICAIM = overlayedImages; minInterval = overlayedImages; 
    maxInterval = overlayedImages;
    
    % loop over fields
    for ii = 1:length(fieldNames)
        [temp1, temp2, temp3, temp4, temp5] = ...
            overlayComp(getfield(icasig, fieldNames{ii}), structuralImage, imageValues, numComp, titleStr);
        % set the field to the corresponding objects
        overlayedImages = setfield(overlayedImages, fieldNames{ii}, temp1);
        maxICAIM = setfield(maxICAIM, fieldNames{ii}, temp2);
        minICAIM = setfield(minICAIM, fieldNames{ii}, temp3);
        maxInterval = setfield(maxInterval, fieldNames{ii}, temp4);
        minInterval = setfield(minInterval, fieldNames{ii}, temp5);
    end
    % end for loop

else
    error('Unknown data type');
end
% end for checking the data type of the component

function [images, maxICAIM, minICAIM, maxInterval, minInterval] = overlayComp(icasig, structuralImage, ...
    imageValues, numComp, titleStr);
% sub-function to overlay component images

if ~exist('titleStr', 'var')
    titleStr = 'Computing Spatial Maps';
end

% - start status bar
StatusHandle = icatb_statusBar('init', 100, titleStr, '', '');
incrementAmount = 100/(numComp);

% get colormap
%cm = icatb_getColormap(1, imageValues, 1);

% -- Overlay images
structDIM = [size(structuralImage, 1), size(structuralImage, 2), size(structuralImage, 3)];
icaDIM = [size(icasig, 2), size(icasig, 3), size(icasig, 4)];

% Loop over components
for ii = 1:numComp
    [im, maxICAIM(ii), minICAIM(ii), minInterval, maxInterval] = icatb_overlayImages(icasig(ii, :, :, :), ...
        structuralImage, icaDIM, structDIM, imageValues);
    images(ii, :, :, :) = im;
    StatusHandle = icatb_statusBar('increment', incrementAmount);
end
clear icasig; clear im;
clear structuralImage;
icatb_statusBar('finished', StatusHandle);