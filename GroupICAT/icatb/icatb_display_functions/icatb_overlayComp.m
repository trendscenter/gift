function [images, maxICAIM, minICAIM, maxInterval, minInterval] = icatb_overlayComp(icasig, structuralImage, ...
    imageValues, numComp, titleStr);


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