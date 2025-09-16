function images = icatb_sliceRange(images, sliceRange)
% gets the slice range of images

if isstruct(images)
    % magnitude images
    magImages = images.mag;
    % phase images
    phaseImages = images.phase;
    
    clear images;
    % slice range of magnitude images
    magImages = magImages(:, :, :, sliceRange);
    
    % slice range of phase images
    phaseImages = phaseImages(:, :, :, sliceRange);
    
    % slice range of images
    images = struct('mag', magImages, 'phase', phaseImages);
    
else
    % slice range of images
    images = images(:, :, :, sliceRange);
end