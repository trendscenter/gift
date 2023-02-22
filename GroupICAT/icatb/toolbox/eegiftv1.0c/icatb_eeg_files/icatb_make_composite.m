function [funcImg, minICAIm, maxICAIm, minInterval, maxInterval] = icatb_make_composite(anat, funcData, ...
    imageValues, anatomicalView, structDIM, VOX)
% creates a composite image from structural and functional image
%
% Input:
% 1. anat - anatomical image
% 2. funcData
% 3. imageValues - 1 means positive and negative, 2 means positive, 3 means
% Absoulte, 4 means Negative
% 4. anatomicalView - Anatomical view
%
% Output:
% 1. funcImg - functional image converted to montage
% 2. minICAIm - Minimum value of functional image
% 3. maxICAIm - Max value of functional image
% 4. minInterval - Minimum interval
% 5. maxInterval - Max interval

% Initialise all vars
maxICAIm = zeros(1, size(funcData, 1));
minICAIm = zeros(1, size(funcData, 1));

%--get image in correct plane
if icatb_findstr(lower(anatomicalView), 'sagittal')
    % sagittal plane
    permuteOrder = [2 3 1];
elseif icatb_findstr(lower(anatomicalView), 'coronal')
    % coronal plane
    permuteOrder = [1 3 2];
else
    permuteOrder = [1 2 3];
end

is2D = (structDIM(3) == 1);

% Loop over number of functional images
for nImages = 1:size(funcData, 1)
    status = 0;
    % Get the current image
    func_im = funcData(nImages, :, :, :);
    
    % Overlay images
    [im, maxICAIm(nImages), minICAIm(nImages), minInterval, maxInterval] = icatb_overlayImages(func_im, anat, structDIM, structDIM, imageValues);
    
    % This will be required to plot the text of colorbar
    minICAIm(nImages) = round(10*minICAIm(nImages))/10;
    maxICAIm(nImages) = round(10*maxICAIm(nImages))/10;
    
    im = reshape(im, structDIM);
    
    % Permute the images
    im = permute(im, permuteOrder);
    DIM = structDIM(permuteOrder);
    if (~is2D)
        % Get the montage
        [montage_im] = icatb_returnMontage(im, [], DIM, VOX, 1);
    else
        montage_im = (fliplr(im))';
    end
    clear im;
    % store in the output cell array
    funcImg(nImages, :, :) = montage_im;
    clear montage_im;
end