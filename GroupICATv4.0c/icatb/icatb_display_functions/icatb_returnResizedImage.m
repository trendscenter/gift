function [structuralImage, compImages, coords, HInfo, text_left_right] = icatb_returnResizedImage(structFile, compFiles, ...
    anatomicalPlane, slices_in_mm, dataType, complexInfo, file_numbers, mask_file)
% return the structural image and the reshaped component images
% if the data type is complex then return complex images
%
% Inputs:
% 1. structFile - Full file path of the structural file
% 2. compFiles - Full file path of the component files
% Optional args:
% 3. anatomicalPlane - anatomical plane
% 4. slices_in_mm - Slices in mm (range)
% 5. dataType - data type ('real' or 'complex')
%
% Output:
% 1. structuralImage - Structural Image
% 2. compIMages - component Images
% 3. coords - Real world coordinates
% 4. HInfo - Header information
% 5. text_left_right - Text for plotting right and left text

icatb_defaults;
% global variable for flipping analyze images
global FLIP_ANALYZE_IMAGES;

% data type
if ~exist('dataType', 'var')
    dataType = 'real';
end

% anatomical plane
if ~exist('anatomicalPlane', 'var')
    anatomicalPlane = 'axial';
end

% check the file numbers
if ~exist('file_numbers', 'var')
    error('Component file numbers should be specified');
end

% make mask file empty if does not exist
if ~exist('mask_file', 'var')
    mask_file = [];
end

% display method
if ~exist('display_method', 'var')
    display_method = 'component';
end

% get the structural volume
%[structVol] = icatb_get_vol_nifti(structFile);

[structuralImage, structHInfo] = icatb_loadData(structFile);

structVol = structHInfo.V(1);


if ~exist('slices_in_mm', 'var')
    % get the slice definitions
    [sliceParameters] = icatb_get_slice_def(structVol, anatomicalPlane);
    % slices in mm
    slices_in_mm = sliceParameters.slices;
end

% if the data type is real
if strcmpi(dataType, 'real')

    % resize the image and return header Info
    [images, coords, HInfo, sliceinMm, text_left_right] = icatb_resizeImage(structVol, compFiles, anatomicalPlane, ...
        slices_in_mm, file_numbers);

    % get the structural image
    structuralImage = squeeze(images(1, :, :, :));

    % get the component images
    compImages = images(2:end, :, :, :);

else
    % search imaginary files and then find real files
    % or find magnitude and phase images

    % check only one file to determine real&imaginary or
    % magnitude&phase

    % File names of complex data
    P = icatb_get_complex_files_naming(compFiles, dataType, complexInfo, 'write');
    complexInfoRead = complexInfo.complexInfoRead;
    PFirst = str2mat(P.first); PSecond = str2mat(P.second);

    % resize the image and return header Info of first set of images
    [images, coords, HInfo, sliceinMm, First_text_left_right] = icatb_resizeImage(structVol, PFirst, anatomicalPlane, ...
        slices_in_mm, file_numbers);
    structuralImage = squeeze(images(1, :, :, :));
    firstImages = images(2:end, :, :, :);
    clear images;

    % resize the image and return header Info of second set of images
    [images, coords, HInfo, sliceinMm, Second_text_left_right] = icatb_resizeImage(structVol, PSecond, anatomicalPlane, ...
        slices_in_mm, file_numbers);
    structuralImage = squeeze(images(1, :, :, :));
    secondImages = images(2:end, :, :, :);
    clear images;

    compImages = struct('firstField', firstImages, 'secondField', secondImages);
    % form object
    compImages = complex_data(compImages);
    clear firstImages; clear secondImages;
    % loop over number of comps
    for ii = 1:size(Second_text_left_right, 1)
        text_left_right(2*ii - 1).tag = First_text_left_right(ii, :);
        text_left_right(2*ii).tag = Second_text_left_right(ii, :);
    end
    % text for left right flipping
    text_left_right = str2mat(text_left_right.tag);
end
% end for checking the data type

% if mask_image exists
if ~isempty(mask_file)
    %disp('Applying mask to the interpolated components ...');
    % resize mask image
    [mImages] = icatb_resizeImage(structVol, mask_file, anatomicalPlane, slices_in_mm, 1);
    % get the mask image
    maskIm = squeeze(mImages(2, :, :, :));
    clear mImages;
    maskvec = find(maskIm(:) > eps);
    clear maskIm;
    % get the masked out images
    compImages = getValues_in_mask(compImages, maskvec);
end
% end for checking mask image


function inputIm = getValues_in_mask(inputIm, maskvec)
% sub function to get the voxel values in the mask only

% size of the input image
size_input = size(inputIm);


% apply this only if maskvec is not empty
if ~isempty(maskvec)
    % loop over number of components
    for ii = 1:size_input(1)
        tempIm = squeeze(inputIm(ii, :, :, :));
        tempIm = tempIm(:);
        % initialise temp vector
        temp = zeros(1, length(tempIm));
        temp(maskvec) = tempIm(maskvec);
        clear tempIm;
        % set the new image to the current component
        inputIm(ii, :, :, :) = reshape(temp, [size_input(2:end)]);
        clear temp;
    end
    % end loop over number of components
end
% end for checking