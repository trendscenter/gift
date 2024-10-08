function [im, numImagesX, numImagesY, slices_in_mm] = icatb_returnMontage(images, compIndex, DIM, VOX, slices_in_mm)
% returns the montage of the image
% 
% Input:
% 1. images - 
% 
% Ouput:
% 1. im - montage
% 2. numImagesX
% 3. numImagesY
% 4. slices_in_mm (order)

% --get montage of image
%a1b = zeros(size(images, 3), size(images, 2), size(images, 4), DIM(3));

numColorDims = 1;
if strcmpi(class(images), 'uint8')
    numColorDims = 3;
end

a1b = zeros([DIM(2), DIM(1), numColorDims, DIM(3)], class(images));
if ~isempty(compIndex)
    temp = reshape(images(compIndex, :, :, :, :), [DIM(1), DIM(2), numColorDims, DIM(3)]);
else
    temp = reshape(images, [DIM(1), DIM(2), numColorDims, DIM(3)]);
end
clear images;
% for k=1:DIM(3)
%     a1b(:,:,:,end-k+1) = (temp(:,end:-1:1,:,k)');
% end
for k=1:DIM(3)
    a1b(:,:,:,end-k+1) = permute(temp(:,end:-1:1,:,k), [2, 1, 3, 4]); %(images(:,end:-1:1,:,k)');
end
clear temp;
temp=reshape(a1b(:,:,:,:), [DIM(2), DIM(1), numColorDims, DIM(3)]);
slices_in_mm = slices_in_mm(end:-1:1);
%[im(:,:) numImagesX numImagesY] = icatb_get_montage(temp, structHInfo.VOX);
[im, numImagesX, numImagesY] = icatb_get_montage(temp, VOX);